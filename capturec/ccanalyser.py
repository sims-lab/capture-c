import argparse
import os
import pysam
import pandas as pd
import time
from functools import wraps
import numpy as np
from datetime import timedelta

p = argparse.ArgumentParser()
p.add_argument('-i', '--input_bam', help='BAM file to parse')
p.add_argument('-a', '--annotations',
               help='Tab-delimited text file containing annotation for each read in bam file')
p.add_argument('--bed_output', help= 'bed output file prefix', default='bed')
p.add_argument('--stats_output', help= 'stats files output prefix', default='stats')
args = p.parse_args()


# assertions - check all input files exist
assert os.path.isfile(args.input_bam), "Input sam file not found"
assert os.path.isfile(args.annotations), "Annotation file  not found"


class SliceFilter():
    def __init__(self, slices):
        self.slices = slices

    @property
    def fragments(self):
        df =  (self.slices.sort_values(['parent_read', 'chrom', 'start'])
                           .groupby('parent_read', as_index=False)
                           .agg({'slice':'nunique',
                                 'pe': 'first',
                                 'mapped': 'sum',
                                 'multimapped':'sum',
                                 'capture':'nunique',
                                 'capture_count':'sum',
                                 'exclusion':'nunique',
                                 'exclusion_count': 'sum',
                                 'restriction_fragment': 'nunique',
                                 'blacklist': 'sum',
                                 'coordinates':'|'.join
                                }))
        df['capture'] = df['capture'] - 1
        df['exclusion'] = df['exclusion'] - 1
        df['reporter_count'] = np.where(df['capture_count'] > 0,
                                        df['mapped'] - (df['exclusion_count'] + df['capture_count'] + df['blacklist']),
                                        0)

        df = df.rename(columns={'capture': 'unique_capture_sites',
                                'exclusion': 'unique_exclusion_sites',
                                'restriction_fragment': 'unique_restriction_fragments',
                                'slice': 'unique_slices',
                                'blacklist': 'blacklisted_slices'
                                 }
                      )
        return df

    @property
    def slice_stats(self):
        return (self.slices.agg({'read_name': 'nunique',
                         'parent_read': 'nunique',
                         'mapped': 'sum',
                         'multimapped':'sum',
                         'capture': 'nunique',
                         'capture_count': lambda col: (col > 0).sum(),
                         'exclusion_count': lambda col: (col > 0).sum(),
                         'blacklist': 'sum'})
                   .rename({'read_name': 'unique_slices',
                            'parent_read': 'unique_fragments',
                            'multimapped': 'multimapping_slices',
                            'capture': 'unique_capture_sites',
                            'capture_count': 'number_of_capture_slices',
                            'exclusion_count': 'number_of_slices_in_exclusion_region',
                            'blacklist': 'number_of_slices_in_blacklisted_region',
                            }))
    @property
    def frag_stats(self):
        return (self.fragments.agg({'parent_read': 'nunique',
                        'mapped': lambda col: (col > 1).sum(),
                        'multimapped': lambda col: (col > 0).sum(),
                        'capture_count': lambda col: (col > 0).sum(),
                        'exclusion_count': lambda col: (col > 0).sum(),
                        'blacklisted_slices': lambda col: (col > 0).sum(),
                        'reporter_count': lambda col: (col > 0).sum()})
                    .rename({'parent_read': 'unique_fragments',
                             'multimapped': 'fragments_with_multimapping_slices',
                             'capture_count': 'fragments_with_capture_sites',
                             'exclusion_count': 'fragments_with_excluded_regions',
                             'blacklisted_slices': 'fragments_with_blacklisted_regions',
                             'reporter_count': 'fragments_with_reporter_slices'}))

    @property
    def reporters(self):
        return self.slices.query('capture == "-"')

    @property
    def captures(self):
        return self.slices.query('~(capture == "-")')


    def remove_unmapped_slices(self):
        self.slices = self.slices.query('mapped == 1')
        return self

    def remove_orphan_slices(self):
        '''Remove fragments with only one aligned slice'''
        fragments = self.fragments
        fragments_multislice = fragments.query('unique_slices > 1')
        self.slices = self.slices[self.slices['parent_read'].isin(fragments_multislice['parent_read'])]
        return self


    def remove_duplicate_re_frags(self):
        self.slices = (self.slices.sort_values('capture_count', ascending=False)
                                  .drop_duplicates(subset=['parent_read', 'restriction_fragment'], keep='first'))

        return self

    def remove_duplicate_slices(self):
        frags_deduplicated = self.fragments.drop_duplicates(subset="coordinates", keep='first')
        self.slices = self.slices[self.slices['parent_read'].isin(frags_deduplicated['parent_read'])]
        return self

    def remove_duplicate_slices_pe(self):
        fragments = self.fragments.assign(read_start=lambda df: df['coordinates'].str.split('|').str[0]
                                                                                 .str.split('-').str[0]
                                                                                 .str.split(':').str[1],
                                          read_end=lambda df: df['coordinates'].str.split('|').str[-1]
                                                                               .str.split('-').str[1])
        frags_duplicated = fragments.loc[(fragments.duplicated(subset=['read_start', 'read_end']))
                                         & (fragments['pe'] == 'pe')]
        self.slices = self.slices[~self.slices['parent_read'].isin(frags_duplicated['parent_read'])]
        return self

    def remove_exluded_and_blacklisted_slices(self):
        self.slices = self.slices.query('blacklist < 1 and exclusion_count < 1')
        return self

    def remove_non_reporter_fragments(self):
        frags_reporter = self.fragments.query('reporter_count > 0')
        self.slices = self.slices[self.slices['parent_read'].isin(frags_reporter['parent_read'])]
        return self

    def remove_multi_capture_fragments(self):
        frags_capture = self.fragments.query('0 < unique_capture_sites < 2')
        self.slices = self.slices[self.slices['parent_read'].isin(frags_capture['parent_read'])]
        return self

    @staticmethod
    def modify_re_frag(frag, adjust=1):
        enzyme, chrom, index = frag.split('_')
        return '_'.join([enzyme, chrom, str(int(index) + adjust)])


    def remove_multicapture_reporters(self, n_adjacent=1):
        '''Deals with an odd situation in which a reporter spanning two capture sites is not removed.
           This is due to the slice not being considered a capture or exclusion. To get around this
           the reporters on restriction fragments adjacent to capture sites are explicitly removed.

           The number of restriction sites to exclude can be modified with n_adjacent '''

        captures = self.captures
        re_frags = captures['restriction_fragment'].unique()

        # Generates a list of restriction fragments to be excluded from further analysis
        excluded_fragments = [modify_re_frag(frag, modifier)
                              for frag in re_frags
                              for modifier in range(-n_adjacent, n_adjacent + 1)]

        # Remove non-capture slices (reporters) in excluded regions
        self.slices = self.slices[(self.slices['capture_count'] > 0) |
                                  (~self.slices['restriction_fragment'].isin(excluded_fragments))]


        return self

def get_timing(task_name=None):
    '''Gets the time taken by the wrapped function'''
    def wrapper(f):
        @wraps(f)
        def wrapped(*args, **kwargs):
            time_start = time.perf_counter()
            result = f(*args, **kwargs)
            time_end = time.perf_counter()

            time_taken = timedelta(seconds = (time_end - time_start))
            print(f'Completed {task_name} in {time_taken} (hh:mm:ss.ms)')
            return result
        return wrapped
    return wrapper

def parse_alignment(aln):
    '''Parses reads from a bam file into a list'''

    read_name = aln.query_name
    parent_read, pe, slice_number = read_name.split("|")
    ref_name = aln.reference_name
    ref_start = aln.reference_start
    ref_end = aln.reference_end
    # Check if read mapped
    if aln.is_unmapped:
        mapped = 0
        multimapped = 0
        ref_name = "unmapped"
        ref_start = ""
        ref_end = ""
        coords = ""
    else:
        mapped = 1
        coords = f'{ref_name}:{ref_start}-{ref_end}'
        # Check if multimapped
        if aln.is_secondary:
            multimapped = 1
        else:
            multimapped = 0
    return [read_name, parent_read, pe, slice_number, mapped, multimapped,
            ref_name, ref_start, ref_end, coords]

@get_timing(task_name='processing BAM file')
def parse_bam(bam):
    '''Uses parse_alignment function to extract all data from the bam file and convert this to a dataframe'''

    df_bam = pd.DataFrame([parse_alignment(aln) for aln in pysam.AlignmentFile(bam, 'rb').fetch(until_eof=True)],
                          columns=['read_name','parent_read', 'pe', 'slice', 'mapped', 'multimapped',
                                   'chrom', 'start', 'end', 'coordinates'
                                   ]
                         )
    df_bam.set_index('read_name', inplace=True)
    return df_bam

@get_timing(task_name='merging annotations with BAM input')
def merge_annotations(df, annotations):
    '''Combines annotations with the parsed bam file output'''
    df_ann = pd.read_csv(annotations, sep='\t', header=0, index_col=0)
    return df.join(df_ann, how='inner')

@get_timing(task_name='filtering slices')
def filter_slices(df_slices):
    '''Performs filtering of slices using the SliceFilter class and outputs statitsics'''

    slice_filterer = SliceFilter(df_slices)
    stats_prefix = args.stats_output

    df_slice_stats = pd.DataFrame()

    # Unfiltered fragments
    slice_filterer = (slice_filterer.remove_unmapped_slices())
    df_slice_stats['mapped'] = slice_filterer.slice_stats

    # Remove fragments that do not have at least one unique capture site
    slice_filterer = (slice_filterer.remove_orphan_slices()
                                    .remove_multi_capture_fragments())
    df_slice_stats['contains_single_capture'] = slice_filterer.slice_stats

    # Filter reporters
    slice_filterer = (slice_filterer.remove_exluded_and_blacklisted_slices()
                                    .remove_non_reporter_fragments()
                                    .remove_multicapture_reporters()
                      )
    df_slice_stats['contains_capture_and_reporter'] = slice_filterer.slice_stats

    # Remove multiple occurences of the same restriction fragment and remove duplicate slices
    slice_filterer = (slice_filterer.remove_duplicate_re_frags()
                                    .remove_duplicate_slices()
                                    .remove_duplicate_slices_pe()
                      )
    df_slice_stats['duplicate_filtered'] = slice_filterer.slice_stats

    # Report statistics
    df_slice_stats.to_csv(f'{stats_prefix}.slice.stats', sep='\t')

    return slice_filterer

@get_timing(task_name='aggregating reporter slices by capture site and outputing .bed files')
def aggregate_by_capture_site(capture, reporter):

    capture = (capture[['parent_read', 'chrom', 'start', 'end', 'capture']]
                      .rename(columns={'chrom':'capture_chrom', 'start': 'capture_start', 'end': 'capture_end'})
                      .set_index('parent_read'))

    reporter = (reporter[['parent_read', 'read_name', 'chrom', 'start', 'end', 'restriction_fragment']]
                        .set_index('parent_read'))

    # Get stats for capture sites
    capture['capture'].value_counts().to_csv(f'{args.stats_output}.capture.stats')

    # Join reporters to captures using the parent read name
    captures_and_reporters = capture.join(reporter)

    # Determine cis/trans interaction statistics
    captures_and_reporters['cis/trans'] = np.where(captures_and_reporters['capture_chrom'] == captures_and_reporters['chrom'],
                                                   'cis', 'trans')
    interactions_by_capture = pd.DataFrame(captures_and_reporters.groupby('capture')['cis/trans']
                                                                 .value_counts())
    interactions_by_capture.to_csv(f'{args.stats_output}.cis_or_trans.reporter.stats')


    return captures_and_reporters.groupby('capture')


@get_timing(task_name='analysis of bam file')
def main():

    # Read bam file and merege annotations
    df_alignment = parse_bam(args.input_bam)
    df_alignment = merge_annotations(df_alignment, args.annotations).reset_index()

    #Filter slices using annotations
    filtered_slices = filter_slices(df_alignment)

    # Get capture and reporter slices from filtered slices
    df_capture_slices, df_reporter_slices = filtered_slices.captures, filtered_slices.reporters

    #Output combined capture and reporter fragments
    (df_capture_slices[['chrom', 'start', 'end', 'read_name']]
                      .to_csv(f'{args.bed_output}.capture.bed.gz', sep='\t', header=False, index=False))

    (df_reporter_slices[['chrom', 'start', 'end', 'read_name']]
                      .to_csv(f'{args.bed_output}.reporter.bed.gz', sep='\t', header=False, index=False))

    #Aggregate reporters by capture site and output these as bed files
    capture_site_aggregated_slices = aggregate_by_capture_site(df_capture_slices, df_reporter_slices)

    # Output the reporter DataFrame for each capture site
    for capture_site, df_rep in capture_site_aggregated_slices:
        df_rep.to_csv(f'{args.bed_output}_{capture_site}.tsv.gz',
                      sep='\t', index=False)





if __name__ == '__main__':
    main()
