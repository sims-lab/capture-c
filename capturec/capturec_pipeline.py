#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Capture-C Pipeline
==================

This script processes data from the capture-c or NG capture-c sequencing
protocols designed to identify 3d interactions in the genome.

It takes Illumina paired-end sequencing reads in fastq format (gz compressed)
as input and performs the following steps:

1: Quality control of raw reads (fastqc, multiqc)
2: Trimming of reads to remove adaptor sequence (trim-galore)
3: Combining overlapping read pairs (flash)
4: In silico digestion of reads (digest_fastq.py)
5: Classification of mapped read slices (ccanalyser.py)
6: Deduplication of fragment (dedup_fragments.py)

@authors: asmith, dsims
"""

# import packages
import sys
import os
import gzip
import seaborn as sns
import matplotlib.colors
from pybedtools import BedTool
import itertools
from pysam import FastxFile
from cgatcore import pipeline as P
from ruffus import mkdir, follows, transform, merge, originate, collate, split, regex, add_inputs, suffix

# Read in parameter file
P.get_parameters('capturec_pipeline.yml')
hub_dir = os.path.join(P.PARAMS["hub_publoc"], P.PARAMS['hub_name'])
assembly_dir = os.path.join(hub_dir, P.PARAMS['hub_genome'])


@follows(mkdir('ccanalyser'), mkdir('ccanalyser/restriction_enzyme_map/'))
@transform(P.PARAMS['genome_fasta'],
           regex(r'.*/(.*).fa.*'),
           r'ccanalyser/restriction_enzyme_map/\1.digest.bed.gz')
def digest_genome(infile, outfile):
    '''Digest genome using restriction enzyme and output fragments in bed file'''
    tmp = outfile.replace('.gz', '')
    statement = '''python %(run_options_scripts_dir)s/digest_genome.py
                   -i %(infile)s -o %(tmp)s -r %(ccanalyser_re)s
                   -l %(tmp)s.log
                   && gzip %(tmp)s'''
    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'])


@follows(mkdir('fastq_pre-processing'), mkdir('fastq_pre-processing/fastqc'))
@transform('*.fastq.gz',
           regex(r'(.*).fastq.gz'),
           r'fastq_pre-processing/fastqc/\1_fastqc.zip')
def qc_reads(infile, outfile):
    '''Quality control of raw sequencing reads'''
    outdir = os.path.dirname(outfile)
    statement = '''fastqc
                   -q
                   -t %(run_options_threads)s
                   --nogroup %(infile)s
                   --outdir %(outdir)s'''
    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'],
          job_threads=P.PARAMS['run_options_threads'])


@follows(mkdir('run_statistics'))
@merge(qc_reads, 'run_statistics/fastqc_report.html')
def multiqc_reads (infile, outfile):
    '''Collate fastqc reports into single report using multiqc'''

    bn = os.path.basename(outfile)
    dn = os.path.dirname(outfile)

    statement = '''export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc fastq_pre-processing/fastqc/
                   -o %(dn)s
                   -n %(bn)s'''
    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'],
          job_memory='2G')


@follows(mkdir('fastq_pre-processing/deduplicated'))
@collate('*.fastq.gz',
         regex(r'(.*)_[12].fastq.gz'),
         r'fastq_pre-processing/deduplicated/\1_1.fastq.gz')
def deduplicate_reads(infiles, outfile):

    '''Checks for duplicate read1/read2 pairs in a pair of fastq files
       any duplicates are discarded'''

    fq1, fq2 = infiles
    out1, out2 = outfile, outfile.replace('_1.fastq.gz', '_2.fastq.gz')
    logfile = out1.replace('_1.fastq.gz', '.log')

    if P.PARAMS['deduplication_pre-dedup']:

        statement = '''python %(run_options_scripts_dir)s/deduplicate_fastq.py
                               -1 %(fq1)s -2 %(fq2)s
                               --out1 %(out1)s --out2 %(out2)s
                               -l %(logfile)s
                               -c %(run_options_compression_level)s
                              '''
    else:
        # If deduplication turned off sylink input files and count number of reads
        statement = '''ln -s $(pwd)/%(fq1)s %(out1)s &&
                       ln -s $(pwd)/%(fq2)s %(out2)s &&
                       lc=$(zcat %(fq1)s | wc -l);
                       logfile=%(logfile)s;
                       echo -e "Read_pairs_processed\\t$(($lc / 4))\\n" > $logfile;
                       echo -e "Read_pairs_unique\\t$(($lc / 4))\\n" >> $logfile;
                       echo -e "Read_pairs_removed\\t0\\n" >> $logfile'''

    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'],
          job_memory='8G')

@follows(mkdir('fastq_pre-processing/trimmed'), deduplicate_reads)
@collate(r'fastq_pre-processing/deduplicated/*.fastq.gz',
         regex(r'fastq_pre-processing/deduplicated/(.*)_[12].fastq.gz'),
         r'fastq_pre-processing/trimmed/\1_1_val_1.fq.gz')
def trim_reads(infiles, outfile):
    '''Trim adaptor sequences using Trim-galore'''
    fastq1, fastq2 = infiles
    outdir = os.path.dirname(outfile)
    statement = '''trim_galore --cores %(run_options_threads)s %(trim_options)s -o %(outdir)s
                     %(fastq1)s %(fastq2)s'''
    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'],
          job_threads=P.PARAMS['run_options_threads'])


@follows(trim_reads, mkdir('fastq_pre-processing/flashed'))
@collate('fastq_pre-processing/trimmed/*.fq.gz',
         regex(r'fastq_pre-processing/trimmed/(.*)_[12]_.*.fq.gz'),
         r'fastq_pre-processing/flashed/\1.extendedFrags.fastq.gz')
def combine_reads(infiles, outfile):
    '''Combine overlapping paired-end reads using flash'''
    fastq1, fastq2 = infiles
    output_prefix = outfile.replace('.extendedFrags.fastq.gz', '')
    statement = '''flash
                   -z
                   -t %(run_options_threads)s
                   -o %(output_prefix)s
                    %(fastq1)s
                    %(fastq2)s'''
    P.run(statement,
      job_queue=P.PARAMS['run_options_queue'],
      job_threads=P.PARAMS['run_options_threads'])


@follows(mkdir('fastq_pre-processing/split'), combine_reads)
@transform('fastq_pre-processing/flashed/*.fastq.gz',
           regex(r'fastq_pre-processing/flashed/(.*).fastq.gz'),
           r'fastq_pre-processing/split/\1_0.fastq.gz')
def split_fastq(infile, outfile):
    '''Splits the combined (flashed) fastq files into chunks for parallel processing'''

    #Small error in function as only processes chunksize - 1 reads
    output_prefix = outfile.replace('_0.fastq.gz', '')
    statement = '''python
                   %(run_options_scripts_dir)s/split_fastq.py
                  -i %(infile)s
                  -n %(output_prefix)s
                  --chunk_size %(split_n_reads)s
                  -c %(run_options_compression_level)s '''
    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'])

@follows(mkdir('fastq_pre-processing/digested'), split_fastq)
@transform('fastq_pre-processing/split/*.fastq.gz',
           regex(r'fastq_pre-processing/split/(.*).extendedFrags_(\d+).fastq.gz'),
           r'fastq_pre-processing/digested/\1.flashed_\2.fastq.gz')
def digest_flashed_reads(infile, outfile):
    '''In silico restriction enzyme digest of combined (flashed) read pairs'''
    statement = '''python %(run_options_scripts_dir)s/digest_fastq.py
                   -o %(outfile)s
                   -l %(outfile)s.log
                   -r %(ccanalyser_re)s
                   -m 18
                   -c %(run_options_compression_level)s
                   -p 1
                   flashed
                   -i %(infile)s'''
    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'],
          job_threads=P.PARAMS['run_options_threads'])

@follows(split_fastq)
@collate('fastq_pre-processing/split/*.fastq.gz',
         regex(r'fastq_pre-processing/split/(.*).notCombined_[12]_(\d+).fastq.gz'),
         r'fastq_pre-processing/digested/\1.pe_\2.fastq.gz')
def digest_pe_reads(infiles, outfile):
    '''In silico restriction enzyme digest of non-combined (non-flashed) read pairs'''

    fq1, fq2 = infiles
    statement = '''python %(run_options_scripts_dir)s/digest_fastq.py
                   -l %(outfile)s.log
                   -r %(ccanalyser_re)s
                   -o %(outfile)s
                   -c %(run_options_compression_level)s
                   -m 18
                   -p 1
                   unflashed
                   -1 %(fq1)s
                   -2 %(fq2)s
                   '''

    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'],
          job_threads=P.PARAMS['run_options_threads'])

@follows(mkdir('aligned'))
@transform([digest_flashed_reads, digest_pe_reads],
           regex(r'fastq_pre-processing/digested/(.*).fastq.gz'),
           r'aligned/\1.bam')
def align_reads(infile, outfile):
    ''' Aligns digested fq files using bowtie2'''
    aligner = P.PARAMS['align_aligner']
    index_flag = P.PARAMS['align_index_flag'] if P.PARAMS['align_index_flag'] else ''
    options = P.PARAMS['align_options'] if P.PARAMS['align_options'] else ''

    statement = '''%(aligner)s %(options)s %(index_flag)s %(align_index)s %(infile)s
                    | samtools view -bS > %(outfile)s
                    && samtools sort %(outfile)s -o %(outfile)s.sorted.bam -m 2G -@ %(run_options_threads)s
                    && mv %(outfile)s.sorted.bam %(outfile)s'''
    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'],
          job_threads=P.PARAMS['run_options_threads'],
          job_memory='4G')

@collate(align_reads,
         regex(r'aligned/(.*)_(\d+).bam'),
         r'aligned/\1.bam')
def merge_bam_files(infiles, outfile):
    '''Combines bam files (by flashed/non-flashed status and sample)'''
    fnames = ' '.join(infiles)

    statement = '''samtools merge %(outfile)s %(fnames)s'''
    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'])

@follows(mkdir('aligned/mapping_statistics'))
@transform(merge_bam_files,
           regex(r'aligned/(.*).bam'),
           r'aligned/mapping_statistics/\1.picard.metrics')
def mapping_qc(infile, outfile):
    '''Uses picard CollectAlignmentSummaryMetrics to get mapping information.'''

    cmd = ['picard',
           'CollectAlignmentSummaryMetrics',
           'R=%(genome_fasta)s',
           'I=%(infile)s',
           'O=%(outfile)s',
           '&> %(outfile)s.log',
           ]

    statement = ' '.join(cmd)

    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'])


@merge(mapping_qc, 'run_statistics/mapping_report.html')
def mapping_multiqc(infiles, outfile):
    '''Combines mapping metrics using multiqc'''

    indir = os.path.dirname(infiles[0])
    out_fn = os.path.basename(outfile)
    out_dn = os.path.dirname(outfile)
    statement = '''export LC_ALL=en_US.UTF-8 &&
                   export LANG=en_US.UTF-8 &&
                   multiqc
                   %(indir)s
                   -o %(out_dn)s
                   -n %(out_fn)s'''
    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'],
          job_memory='8G')


@follows(mkdir('ccanalyser/annotations'))
@transform(align_reads,
           regex(r'aligned/(.*).bam'),
           r'ccanalyser/annotations/\1.bam.bed.gz')
def bam_to_bed(infile, outfile):
    '''Converts bam files to bed for faster intersection'''
    tmp = outfile.replace('.gz', '')
    statement =  '''bedtools bamtobed -i %(infile)s
                    | awk '{OFS="\\t"; print $1, $2, $3, $4"|"NR}'
                    | gzip > %(outfile)s'''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'])

@transform(P.PARAMS["ccanalyser_capture"],
           regex(P.PARAMS["ccanalyser_capture"]),
           r'ccanalyser/annotations/exclude.bed.gz')
def build_exclusion_bed(infile, outfile):
    '''Generates exclusion window around each capture site'''

    statement =  '''bedtools slop
                    -i %(infile)s -g %(genome_fai)s -b %(ccanalyser_exclude_window)s
                    | bedtools subtract -a - -b %(ccanalyser_capture)s
                    | gzip > %(outfile)s'''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'])

@transform(bam_to_bed,
           regex(r'ccanalyser/annotations/(.*).bam.bed.gz'),
           r'ccanalyser/annotations/\1.annotation.capture.count.gz')
def capture_intersect_count(infile, outfile):
    '''Report count of overlaps for the capture sites'''

    statement =  '''bedtools intersect -c -f 1
                    -a %(infile)s -b %(ccanalyser_capture)s
                    | awk 'BEGIN {OFS = "\\t"} {if ($7 != "0") {print $4, $7}}'
                    | sed '1i read_name\\tcapture_count'
                    | gzip > %(outfile)s'''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


@follows(build_exclusion_bed)
@transform(bam_to_bed,
           regex(r'ccanalyser/annotations/(.*).bam.bed.gz'),
           r'ccanalyser/annotations/\1.annotation.exclude.count.gz')
def exclusion_intersect_count(infile, outfile):
    '''Intersect reads with exclusion files.
    report count of overlaps for each input bed using -c'''

    statement =  '''bedtools intersect -c
                    -a %(infile)s -b ccanalyser/annotations/exclude.bed.gz
                    | awk 'BEGIN {OFS = "\\t"} {if ($7 != "0") {print $4, $7}}'
                    | sed '1i read_name\\texclusion_count'
                    | gzip > %(outfile)s'''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


@transform(bam_to_bed,
           regex(r'ccanalyser/annotations/(.*).bam.bed.gz'),
           r'ccanalyser/annotations/\1.annotation.capture.gz')
def capture_intersect(infile, outfile):
    '''Intersect reads with capture files.
    Capture slices must be fully contained within the capture restriction_fragment'''

    statement =  '''bedtools intersect -loj -f 1
                    -a %(infile)s -b %(ccanalyser_capture)s
                    | awk 'BEGIN {OFS = "\\t"} {if ($10 != ".") {print $4, $10}}'
                    | sed '1i read_name\\tcapture'
                    | gzip > %(outfile)s'''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


@follows(build_exclusion_bed)
@transform(bam_to_bed,
           regex(r'ccanalyser/annotations/(.*).bam.bed.gz'),
           r'ccanalyser/annotations/\1.annotation.exclude.gz')
def exclusion_intersect(infile, outfile):
    '''Intersect reads with exclusion files'''

    statement =  '''bedtools intersect -loj
                    -a %(infile)s -b ccanalyser/annotations/exclude.bed.gz
                    | awk 'BEGIN {OFS = "\\t"} {if ($10 != ".") {print $4, $10}}'
                    | sed '1i read_name\\texclusion'
                    | gzip > %(outfile)s'''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


@transform(bam_to_bed,
           regex(r'ccanalyser/annotations/(.*).bam.bed.gz'),
           add_inputs(digest_genome),
           r'ccanalyser/annotations/\1.annotation.re.gz')
def re_intersect(infiles, outfile):
    '''Intersect reads with restriction fragment map'''
    bam, genome = infiles
    statement =  '''bedtools intersect -loj
                    -a %(bam)s -b %(genome)s
                    | awk 'BEGIN {OFS = "\\t"} {if ($10 != ".") {print $4, $10}}'
                    | sed '1i read_name\\trestriction_fragment'
                    | gzip > %(outfile)s'''
    P.run(statement, job_queue=P.PARAMS['run_options_queue'])

@transform(bam_to_bed,
           regex(r'ccanalyser/annotations/(.*).bam.bed.gz'),
           r'ccanalyser/annotations/\1.annotation.blacklist.count.gz')
def blacklist_intersect_count(infile, outfile):
    '''Intersect reads with blacklisted regions.
    report count of overlaps for each input bed using -c '''

    if P.PARAMS['ccanalyser_blacklist']:
        statement =  '''bedtools intersect -c
                        -a %(infile)s -b %(ccanalyser_blacklist)s
                        | awk 'BEGIN {OFS = "\\t"} {if ($7 != "0") {print $4, $7}}'
                        | sed '1i read_name\\tblacklist'
                        | gzip > %(outfile)s'''
    else:
        statement = '''touch %(outfile)s'''

    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


@collate([capture_intersect_count, exclusion_intersect_count,
         capture_intersect, exclusion_intersect,
         re_intersect, blacklist_intersect_count],
    regex(r'ccanalyser/annotations/(.*).annotation.*'),
         r"ccanalyser/annotations/\1.annotations.tsv.gz")
def merge_annotations(infiles, outfile):
    '''Merge all intersections into a single file '''
    inlist = " ".join(infiles)
    statement = '''python %(run_options_scripts_dir)s/join_tsv.py
                   -f read_name
                   -o %(outfile)s
                   -i %(inlist)s'''
    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'])


@follows(merge_annotations,
         mkdir('ccanalyser/captures_and_reporters'),
         mkdir('ccanalyser/stats'))
@transform(align_reads,
           regex(r'aligned/(.*).bam'),
           add_inputs(r'ccanalyser/annotations/\1.annotations.tsv.gz'),
           r'ccanalyser/captures_and_reporters/\1.reporter.bed.gz')
def ccanalyser(infiles, outfile):
    '''Processes bam files and annotations, filteres slices and outputs
       reporter slices for each capture site'''

    bam, annotations = infiles
    bed_out = outfile.replace('.reporter.bed.gz', '')
    stats_out = bed_out.replace('captures_and_reporters', 'stats')
    statement =  '''python %(run_options_scripts_dir)s/ccanalyser.py
                    -i %(bam)s
                    -a %(annotations)s
                    --bed_output %(bed_out)s
                    --stats_out %(stats_out)s'''
    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'],
          job_memory=P.PARAMS['run_options_memory'])

@follows(ccanalyser, mkdir('ccanalyser/reporters_aggregated'))
@collate('ccanalyser/captures_and_reporters/*.tsv.gz',
         regex(r'ccanalyser/captures_and_reporters/(.*)\..*_\d+_(.*).tsv.gz'),
         r'ccanalyser/reporters_aggregated/\1_\2.tsv.gz')
def collate_ccanalyser_output(infiles, outfile):
    '''Combines multiple capture site bed files'''

    infiles = ' '.join(infiles)
    statement = '''cat %(infiles)s > %(outfile)s'''
    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'])

@follows(mkdir('ccanalyser/bedgraphs'))
@transform(collate_ccanalyser_output,
           regex(r'ccanalyser/reporters_aggregated/(.*).tsv.gz'),
           add_inputs(digest_genome),
           r'ccanalyser/bedgraphs/\1.bedgraph.gz')
def make_bedgraph(infiles, outfile):
    '''Intersect reporters with genome restriction fragments to create bedgraph'''
    tsv_fn = infiles[0]
    re_map = infiles[1]
    statement = '''python
                   %(run_options_scripts_dir)s/convert_tsv_to_bedgraph.py
                   --reporter_tsv %(tsv_fn)s
                   --re_map %(re_map)s
                   --output %(outfile)s'''

    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'])


@follows(mkdir('visualise'))
@transform(make_bedgraph,
           regex(r'ccanalyser/bedgraphs/(.*).bedgraph.gz'),
           r'visualise/\1.bigWig')
def make_bigwig(infile, outfile):
    '''Uses UCSC tools bedGraphToBigWig to generate bigWigs for each bedgraph'''

    tmp = infile.replace('.gz', '')
    statement = '''zcat %(infile)s > %(tmp)s
                   && bedGraphToBigWig %(tmp)s %(genome_chrom_sizes)s %(outfile)s
                   && rm %(tmp)s'''

    P.run(statement,
          job_queue=P.PARAMS['run_options_queue'])


@follows(mkdir(hub_dir))
@originate(os.path.join(hub_dir, 'hub.txt'))
def generate_hub_metadata(outfile):

    content = {'hub': P.PARAMS['hub_name'],
               'shortLabel': P.PARAMS['hub_short'] if P.PARAMS['hub_short'] else P.PARAMS['hub_name'],
               'longLabel': P.PARAMS['hub_long'] if P.PARAMS['hub_long'] else P.PARAMS['hub_name'],
               'genomesFile': 'genomes.txt',
               'email': P.PARAMS['hub_email'],
               'descriptionUrl': f'{P.PARAMS["hub_url"].rstrip("/")}/{P.PARAMS["hub_publoc"].strip("/")}',
               }

    with open(outfile, 'w') as w:
        for label, info in content.items():
            w.write(f'{label} {info}\n')


@follows(generate_hub_metadata)
@originate(os.path.join(hub_dir, 'genomes.txt'))
def generate_assembly_metadata(outfile):

    content = {'genome': P.PARAMS['hub_genome'],
               'trackDb': os.path.join(P.PARAMS['hub_genome'], 'trackDb.txt'),
              }

    with open(outfile, 'w') as w:
        for label, info in content.items():
            w.write(f'{label} {info}\n')


@follows(generate_hub_metadata, mkdir(assembly_dir))
@merge(make_bigwig, f'{assembly_dir}/trackDb.txt')
def generate_trackdb_metadata(infiles, outfile):
    def get_track_data(fn):
        return {'track': fn,
                'bigDataUrl': f'{P.PARAMS["hub_url"].rstrip("/")}/{(os.path.join(assembly_dir, fn)).lstrip("/")}',
                'shortLabel': fn,
                'longLabel': fn,
                'type': f'{fn.split(".")[-1]}',
                }

    # Generate all separate tracks
    bigwig_tracks_all = [get_track_data(os.path.basename(fn)) for fn in infiles]

    # Add colours to tracks
    if not P.PARAMS['hub_colors']:
        colors = sns.color_palette('husl', len(bigwig_tracks_all))
        for track, color in zip(bigwig_tracks_all, colors):
            track['color'] = ','.join([str(c * 255) for c in color])
    else:
        for track, color in zip(bigwig_tracks_all,
                                itertools.cycle(P.PARAMS['hub_colors'].split(' '))):

            track['color'] = ','.join([str(c * 255)
                                      for c in matplotlib.colors.to_rgb(color)])


    # Write track data separated
    with open(outfile, 'w') as w:
        for track in bigwig_tracks_all:
            for label, data in track.items():
                w.write(f'{label} {data}\n')
            # Need to separate each track with a new line
            w.write('\n')

        #Group tracks by sample name and make separate combined tracks for each
        sample_key = lambda d: d['track'].split('.')[0].split('_')[0]
        bigwig_tracks_grouped = {sample: list(track) for sample, track in
                                 itertools.groupby(sorted(bigwig_tracks_all, key=sample_key), key=sample_key)}

        for sample, grouped_tracks in bigwig_tracks_grouped.items():

            # Generate overlay track
            combined_track_details = {'track': f'{sample}_combined',
                                      'container': 'multiWig',
                                      'aggregate': 'transparentOverlay',
                                      'showSubtrackColorOnUi': 'on',
                                      'type': 'bigWig 0 250',
                                      'shortLabel': f'{sample}_combined',
                                      'longLabel': f'{sample}_all_capture_probes_combined',
                                     }

            # Write overlay track
            for label, data in combined_track_details.items():
                w.write(f'{label} {data}\n')
            w.write('\n')


            # Write sub-tracks
            for track in grouped_tracks:
                track['track'] = track['track'].replace('.bigWig', '_subtrack.bigWig')
                for label, data in track.items():
                    w.write(f'\t{label} {data}\n')
                w.write(f'\tparent {combined_track_details["track"]}\n')
                # Need to separate each track with a new line
                w.write('\n')


@follows(generate_trackdb_metadata)
@transform(make_bigwig,
          regex('visualise/(.*).bigWig'),
          f'{assembly_dir}' + r'/\1.bigWig')
def link_bigwigs(infile, outfile):
    try:
        infile_fp = os.path.abspath(infile)
        os.symlink(infile_fp, outfile)
    except Exception as e:
        print(e)

@follows(ccanalyser)
@merge(['fastq_pre-processing/deduplicated/*.log',
        'fastq_pre-processing/digested/*.log',
        'ccanalyser/stats/*.slice.stats',
        'ccanalyser/stats/*.reporter.stats'
        ],
        'run_statistics/combined_stats.tsv')
def aggregate_stats(infiles, outfile):

    dedup = ' '.join([fn for fn in infiles if 'deduplicated/' in fn])
    digestion = ' '.join([fn for fn in infiles if 'digested/' in fn])
    slices = ' '.join([fn for fn in infiles if '.slice.stats' in fn])
    reporters = ' '.join([fn for fn in infiles if '.reporter.stats' in fn])
    outdir = os.path.dirname(outfile)

    statement =  '''python %(run_options_scripts_dir)s/aggregate_statistics.py
                    --deduplication_stats %(dedup)s
                    --digestion_stats %(digestion)s
                    --ccanalyser_stats %(slices)s
                    --reporter_stats %(reporters)s
                    --output_dir %(outdir)s
                '''

    P.run(statement, job_queue=P.PARAMS['run_options_queue'])

@follows(aggregate_stats)
@merge('run_statistics/*.tsv', r'run_statistics/visualise_run_statistics.html')
def build_report(infile, outfile):
    '''Run jupyter notebook for reporting and plotting. First moves the notebook
       then converts to html'''
    statement = '''papermill
                   %(run_options_scripts_dir)s/visualise_capture-c_stats.ipynb
                   run_statistics/visualise_run_statistics.ipynb
                   -p directory $(pwd)/run_statistics/
                   &&
                   jupyter
                   nbconvert
                   --no-input
                   run_statistics/visualise_run_statistics.ipynb
                   run_statistics/visualise_run_statistics.html
                   '''

    P.run(statement, job_queue=P.PARAMS['run_options_queue'])


if __name__ == "__main__":
        sys.exit( P.main(sys.argv))
