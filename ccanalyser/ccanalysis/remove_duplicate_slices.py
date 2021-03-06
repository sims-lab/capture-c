import argparse
import os
import sys
import pandas as pd
import dask
from dask.distributed import Client, LocalCluster
import dask.array as da
import dask.dataframe as dd
import random
from ccanalyser.utils.helpers import hash_column
import numpy as np


def main(
    input_files,
    deduplicated_fragments="deduplicated.hdf5",
    mode="fragments",
    output=None,
    shuffle=False,
    n_cores=8,
    max_memory="64GB",
    sample_name="",
    read_type="",
    stats_prefix="",
):

    if mode == "fragments":

        if shuffle:
            random.shuffle(input_files)

        read_csv_options = dict()

        if any(".gz" in fn for fn in input_files):
            read_csv_options["compression"] = "gzip"
            read_csv_options["blocksize"] = None

        deduplicated_ids = (
            dd.read_csv(
                input_files,
                sep="\t",
                usecols=["parent_read", "coordinates"],
                **read_csv_options,
            )
            .map_partitions(
                lambda df: df[["parent_read", "coordinates"]].apply(hash_column).astype(np.int64)
            )
            .drop_duplicates(subset="coordinates")
            .drop(columns='coordinates')
            .assign(duplicated=0)
            .set_index("parent_read")
            .to_hdf(deduplicated_fragments, key="deduplicated", mode="w")
        )

    elif mode == "slices":

        df_slices = (
            pd.read_csv(input_files, sep="\t")
            .assign(parent_read_hashed=lambda df: hash_column(df["parent_read"]))
            .set_index("parent_read_hashed")
        )

        n_fragments_total = df_slices["parent_read"].nunique()

        dd_fragments_deduplicated = dd.read_hdf(
            deduplicated_fragments, key="deduplicated", mode="r"
        )
        df_deduplicated = dd_fragments_deduplicated.join(
            df_slices, how="inner"
        ).compute()
        df_deduplicated.reset_index(drop=True).to_csv(output, sep="\t", index=False)

        n_fragments_unique = df_deduplicated["parent_read"].nunique()

        # Sort stats
        df_stats = pd.DataFrame()
        df_stats["stat_type"] = ["not-deduplicated", "deduplicated"]
        df_stats["stat"] = [n_fragments_total, n_fragments_unique]
        df_stats["sample"] = sample_name
        df_stats["read_type"] = read_type
        df_stats['read_number'] = 0
        df_stats["stage"] = "deduplicate_slices"
        df_stats.to_csv(f"{stats_prefix}.read.stats.csv", index=False)


# if __name__ == "__main__":


    # parser = argparse.ArgumentParser()
    # subparser = parser.add_subparsers(dest="mode")

    # parser_fragments = subparser.add_parser("fragments")
    # parser_fragments.add_argument("-i", "--input_files", nargs="+", required=True)
    # parser_fragments.add_argument("-f", "--deduplicated_fragments", required=True)
    # parser_fragments.add_argument(
    #     "--shuffle",
    #     help="shuffles the input files to randomise the deduplication",
    #     action="store_true",
    # )
    # parser_fragments.add_argument("-p", "--n_cores", default=8, type=int)
    # parser_fragments.add_argument("-m", "--max_memory", default="64GB", type=str)

    # parser_slices = subparser.add_parser("slices")
    # parser_slices.add_argument("-i", "--input_files", required=True)
    # parser_slices.add_argument("-f", "--deduplicated_fragments", required=True)
    # parser_slices.add_argument("-o", "--output", default="deduplicated.tsv.gz")
    # parser_slices.add_argument("-p", "--n_cores", default=8, type=int)
    # parser_slices.add_argument("-m", "--max_memory", default="64GB", type=str)
    # parser_slices.add_argument('-n', '--sample_name', help='sample name for stats')
    # parser_slices.add_argument('--read_type', help='read type for stats', default='pe')


    # args = parser.parse_args()

    #main(**vars(args))
