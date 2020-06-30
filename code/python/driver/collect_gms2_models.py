# Author: karl
# Created: 2020-06-21, 8:29 a.m.

import logging
import argparse
import pandas as pd
from typing import *
from tqdm import tqdm

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_argparse.parallelization import add_parallelization_options
from mg_io.general import save_obj
from mg_general.general import os_join, get_value
from mg_parallelization.pbs import PBS
from mg_container.gms2_mod import GMS2Mod
from mg_pbs_data.splitters import split_gil
from mg_pbs_data.mergers import merge_dataframes
from mg_bio.general import compute_single_gc_from_file
from mg_general import Environment, add_env_args_to_parser
from mg_container.genome_list import GenomeInfoList, GenomeInfo
from mg_options.parallelization import ParallelizationOptions

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Collect GeneMarkS-2 model files from existing runs, and pickle them"
                                 "into single dataframe.")

parser.add_argument('--pf-gil', required=True, help="File containing genome information list")
parser.add_argument('--pf-output', required=True, help="Path to output file")
add_parallelization_options(parser)

add_env_args_to_parser(parser)
parsed_args = parser.parse_args()

# ------------------------------ #
#           Main Code            #
# ------------------------------ #

# Load environment variables
my_env = Environment.init_from_argparse(parsed_args)

# Setup logger
logging.basicConfig(level=parsed_args.loglevel)
logger = logging.getLogger("logger")  # type: logging.Logger


def collect_start_info_from_gi(env, gi):
    # type: (Environment, GenomeInfo) -> Dict[str, Any]
    pd_genome = os_join(env["pd-data"], gi.name)
    pf_sequence = os_join(pd_genome, "sequence.fasta")

    gc = gi.attributes.get("gc")
    if gc is None or float(gc) == 0:
        gc = compute_single_gc_from_file(pf_sequence)

    pd_genome_run = os_join(env["pd-runs"], gi.name)
    pd_gms2 = os_join(pd_genome_run, "gms2")
    pf_mod = os_join(pd_gms2, "GMS2.mod")

    mod = GMS2Mod.init_from_file(pf_mod)

    # clean up some things
    del mod.items["COD_MAT"]

    return {
        "Genome": gi.name,
        "GC": gc,
        "Mod": mod,
        # # **{
        # #     x: mod.items[x] for x in {
        # #         "GENOME_TYPE", "RBS_MAT", "RBS_MAT", "PROMOTER_MAT", "PROMOTER_WIDTH", "RBS_WIDTH",
        # #         "RBS_POS_DISTR", "PROMOTER_POS_DISTR", "ATG", "GTG", "TTG", "TAA", "TGA", "TAG",
        # #         "NON_MAT"
        # #     } if x in mod.items.keys()
        # }
    }


def collect_start_info_from_gil(env, gil, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, Any]) -> pd.DataFrame
    pf_output = get_value(kwargs, "pf_output", None, valid_type=str)

    list_entries = list()
    for gi in tqdm(gil, total=len(gil)):
        entry = collect_start_info_from_gi(env, gi)
        list_entries.append(entry)

    df = pd.DataFrame(list_entries)
    if pf_output is not None:
        df.to_csv(pf_output, index=False)

    return df

def collect_start_info_from_gil_and_print_to_file(env, gil, pf_output):
    # type: (Environment, GenomeInfoList, str) -> str

    list_entries = list()
    for gi in tqdm(gil, total=len(gil)):
        entry = collect_start_info_from_gi(env, gi)
        list_entries.append(entry)

    df = pd.DataFrame(list_entries)
    df.to_csv(pf_output, index=False)
    return pf_output




def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_gil)
    prl_options = ParallelizationOptions.init_from_dict(env, args.pf_parallelization_options, vars(args))

    if prl_options["use-pbs"]:
        pbs = PBS(env, prl_options,
                  splitter=split_gil,
                  merger=merge_dataframes
                  )

        df = pbs.run(
            data=gil,
            func=collect_start_info_from_gil,
            func_kwargs={
                "env": env,
            }
        )

    else:
        df = collect_start_info_from_gil(env, gil)

    save_obj(df, args.pf_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
