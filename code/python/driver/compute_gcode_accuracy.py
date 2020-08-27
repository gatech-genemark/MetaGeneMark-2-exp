# Author: Karl Gemayel
# Created: 8/25/20, 8:55 AM

import os
import logging
import argparse
import re

import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.genome_list import GenomeInfoList, GenomeInfo
from mg_general import Environment, add_env_args_to_parser
import mg_argparse.parallelization
from mg_general.general import get_value, os_join
from mg_general.genome_splitter import GenomeSplitter
from mg_io.general import mkdir_p
from mg_models.shelf import run_tool
from mg_options.parallelization import ParallelizationOptions
from mg_parallelization.generic_threading import run_n_per_thread
from mg_parallelization.pbs import PBS
from mg_pbs_data.mergers import merge_identity
from mg_pbs_data.splitters import split_gil
from mg_viz.shelf import mkstemp_closed
from mg_io.shelf import read_sequences_for_gi, read_labels_for_gi

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #


parser = argparse.ArgumentParser("DRIVER DESCRIPTION.")

parser.add_argument('--pf-gil', required=True)
parser.add_argument('--tools', required=True, nargs="+", choices=["mgm2", "mprodigal"], type=str.lower)

parser.add_argument('--dn_tools', nargs="+")
parser.add_argument('--dn-prefix', default=None, help="Applies prefix to all run directories")
parser.add_argument('--pf-summary', required=True, help="Output file that will contain summary of runs")
parser.add_argument('--force-split-in-intergenic', action='store_true')
parser.add_argument('--skip-if-exists', action='store_true')

parser.add_argument('--pf-mgm2-mod', type=os.path.abspath)
parser.add_argument('--pf-mgm-mod', type=os.path.abspath)
parser.add_argument('--chunk-sizes-nt', nargs="+",
                    default=[250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 5000], type=int)
                    # default=[250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 5000], type=int)

mg_argparse.parallelization.add_parallelization_options(parser)

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


def get_gcode_per_contig_for_mgm2(pf_prediction):
    # type: (str) -> Dict[str, str]

    pattern = re.compile(r"^# seqid:\s*([^\s ]).+genetic code:\s*(\d+)\s*.*$")
    gcode_per_contig = dict()

    contig = 1
    with open(pf_prediction, "r") as f:
        for line in f:
            line = line.strip()
            if "seqid" in line and "genetic code" in line:
                m = pattern.match(line)

                if m:
                    gcode = str(m.group(2))
                    gcode_per_contig[str(contig)] = gcode
                    contig += 1

    return gcode_per_contig


def get_gcode_per_contig_for_mprodigal(pf_prediction):
    # type: (str) -> Dict[str, str]
    pattern = re.compile(r"^# Model Data.+transl_table=\s*(\d+).*$")
    gcode_per_contig = dict()

    contig = 1
    with open(pf_prediction, "r") as f:
        for line in f:
            line = line.strip()
            if "Model Data" in line:
                m = pattern.match(line)

                if m:
                    gcode = str(m.group(1))
                    gcode_per_contig[str(contig)] = gcode
                    contig += 1

    return gcode_per_contig


def get_accuracy_gcode_predicted(tool, pf_prediction, gcode_true):
    # type: (str, str, str) -> Dict[str, Any]

    if tool == "mgm2":
        gcode_per_contig = get_gcode_per_contig_for_mgm2(pf_prediction)
    elif tool == "mprodigal" or tool == "prodigal":
        gcode_per_contig = get_gcode_per_contig_for_mprodigal(pf_prediction)
    else:
        raise ValueError("Unknown tool")

    num_matches = sum([1 for v in gcode_per_contig.values() if str(v) == gcode_true])
    num_mismatches = sum([1 for v in gcode_per_contig.values() if str(v) != gcode_true])

    total = len(gcode_per_contig)

    return {
        "Matches": num_matches,
        "Mismatches": num_mismatches,
        "Match Rate": num_matches / float(total) if total > 0 else 0,
        "Total Contigs": total,
        "Tool": tool
    }


def compute_gcode_accuracy_for_tool_on_sequence(env, tool, pf_sequences, pf_prediction, **kwargs):
    # type: (Environment, str, str, str, Dict[str, Any]) -> pd.Series

    gcode_true = get_value(kwargs, "gcode_true", required=True, type=str)
    skip_if_exists = get_value(kwargs, "skip_if_exists", False)

    if not skip_if_exists or (skip_if_exists and not os.path.isfile(pf_prediction)):
        run_tool(env, pf_sequences, pf_prediction, tool + "_autogcode", **kwargs)

    dict_entries = get_accuracy_gcode_predicted(tool, pf_prediction, gcode_true)
    dict_entries["True Gcode"] = gcode_true

    return pd.Series(dict_entries)


def compute_gcode_accuracy_for_tools_on_chunk(env, gi, tools, chunk, **kwargs):
    # type: (Environment, GenomeInfo, List[str], int, Dict[str, Any]) -> pd.DataFrame

    dn_tools = get_value(kwargs, "dn_tools", tools)
    dn_prefix = get_value(kwargs, "dn_prefix", "")

    gcode_true = int(gi.genetic_code)

    # split genome into chunks
    gs = GenomeSplitter(
        read_sequences_for_gi(env, gi), chunk,
        labels=read_labels_for_gi(env, gi),
        allow_splits_in_cds=kwargs.get("allow_splits_in_cds")
    )

    pf_chunks = mkstemp_closed(dir=env["pd-work"], suffix=".fasta")
    gs.write_to_file(pf_chunks)

    list_entries = list()

    for t, dn in zip(tools, dn_tools):
        pd_run = os_join(env["pd-work"], gi.name, f"{dn_prefix}{dn}_{chunk}")
        mkdir_p(pd_run)

        pf_prediction = os_join(pd_run, "prediction.gff")
        results = compute_gcode_accuracy_for_tool_on_sequence(env, t, pf_chunks, pf_prediction,
                                                              gcode_true=gcode_true, **kwargs)

        results["Genome"] = gi.name
        results["Chunk Size"] = chunk
        list_entries.append(results)

    return pd.DataFrame(list_entries)


def compute_gcode_accuracy_for_gi(env, gi, tools, chunks, **kwargs):
    # type: (Environment, GenomeInfo, List[str], List[int], Dict[str, Any]) -> pd.DataFrame
    list_df = list()
    num_processors = get_value(kwargs, "num_processors", 1, valid_type=int)

    if num_processors > 1:
        list_df = run_n_per_thread(
            chunks, compute_gcode_accuracy_for_tools_on_chunk, "chunk",
           {
               "env": env, "gi": gi, "tools": tools, **kwargs
           })


    else:

        list_df = list()
        for chunk in chunks:
            logger.debug(f"{gi.name};{chunk}")
            curr = compute_gcode_accuracy_for_tools_on_chunk(env, gi, tools, chunk, **kwargs)
            list_df.append(curr)

    return pd.concat(list_df, ignore_index=True, sort=False)


def compute_gcode_accuracy(env, gil, tools, chunks, **kwargs):
    # type: (Environment, GenomeInfoList, List[str], List[int], Dict[str, Any]) -> pd.DataFrame
    list_df = list()

    for gi in gil:
        list_df.append(
            compute_gcode_accuracy_for_gi(env, gi, tools, chunks, **kwargs)
        )

    return pd.concat(list_df, ignore_index=True, sort=False)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    gil = GenomeInfoList.init_from_file(args.pf_gil)
    prl_options = ParallelizationOptions.init_from_dict(env, args.pf_parallelization_options, vars(args))

    tools = args.tools
    chunks = args.chunk_sizes_nt
    dn_tools = args.dn_tools if args.dn_tools is not None else tools

    # check that both have the same length
    if len(tools) != len(dn_tools):
        raise ValueError(f"The 'tools' and 'dn-tools' arguments"
                         f" must have equal lengths: {len(tools)} != {len(dn_tools)}")

    if prl_options["use-pbs"]:
        pbs = PBS(env, prl_options, splitter=split_gil, merger=merge_identity)
        list_df = pbs.run(
            gil, compute_gcode_accuracy,
            {
                "env": env, "tools": tools, "chunks": chunks, "dn_tools": dn_tools,
                "pf_mgm2_mod": args.pf_mgm2_mod,
                "pf_mgm_mod": args.pf_mgm_mod,
                "num_processors": prl_options.safe_get("pbs-ppn"),
                "allow_splits_in_cds": not args.force_split_in_intergenic,
                "dn_prefix": args.dn_prefix,
                "skip_if_exists": args.skip_if_exists
            }
        )
        df = pd.concat(list_df, ignore_index=True, sort=False)

    else:
        list_df = run_n_per_thread(
            list(gil), compute_gcode_accuracy_for_gi, "gi", {
                "env": env, "chunks": chunks, "tools": tools, "dn_tools": dn_tools,
                "pf_mgm2_mod": args.pf_mgm2_mod,
                "pf_mgm_mod": args.pf_mgm_mod,
                "num_processors": 1,
                "allow_splits_in_cds": not args.force_split_in_intergenic,
                "dn_prefix": args.dn_prefix,
                "skip_if_exists": args.skip_if_exists
            }, simultaneous_runs=7
        )

        df = pd.concat(list_df, sort=False, ignore_index=True)

    df.to_csv(args.pf_summary, index=False)


if __name__ == "__main__":
    main(my_env, parsed_args)
