# Author: Karl Gemayel
# Created: 8/26/20, 8:53 PM

import logging
import argparse
import pandas as pd
from typing import *

# noinspection All
from Bio import SeqIO

import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.genome_list import GenomeInfoList, GenomeInfo
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import os_join, get_value
from mg_general.genome_splitter import GenomeSplitter
from mg_io.general import mkdir_p, remove_p
from mg_io.labels import read_labels_from_file
from mg_io.shelf import read_sequences_for_gi, read_labels_for_gi
from mg_models.shelf import run_tool
from mg_options.parallelization import ParallelizationOptions
from mg_parallelization.generic_threading import run_n_per_thread
from mg_parallelization.pbs import PBS
from mg_pbs_data.mergers import merge_identity
from mg_pbs_data.splitters import split_gil
from mg_viz.shelf import mkstemp_closed

parser = argparse.ArgumentParser("DRIVER DESCRIPTION.")

parser.add_argument('--pf-gil', required=True)
parser.add_argument('--tool', required=True, choices=["mgm", "mgm2"])
parser.add_argument('--pf-output', required=True)

parser.add_argument('--dn-prefix', default="gct", help="Applies prefix to all run directories")
parser.add_argument('--force-split-in-intergenic', action='store_true')
parser.add_argument('--skip-if-exists', action='store_true')
parser.add_argument('--pf-parallelization-options')
parser.add_argument('--chunk-sizes-nt', nargs="+", default=[50000, 100000])# default=[250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 5000 ], type=int)


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


def get_features_from_prediction(tool, pf_prediction, gcode_true, tag):
    # type: (str, str, str, str) -> Dict[str, Dict[str, Any]]
    entries = dict()
    key_value_delimiters_gff = {
        "mgm": "=",
        "mgm2": " ",
        "gms2": " ",
        "mprodigal": "=",
        "prodigal": "=",
        "fgs": "=",
        "mga": "="
    }

    attribute_delimiter_gff = {
        "mgm": ","
    }

    # update labels file based on offset
    labels = read_labels_from_file(pf_prediction, shift=0, key_value_delimiter=key_value_delimiters_gff.get(
        tool.lower(), "="
    ), attribute_delimiter=attribute_delimiter_gff.get(tool.lower()), ignore_partial=False)

    labels_per_seqname = dict()
    for lab in labels:
        if lab.seqname() not in labels_per_seqname:
            labels_per_seqname[lab.seqname()] = list()

        labels_per_seqname[lab.seqname()].append(lab)

    counter = 0
    for seqname in labels_per_seqname:
        entries[seqname] = dict()

        total_score = 0
        avg_gene_length = 0
        avg_gc = 0
        num_genes = 0

        for lab in labels_per_seqname[seqname]:
            score = lab.get_attribute_value("score")
            try:
                score = float(score)
                total_score += score
            except (ValueError, TypeError):
                pass

            try:
                avg_gc += float(lab.get_attribute_value("gc"))
            except (ValueError, TypeError):
                pass

            num_genes += 1
            avg_gene_length += abs(lab.right() - lab.left() + 1)


        avg_gene_length /= num_genes if num_genes > 0 else 0
        avg_gc /= num_genes if num_genes > 0 else 0
        entries[seqname] = {
            f"{tag}: Total Score": total_score,
            f"{tag}: Average Gene Length": avg_gene_length,
            f"{tag}: Average Gene GC": avg_gc,
            f"{tag}: Number of Genes": num_genes
        }
        counter += 1
        # if counter > 5:
        #     break
    return entries


def build_gcode_features_for_sequence(env, tool, pf_sequences, pf_prediction, **kwargs):
    # type: (Environment, str, str, str, Dict[str, Any]) -> pd.DataFrame

    gcode_true = get_value(kwargs, "gcode_true", required=True, type=str)
    skip_if_exists = get_value(kwargs, "skip_if_exists", False)

    run_tool(env, pf_sequences, pf_prediction, tool, gcode=4, pf_mgm=None, fmt="ext", **kwargs)
    dict_entries_4 = get_features_from_prediction(tool, pf_prediction, gcode_true, tag="4")

    run_tool(env, pf_sequences, pf_prediction, tool, gcode=11, pf_mgm=None, fmt="ext", **kwargs)
    dict_entries_11 = get_features_from_prediction(tool, pf_prediction, gcode_true, tag="11")

    result = dict()
    for seqname in set(dict_entries_4.keys()).union(dict_entries_11.keys()):
        d4 = dict() if seqname not in dict_entries_4 else dict_entries_4[seqname]
        d11 = dict() if seqname not in dict_entries_11 else dict_entries_11[seqname]


        result[seqname] = d4
        result[seqname].update(d11)
        result[seqname]["True Gcode"] = gcode_true

    # sequence stats
    record_dict = SeqIO.to_dict(SeqIO.parse(pf_sequences, "fasta"))
    for seqname, r in record_dict.items():
        if seqname not in result:
            continue
        result[seqname]["Sequence Length"] = len(r)
        result[seqname]["4: Gene Density"] = result[seqname]["4: Number of Genes"] / len(r.seq) if "4: Number of Genes" in result[seqname] else 0
        result[seqname]["11: Gene Density"] = result[seqname]["11: Number of Genes"] / len(r.seq) if "11: Number of Genes" in result[seqname] else 0

    return pd.DataFrame(result.values())


def build_gcode_features_for_gi_for_chunk(env, gi, tool, chunk, **kwargs):
    # type: (Environment, GenomeInfo, str, int, Dict[str, Any]) -> pd.DataFrame

    dn_tool = get_value(kwargs, "dn_tool", tool)
    dn_prefix = get_value(kwargs, "dn_prefix", "")
    dn = tool

    gcode_true = gi.genetic_code

    # split genome into chunks
    gs = GenomeSplitter(
        read_sequences_for_gi(env, gi), chunk,
        labels=read_labels_for_gi(env, gi),
        allow_splits_in_cds=kwargs.get("allow_splits_in_cds")
    )

    pf_chunks = mkstemp_closed(dir=env["pd-work"], suffix=".fasta")
    gs.write_to_file(pf_chunks)

    list_entries = list()


    pd_run = os_join(env["pd-work"], gi.name, f"{dn_prefix}{dn}_{chunk}")
    mkdir_p(pd_run)

    pf_prediction = os_join(pd_run, "prediction.gff")
    results = build_gcode_features_for_sequence(env, tool, pf_chunks, pf_prediction,
                                                          gcode_true=gcode_true, **kwargs)

    results["Genome"] = gi.name
    list_entries.append(results)

    remove_p(pf_prediction)
    remove_p(pf_chunks)

    return pd.concat(list_entries, ignore_index=True, sort=False)


def build_gcode_features_for_gi(env, gi, tool, chunks, **kwargs):
    # type: (Environment, GenomeInfo, str, List[int], Dict[str, Any]) -> pd.DataFrame
    list_df = list()
    num_processors = get_value(kwargs, "num_processors", 1, valid_type=int)

    if num_processors > 1:
        list_df = run_n_per_thread(
            chunks, build_gcode_features_for_gi_for_chunk, "chunk",
            {
                "env": env, "gi": gi, "tool": tool, **kwargs
            }
        )

    else:
        list_df = list()
        for chunk in chunks:
            logger.debug(f"{gi.name};{chunk}")
            curr = build_gcode_features_for_gi_for_chunk(env, gi, tool, chunk, **kwargs)
            list_df.append(curr)

    return pd.concat(list_df, ignore_index=True, sort=False)


def build_gcode_features(env, gil, tool, chunks, **kwargs):
    # type: (Environment, GenomeInfoList, str, List[int], Dict[str, Any]) -> pd.DataFrame
    list_df = list()

    for gi in gil:
        list_df.append(
            build_gcode_features_for_gi(env, gi, tool, chunks, **kwargs)
        )

    return pd.concat(list_df, ignore_index=True, sort=False)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    gil = GenomeInfoList.init_from_file(args.pf_gil)
    prl_options = ParallelizationOptions.init_from_dict(env, args.pf_parallelization_options, vars(args))

    tool = args.tool
    dn_tool = tool

    chunks = args.chunk_sizes_nt

    if prl_options["use-pbs"]:
        pbs = PBS(env, prl_options, splitter=split_gil, merger=merge_identity)
        list_df = pbs.run(
            gil, build_gcode_features,
            {
                "env": env, "chunks": chunks, "tool": tool,
                "num_processors": prl_options.safe_get("pbs-ppn"),
                "allow_splits_in_cds": True,
                "dn_prefix": args.dn_prefix,
                "skip_if_exists": args.skip_if_exists
            }
        )
        df = pd.concat(list_df, ignore_index=True, sort=False)

    else:
        list_df = run_n_per_thread(
            list(gil), build_gcode_features_for_gi, "gi", {
                "env": env, "chunks": chunks,  "tool": tool, "dn_tool": dn_tool,
                "num_processors": 1,
                "allow_splits_in_cds": True,
                "dn_prefix": args.dn_prefix,
                "skip_if_exists": args.skip_if_exists
            }, simultaneous_runs=7
        )

        df = pd.concat(list_df, sort=False, ignore_index=True)

    df.to_csv(args.pf_output, index=False)


if __name__ == "__main__":
    main(my_env, parsed_args)
