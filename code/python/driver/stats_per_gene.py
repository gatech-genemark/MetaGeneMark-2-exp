# Author: Karl Gemayel
# Created: 6/29/20, 10:58 AM

import logging
import argparse
import os
from Bio import SeqIO
import pandas as pd
from typing import *

# noinspection All
import pathmagic


# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.genome_list import GenomeInfoList, GenomeInfo
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import get_value, os_join
from mg_general.labels import Labels, create_key_3prime_from_label
from mg_general.shelf import compute_gc
from mg_io.labels import read_labels_from_file
from mg_options.parallelization import ParallelizationOptions
from mg_parallelization.generic_threading import run_n_per_thread
from mg_parallelization.pbs import PBS
from mg_pbs_data.mergers import merge_identity
from mg_pbs_data.splitters import split_gil

parser = argparse.ArgumentParser("Collect comparisons of tools per gene.")

parser.add_argument('--pf-gil', required=True, help="File containing genome list.")
parser.add_argument('--tools', required=True, nargs="+", help="Name of tools")
parser.add_argument('--dn-tools', nargs="+", help="Names of run directories for the mentioned tools. If set, "
                                                  " should have a length equal to 'tools'.")
parser.add_argument('--pf-output', required=True, help="Output file")

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


def get_unique_gene_keys(*args):
    # type: (List[Labels]) -> Set[str]

    keys = set()
    for labels in args:
        keys = keys.union(set(create_key_3prime_from_label(lab) for lab in labels))

    return keys


def stats_per_gene_for_gi(env, gi, tools, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, str], Dict[str, Any]) -> pd.DataFrame

    pf_predictions = {
        t: os_join(env["pd-runs"], gi.name, tools[t], "prediction.gff") for t in tools.keys()
    }

    name_to_labels = {
        t: read_labels_from_file(pf_predictions[t], shift=0, name=t) for t in tools.keys()
    }  # type: Dict[str, Labels]

    keys_3prime = get_unique_gene_keys(*name_to_labels.values())

    pf_sequence = os_join(env["pd-data"], gi.name, "sequence.fasta")
    sequences = SeqIO.to_dict(SeqIO.parse(pf_sequence, "fasta"))

    genome_gc = compute_gc(sequences)

    # Each gene key will have a row in the dataframe
    # Columns will indicate whether it was 3p and 5p were predicted by each tool
    list_entries = list()
    for key in keys_3prime:
        entry = dict()

        shortest_label = None
        for t in tools.keys():

            label = name_to_labels[t].get_by_3prime_key(key)
            if label is None:
                entry[f"5p-{t}"] = None  # 5prime end
            else:
                entry[f"5p-{t}"] = label.get_5prime()
                if shortest_label is None:
                    shortest_label = label
                elif shortest_label.length() < label.length():
                    shortest_label = label

        # compute GC of label
        gene_gc = compute_gc(sequences, shortest_label)

        list_entries.append({
            "3p-key": key,
            "Genome": gi.name,
            "Genome GC": genome_gc,
            "Gene GC": gene_gc,
            **entry
        })

    df = pd.DataFrame(list_entries)
    # for t in tools.keys():
    #     df[f"5p-{t}"] = df[f"5p-{t}"].astype(int)

    return df


def helper_stats_per_gene(env, gil, tools, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, str], Dict[str, Any]) -> pd.DataFrame
    """

    :param env:
    :param gil:
    :param tools:
    :param kwargs:
        suppress_return: Default is True. If set to False, nothing is returned.
    :return: Either nothing or the dataframe, depending on the value of suppress_return.
    """

    list_df = list()
    for gi in gil:
        list_df.append(stats_per_gene_for_gi(env, gi, tools, **kwargs))

    if len(list_df) == 0:
        return pd.DataFrame()

    return pd.concat(list_df, ignore_index=True, sort=False)


def stats_per_gene(env, gil, tools, pf_output, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, str], str, Dict[str, Any]) -> None

    prl_options = get_value(kwargs, "prl_options", None)  # type: ParallelizationOptions

    # no parallelization
    if prl_options is None:
        df = helper_stats_per_gene(env, gil, tools, **kwargs)
    else:
        # PBS parallelization
        if prl_options.safe_get("use-pbs"):
            pbs = PBS(env, prl_options, splitter=split_gil, merger=merge_identity)
            raise NotImplementedError()
            # pbs.run(
            #     data={},
            #     , suppress_return = True
            # )

        # threading
        else:
            list_df = run_n_per_thread(
                list(gil), stats_per_gene_for_gi,
                data_arg_name="gi",
                func_kwargs={
                    "env": env, "tools": tools, **kwargs
                },
                simultaneous_runs=prl_options.safe_get("num-processors")
            )

            df = pd.concat(list_df, ignore_index=True, sort=False)

    df.to_csv(pf_output, index=False)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    gil = GenomeInfoList.init_from_file(os.path.abspath(args.pf_gil))

    tools = args.tools
    dn_tools = args.dn_tools if args.dn_tools is not None else tools

    # check that both have the same length
    if len(tools) != len(dn_tools):
        raise ValueError(f"The 'tools' and 'dn-tools' arguments"
                         f" must have equal lengths: {len(tools)} != {len(dn_tools)}")

    # collect tool name and directory together
    tool_to_dir = {a: b for a, b in zip(tools, dn_tools)}

    stats_per_gene(env, gil, tool_to_dir, args.pf_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
