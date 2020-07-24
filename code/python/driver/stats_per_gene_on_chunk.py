# Author: Karl Gemayel
# Created: 7/23/20, 9:36 AM
import copy
import logging
import argparse
import pandas as pd
from typing import *

# noinspection All
from Bio import SeqIO
from intervaltree import IntervalTree, Interval

import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import get_value, os_join
from mg_general.genome_splitter import GenomeSplitter
from mg_general.labels import Labels, Label, get_unique_gene_keys
from mg_general.labels_comparison_detailed import LabelsComparisonDetailed
from mg_general.shelf import compute_gc
from mg_io.labels import read_labels_from_file
from mg_options.parallelization import ParallelizationOptions

parser = argparse.ArgumentParser("Collects gene-level stats for runs with chunking.")

parser.add_argument('--pf-summary', required=True)
parser.add_argument('--pf-output', required=True)
parser.add_argument('--reference-tools', nargs="+")
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

def stats_per_gene_on_chunks_for_info(env, info, **kwargs):
    # type: (Environment, Dict[str, Any], Dict[str, Any]) -> pd.DataFrame

    gcfid = get_value(info, "Genome", required=True)
    pf_prediction = get_value(info, "Prediction", required=True)
    chunk_size_nt = get_value(info, "Chunk Size", required=True)
    tool = get_value(info, "Tool", required=True)

    key_value_delimiters_gff = {
        "mgm": " ",
        "mgm2": " ",
        "gms2": " ",
        "mprodigal": "=",
        "prodigal": "=",
    }

    labels = read_labels_from_file(
        pf_prediction, shift=0, ignore_partial=False,
        key_value_delimiter=key_value_delimiters_gff.get(tool.lower(), "=")
    )

def apply_genome_splitter_to_labels(seqname_to_info, labels):
    # type: (Dict[str, Iterable[[int, int]]], Labels) -> Labels

    # sketch: find which chunk(s) each label lies in
    # seqname_to_info = dict()
    # for x in genome_splitter.split_sequences_:
    #     seqname = x[0].id.split("_offset")[0]
    #     if seqname not in seqname_to_info:
    #         seqname_to_info = list()
    #     seqname_to_info[seqname].append(x)

    seqname_to_interval_tree = {
        s: IntervalTree(
            [Interval(x[0], x[1], (x[0], x[1])) for x in list_info if x[1] > x[0]]
        )
        for s, list_info in seqname_to_info.items()
    }

    list_labels = list()

    for lab in labels:
        left_chunk = seqname_to_interval_tree[lab.seqname()][lab.left()]
        right_chunk = seqname_to_interval_tree[lab.seqname()][lab.right()]

        if len(left_chunk) == 0:
            left_chunk = None
        else:
            left_chunk = left_chunk.pop().data

        if len(right_chunk) == 0:
            right_chunk = None
        else:
            right_chunk = right_chunk.pop().data


        if left_chunk == right_chunk:
            # no splitting
            list_labels.append(lab)
        else:
            # splitting, create two partial labels
            lab_left = copy.deepcopy(lab)
            lab_right = copy.deepcopy(lab)


            # label in left chunk has partial on right
            if left_chunk is not None:
                lab_left.coordinates().right = left_chunk[1]
                lab_left.coordinates().right -= (lab_left.length() % 3) # make multiple of 3
                if lab_left.get_attribute_value("partial") in {"10", "11"}:
                    lab_left.set_attribute_value("partial", "11")
                else:
                    lab_left.set_attribute_value("partial", "01")

                list_labels.append(lab_left)

            # label in right chunk has partial on left
            if right_chunk is not None:
                lab_right.coordinates().left = right_chunk[0]
                lab_right.coordinates().left += (lab_left.length() % 3) # make multiple of 3
                if lab_right.get_attribute_value("partial") in {"01", "11"}:
                    lab_right.set_attribute_value("partial", "11")
                else:
                    lab_right.set_attribute_value("partial", "10")

                list_labels.append(lab_right)

    return Labels(list_labels)

def stats_per_gene_on_chunks_for_genome(env, df_summary_genome, **kwargs):
    # type: (Environment, pd.DataFrame, Dict[str, Any]) -> pd.DataFrame
    """Input is summary dataframe for all entries of a given genome (including all chunks)
    Each row has the following:
    Genome: Label of genome
    Seqname: Sequence name/accession number
    Left-Chunk: left location of chunk in original sequence  (inclusive)
    Right-Chunk right location of chunk in original sequence (inclusive)
    Key: made from (possibly partial) 3' end, seqname name, strand
    For each tool:
        Left, right, partial-3', partial-5'

    """

    list_entries = list()
    for genome, df_genome in df_summary_genome.groupby("Genome", as_index=False):

        pf_sequence = os_join(env["pd-data"], genome, "sequence.fasta")
        sequences = SeqIO.to_dict(SeqIO.parse(pf_sequence, "fasta"))

        genome_gc = compute_gc(sequences)


        reference_tools = get_value(kwargs, "reference_tools", None)        # type: Union[List[str], None]

        # if multiple reference, merge them first by 5
        list_reference_labels = list()
        if reference_tools is not None:
            list_reference_labels = [read_labels_from_file(os_join(env["pd-runs"], genome, rt, "prediction.gff"),
                                                       ignore_partial=False, shift=-1)
                                 for rt in reference_tools]

        # merge labels by 5'
        merged_reference_labels = None
        if len(list_reference_labels) > 1:
            merged_reference_labels = LabelsComparisonDetailed(
                list_reference_labels[0], list_reference_labels[1]
            ).match_3p_5p("a")

            for i in range(2, len(list_reference_labels)):
                merged_reference_labels = LabelsComparisonDetailed(
                    merged_reference_labels, list_reference_labels[i]
                ).match_3p_5p("a")

        for chunk_size, df_chunk in df_genome.groupby("Chunk Size", as_index=False):        # type: int, pd.DataFrame

            # read all label files for chunk
            tool_to_labels = {
                df_chunk.at[idx, "Tool"]: read_labels_from_file(
                    df_chunk.at[idx, "Predictions"],
                    shift=-1, ignore_partial=False
                )
                for idx in df_chunk.index
            }

            # get all chunk intervals
            seqname_to_chunks = dict()
            for labels in tool_to_labels.values():
                for lab in labels:          # type: Label
                    seqname = lab.seqname()
                    if seqname not in seqname_to_chunks:
                        seqname_to_chunks[seqname] = set()

                    seqname_to_chunks[seqname].add(
                        (int(lab.get_attribute_value("chunk_left_in_original")),
                         int(lab.get_attribute_value("chunk_right_in_original")))
                    )

            # get all possible chunks, and map reference labels to it
            chunked_reference_labels = dict()
            if len(list_reference_labels) > 0:
                chunked_reference_labels = {
                    ref:
                    apply_genome_splitter_to_labels(seqname_to_chunks, labels) for ref, labels in zip(
                        reference_tools,
                        list_reference_labels
                    )
                }

            if merged_reference_labels is not None:
                chunked_reference_labels["=".join(reference_tools)] = apply_genome_splitter_to_labels(seqname_to_chunks,
                                                                                                      merged_reference_labels)


            # Add prediction info to dataframe

            name_to_labels = copy.copy(tool_to_labels)
            name_to_labels.update(chunked_reference_labels)

            all_labels = list(chunked_reference_labels.values()) + list(tool_to_labels.values())
            keys_3prime = get_unique_gene_keys(*all_labels)        # all 3prime keys


            # Each gene key will have a row in the dataframe
            # Columns will indicate whether it was 3p and 5p were predicted by each tool
            for key in keys_3prime:
                entry = dict()



                shortest_label = None
                for t in name_to_labels.keys():

                    label = name_to_labels[t].get_by_3prime_key(key)
                    if label is None:
                        entry[f"5p-{t}"] = None  # 5prime end
                    else:
                        entry[f"5p-{t}"] = label.get_5prime()
                        if shortest_label is None:
                            shortest_label = label
                        elif shortest_label.length() < label.length():
                            shortest_label = label

                        entry[f"Partial3p-{t}"] = label.incomplete_at_3prime()
                        entry[f"Partial5p-{t}"] = label.incomplete_at_5prime()



                # compute GC of label
                gene_gc = 0 # compute_gc(sequences, shortest_label)

                list_entries.append({
                    "3p-key": key,
                    "Genome": genome,
                    "Genome GC": genome_gc,
                    "Gene GC": gene_gc,
                    "Chunk Size": chunk_size,
                    **entry
                })

    return pd.DataFrame(list_entries)







def stats_per_gene_on_chunks(env, df_summary, pf_output, **kwargs):
    prl_options = get_value(kwargs, "prl_options", None)  # type: ParallelizationOptions

    # no parallelization
    if prl_options is None:
        df = stats_per_gene_on_chunks_for_genome(env, df_summary, **kwargs)
    else:
        pass
        # # PBS parallelization
        # if prl_options.safe_get("use-pbs"):
        #     pbs = PBS(env, prl_options, splitter=split_gil, merger=merge_identity)
        #     list_df = pbs.run(
        #         gil, helper_stats_per_gene,
        #         {
        #             "env": env, "tools": tools, **kwargs
        #         }
        #     )
        #     df = pd.concat(list_df, ignore_index=True, sort=False)
        #
        # # threading
        # else:
        #     list_df = run_n_per_thread(
        #         list(gil), stats_per_gene_for_gi,
        #         data_arg_name="gi",
        #         func_kwargs={
        #             "env": env, "tools": tools, **kwargs
        #         },
        #         simultaneous_runs=prl_options.safe_get("num-processors")
        #     )
        #
        #     df = pd.concat(list_df, ignore_index=True, sort=False)

    df.to_csv(pf_output, index=False)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = pd.read_csv(args.pf_summary)
    stats_per_gene_on_chunks(env, df, args.pf_output,
                             reference_tools=args.reference_tools)


if __name__ == "__main__":
    main(my_env, parsed_args)
