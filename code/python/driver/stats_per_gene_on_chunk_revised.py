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
from mg_argparse.parallelization import add_parallelization_options
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
from mg_parallelization.generic_threading import run_n_per_thread
from mg_parallelization.pbs import PBS
from mg_pbs_data.mergers import merge_identity, merge_dataframes_to_file
from mg_pbs_data.splitters import split_list

parser = argparse.ArgumentParser("Collects gene-level stats for runs with chunking.")

parser.add_argument('--pf-summary', required=True)
parser.add_argument('--pf-output', required=True)
parser.add_argument('--reference-tools-fn', nargs="+")
parser.add_argument('--reference-tools-fp', nargs="+")
parser.add_argument('--batch-size', default=None, type=int)
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

def convert_ref_labels_to_intervals(labels):
    # type: (Labels) -> Dict[str, Dict[str, IntervalTree]]

    result = dict()

    for lab in labels:
        accession = lab.seqname()
        strand = lab.strand()

        if accession not in result.keys():
            result[accession] = dict()
        if strand not in result[accession].keys():
            result[accession][strand] = IntervalTree()

        # create interval for label
        try:
            if lab.right() > lab.left():
                result[accession][strand].add(Interval(lab.left(), lab.right(), lab))
        except ValueError:
            pass


    return result



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
        overlapping = seqname_to_interval_tree[lab.seqname()].overlap(lab.left(), lab.right() + 1)

        left_chunk = seqname_to_interval_tree[lab.seqname()][lab.left()]
        right_chunk = seqname_to_interval_tree[lab.seqname()][lab.right()]

        if len(left_chunk) > 0:
            left_chunk = left_chunk.pop()
            overlapping.update([left_chunk])

        if len(right_chunk) > 0:
            right_chunk = right_chunk.pop()

            overlapping.update([right_chunk])




        if len(overlapping) > 2:
            for chunk in overlapping:

                # partial both
                if lab.left() < chunk[0] and lab.right() > chunk[1]:
                    lab_partial = copy.deepcopy(lab)
                    lab_partial.coordinates().right = chunk[1]
                    lab_partial.coordinates().left = chunk[0] + 1
                    lab_partial.coordinates().right -= (lab_partial.length() % 3)
                    lab_partial.set_attribute_value("partial", "11")
                    list_labels.append(lab_partial)

                # partial left
                elif lab.left() < chunk[0]:
                    lab_partial = copy.deepcopy(lab)
                    lab_partial.coordinates().left = chunk[0]
                    lab_partial.coordinates().left += (lab_partial.length() % 3) # make multiple of 3
                    if lab_partial.get_attribute_value("partial") in {"01", "11"}:
                        lab_partial.set_attribute_value("partial", "11")
                    else:
                        lab_partial.set_attribute_value("partial", "10")

                    list_labels.append(lab_partial)
                # partial right
                elif lab.right() > chunk[1]:
                    lab_partial = copy.deepcopy(lab)
                    lab_partial.coordinates().right = chunk[1]
                    lab_partial.coordinates().right -= (lab_partial.length() % 3)  # make multiple of 3
                    if lab_partial.get_attribute_value("partial") in {"10", "11"}:
                        lab_partial.set_attribute_value("partial", "11")
                    else:
                        lab_partial.set_attribute_value("partial", "01")

                    list_labels.append(lab_partial)




        if len(overlapping) <= 2:
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
                    assert (lab_left.length() % 3 == 0)
                    if lab_left.get_attribute_value("partial") in {"10", "11"}:
                        lab_left.set_attribute_value("partial", "11")
                    else:
                        lab_left.set_attribute_value("partial", "01")

                    list_labels.append(lab_left)

                # label in right chunk has partial on left
                if right_chunk is not None:
                    lab_right.coordinates().left = right_chunk[0]
                    lab_right.coordinates().left += (lab_right.length() % 3) # make multiple of 3
                    assert(lab_right.length() % 3 == 0)
                    if lab_right.get_attribute_value("partial") in {"01", "11"}:
                        lab_right.set_attribute_value("partial", "11")
                    else:
                        lab_right.set_attribute_value("partial", "10")

                    list_labels.append(lab_right)

    return Labels(list_labels)

def get_interval_with_largest_overlap(label, overlaps):
    # type: (Label, Set[Interval]) -> [Interval, float]

    best_interval = None
    score = None

    for o in overlaps:
        overlap_size = min(label.right(), o[1]) - max(label.left(), o[0])

        if best_interval is None or score < overlap_size:
            best_interval = o
            score = overlap_size

    return best_interval, score


def has_same_3prime_end(partial, full):
    # type: (Label, Label) -> bool

    if partial.strand() != full.strand():
        return False

    strand = partial.strand()
    buffer = 3      # used because sometimes we don't know exact frame

    if strand == "+":
        # partial at 3prime end
        if partial.incomplete_at_3prime():
            if full.right() >= partial.right()-buffer:
                return True
            else:
                return False
        elif full.incomplete_at_3prime():
            if partial.right() >= full.right()-buffer:
                return True
            else:
                return False
        else:
            # both complete
            return partial.right() == full.right()
    else:
        if partial.incomplete_at_3prime():
            return full.left() <= partial.left()+buffer
        elif full.incomplete_at_3prime():
            return partial.left() <= full.left()+buffer
        else:
            return partial.left() == full.left()

# def has_same_5prime_and_3prime_end(partial, full):
#     # type: (Label, Label) -> bool
#     if not has_same_3prime_end(partial, full):
#         return False
#
#     strand = partial.strand()
#
#




def compare_chunked_prediction_to_annotation(env, labels_pred, labels_ref_fp, labels_ref_fn):
    # type: (Environment, Labels, Labels) -> Dict[str, Any]


    pred_intervals = convert_ref_labels_to_intervals(labels_pred)

    # for each prediction, count it as false positive/negative, and find gene-start
    result = {
        "FP": 0,
        "FN": 0,
        "TP": 0,
        "Total Reference": 0
    }

    ref_intervals = convert_ref_labels_to_intervals(labels_ref_fp)
    for lab in labels_pred:

        # first, get interval where this label exists
        if lab.seqname() in ref_intervals and lab.strand() in ref_intervals[lab.seqname()]:
            overlaps = ref_intervals[lab.seqname()][lab.strand()].overlap(lab.left(), lab.right())
            if len(overlaps) == 0:
                result["FP"] += 1
            else:
                largest_overlap, score = get_interval_with_largest_overlap(lab, overlaps)

                if has_same_3prime_end(lab, largest_overlap.data):
                    result["TP"] += 1
                else:
                    result["FP"] += 1

        else:
            result["FP"] += 1

    ref_intervals = convert_ref_labels_to_intervals(labels_ref_fn)
    # compute false negatives
    for lab in labels_ref_fn:
        if lab.seqname() in pred_intervals and lab.strand() in pred_intervals[lab.seqname()]:
            overlaps = ref_intervals[lab.seqname()][lab.strand()].overlap(lab.left(), lab.right())

            if len(overlaps) == 0:
                result["FN"] += 1
            else:
                largest_overlap, score = get_interval_with_largest_overlap(lab, overlaps)

                if not has_same_3prime_end(lab, largest_overlap.data):
                    result["FN"] += 1
        else:
            result["FN"] += 1

    result["Total Reference"] = len(labels_ref_fn)

    return result



def merge_labels_by_5p(list_labels):
    # type: (List[Labels]) -> Labels
    # merge labels by 5'
    merged_reference_labels = None
    if len(list_labels) > 1:
        merged_reference_labels = LabelsComparisonDetailed(
            list_labels[0], list_labels[1]
        ).match_3p_5p("a")

        for i in range(2, len(list_labels)):
            merged_reference_labels = LabelsComparisonDetailed(
                merged_reference_labels, list_labels[i]
            ).match_3p_5p("a")
    else:
        merged_reference_labels = list_labels[0]

    return merged_reference_labels


def filter_labels_shorter_than(labels, threshold_nt, threshold_nt_partial=None):
    # type: (Labels, int, int) -> Labels

    if threshold_nt_partial is None:
        threshold_nt_partial = threshold_nt

    list_labels = list()
    for lab in labels:
        if lab.is_partial():
            if lab.length() >= threshold_nt_partial:
                list_labels.append(lab)
        else:
            if lab.length() >= threshold_nt:
                list_labels.append(lab)

    return Labels(list_labels, name=labels.name)



def stats_per_gene_on_chunks_for_genome(env, df_summary_genome, reference_tools_fn, reference_tools_fp, **kwargs):
    # type: (Environment, pd.DataFrame, List[str], List[str], Dict[str, Any]) -> pd.DataFrame
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

    reference_tools_fp = sorted(reference_tools_fp) if isinstance(reference_tools_fp, list) else reference_tools_fp
    reference_tools_fn = sorted(reference_tools_fn) if isinstance(reference_tools_fn, list) else reference_tools_fn

    list_entries = list()
    for genome, df_genome in df_summary_genome.groupby("Genome", as_index=False):

        pf_sequence = os_join(env["pd-data"], genome, "sequence.fasta")
        sequences = SeqIO.to_dict(SeqIO.parse(pf_sequence, "fasta"))
        clade = None
        if "Clade" in df_genome.columns:
            clade = df_genome.at[df_genome.index[0], "Clade"]

        genome_gc = compute_gc(sequences)

        list_reference_labels_fn = [
            read_labels_from_file(os_join(env["pd-runs"], genome, rt, "prediction.gff"),
                                  ignore_partial=False, shift=-1)
            for rt in reference_tools_fn
        ]

        list_reference_labels_fp = [
            read_labels_from_file(os_join(env["pd-runs"], genome, rt, "prediction.gff"),
                                  ignore_partial=False, shift=-1)
            for rt in reference_tools_fp
        ]

        ref_labels_fn = merge_labels_by_5p(list_reference_labels_fn)
        ref_labels_fp = merge_labels_by_5p(list_reference_labels_fp)

        ref_labels_fn = filter_labels_shorter_than(ref_labels_fn, 90)
        ref_labels_fp = filter_labels_shorter_than(ref_labels_fp, 90)


        for chunk_size, df_chunk in df_genome.groupby("Chunk Size", as_index=False):        # type: int, pd.DataFrame

            # read all label files for chunk
            tool_to_labels = {
                df_chunk.at[idx, "Tool"]: filter_labels_shorter_than(read_labels_from_file(
                    df_chunk.at[idx, "Predictions"],
                    shift=-1, ignore_partial=False
                ), 90)
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

            # # get all possible chunks, and map reference labels to it
            # chunked_reference_labels = dict()
            # if len(list_reference_labels) > 0:
            #     chunked_reference_labels = {
            #         ref:
            #         apply_genome_splitter_to_labels(seqname_to_chunks, labels) for ref, labels in zip(
            #             reference_tools,
            #             list_reference_labels
            #         )
            #     }
            #
            # chunked_reference_labels["=".join(reference_tools)] = apply_genome_splitter_to_labels(seqname_to_chunks,
            #                                                                                           merged_reference_labels)

            chunked_ref_labels_fn = apply_genome_splitter_to_labels(seqname_to_chunks, ref_labels_fn)
            chunked_ref_labels_fp = apply_genome_splitter_to_labels(seqname_to_chunks, ref_labels_fp)

            for t in tool_to_labels:
                list_entries.append(
                    {
                        "Genome": genome,
                        "Tool": t,
                        "Genome GC": genome_gc,
                        "Chunk Size": chunk_size,
                        **compare_chunked_prediction_to_annotation(env, tool_to_labels[t], chunked_ref_labels_fp,
                                                                   chunked_ref_labels_fn)
                    }
                )

            #
            # # Add prediction info to dataframe
            #
            # name_to_labels = copy.copy(tool_to_labels)
            # name_to_labels.update(chunked_reference_labels)
            #
            # all_labels = list(chunked_reference_labels.values()) + list(tool_to_labels.values())
            # keys_3prime = get_unique_gene_keys(*all_labels)        # all 3prime keys
            #




            # # Each gene key will have a row in the dataframe
            # # Columns will indicate whether it was 3p and 5p were predicted by each tool
            # for key in keys_3prime:
            #     entry = dict()
            #
            #
            #
            #     shortest_label = None
            #     for t in name_to_labels.keys():
            #
            #         label = name_to_labels[t].get_by_3prime_key(key)
            #         if label is None:
            #             entry[f"5p-{t}"] = None  # 5prime end
            #             entry[f"3p-{t}"] = None
            #         else:
            #             entry[f"5p-{t}"] = label.get_5prime()
            #             entry[f"3p-{t}"] = label.get_3prime()
            #             if shortest_label is None:
            #                 shortest_label = label
            #             elif shortest_label.length() < label.length():
            #                 shortest_label = label
            #
            #             entry[f"Partial3p-{t}"] = label.incomplete_at_3prime()
            #             entry[f"Partial5p-{t}"] = label.incomplete_at_5prime()
            #
            #
            #
            #     # compute GC of label
            #     gene_gc = 0 # compute_gc(sequences, shortest_label)
            #
            #     list_entries.append({
            #         "3p-key": key,
            #         "Genome": genome,
            #         "Genome GC": genome_gc,
            #         "Gene GC": gene_gc,
            #         "Chunk Size": chunk_size,
            #         "Runtime": df_chunk["Runtime"].mean(),
            #         "Clade": clade,
            #         **entry
            #     })


    return pd.DataFrame(list_entries)







def stats_per_gene_on_chunks(env, df_summary, pf_output, **kwargs):
    prl_options = get_value(kwargs, "prl_options", None)  # type: ParallelizationOptions
    append = get_value(kwargs, "append", False)
    mode = "w" if not append else "a"

    # no parallelization
    if prl_options is None:
        df = stats_per_gene_on_chunks_for_genome(env, df_summary, **kwargs)
    else:

        # stats need to be collected on from all entries of a genome. So make sure to pass
        # all entries of a given genome together
        list_df = [df_genome for _, df_genome in df_summary.groupby("Genome", as_index=False)]


        # PBS parallelization
        if prl_options.safe_get("use-pbs"):
            pbs = PBS(env, prl_options, splitter=split_list, merger=merge_dataframes_to_file)
            list_gen = pbs.run_on_generator(
                list_df, stats_per_gene_on_chunks_for_genome,
                {
                    "env": env, **kwargs
                },
                split_kwargs={
                    "arg_name_data": "df_summary_genome"
                },
                merge_kwargs={"pf_output": pf_output, "append": append}
            )
            df = None # = pd.concat([x for x in list_gen], ignore_index=True, sort=False)

        # threading
        else:
            list_df = run_n_per_thread(
                list_df, stats_per_gene_on_chunks_for_genome,
                data_arg_name="df_summary_genome",
                func_kwargs={
                    "env": env, **kwargs
                },
                simultaneous_runs=1#prl_options.safe_get("num-processors")
            )

            df = pd.concat(list_df, ignore_index=True, sort=False)

    if df is not None:
        df.sort_index(axis=1, inplace=True)
        df.to_csv(pf_output, index=False, mode=mode, header=mode=="w")


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = pd.read_csv(args.pf_summary)
    prl_options = ParallelizationOptions.init_from_dict(env, args.pf_parallelization_options,
                                                        vars(args))
    if args.batch_size is None:
        stats_per_gene_on_chunks(env, df, args.pf_output,
                                 reference_tools_fn=args.reference_tools_fn,
                                 reference_tools_fp=args.reference_tools_fp,
                                 prl_options=prl_options)
    else:
        df.reset_index(inplace=True)
        bs = args.batch_size

        list_df = list([x[1] for x in df.groupby("Genome", as_index=False)])

        start = 0
        end = min(bs, len(list_df))

        while start < len(list_df):
            curr_df = pd.concat(list_df[start:end], sort=False, ignore_index=True)

            stats_per_gene_on_chunks(env, curr_df, args.pf_output,
                                     reference_tools_fn=args.reference_tools_fn,
                                     reference_tools_fp=args.reference_tools_fp,
                                     append=start > 0,
                                     prl_options=prl_options)
            start = end
            end = min(start + bs, len(list_df))



if __name__ == "__main__":
    main(my_env, parsed_args)
