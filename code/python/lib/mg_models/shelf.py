# Author: Karl Gemayel
# Created: 2020-06-21, 6:27 p.m.

import logging
import numpy as np
import pandas as pd
from typing import *

from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Seq import Seq

from mg_container.msa import MSAType
from mg_general import Environment
from mg_general.general import get_value, os_join
from mg_io.general import remove_p
from mg_io.shelf import write_sequence_list_to_fasta_file

log = logging.getLogger(__name__)


def bin_by_gc(df, step=1):
    # type: (pd.DataFrame, int) -> List[Tuple[float, float, pd.DataFrame]]

    gc_ranges = range(20, 80, step)
    result = list()
    a = 0
    for b in gc_ranges:
        result.append(
            (a, b, df[(df["GC"] >= a) & (df["GC"] < b)])
        )
        a = b
    return result


def get_consensus_sequence(dict_mat):
    # type: (Dict[str, List[float]]) -> str

    num_positions = len(next(iter(dict_mat.values())))
    out = ""
    for n in range(num_positions):
        best_letter = None
        best_val = None

        for letter in dict_mat.keys():
            if best_letter is None:
                best_letter = letter
                best_val = dict_mat[letter][n]
            else:
                if dict_mat[letter][n] > best_val:
                    best_letter = letter
                    best_val = dict_mat[letter][n]

        out += best_letter

    return out


def get_position_distributions_by_shift(df, col, shifts):
    # type: (pd.DataFrame, str, List[int]) -> Dict[int, List[Dict[int,str]]]

    result = dict()
    for n in range(len(df.index)):
        idx = df.index[n]
        s = shifts[n]

        if s not in result:
            result[s] = list()

        result[s].append(df.at[idx, col])

    return result


def sort_sequences_by_first_non_gap_and_consensus(list_seqs):
    # type: (List[str]) -> List[str]

    def first_non_gap(l_seq):
        # type: (str) -> int
        for p in range(len(l_seq)):
            if l_seq[p] != "-":
                return p

        raise ValueError("Sequence is all gaps")

    pos_to_list_seqs = dict()
    for l in list_seqs:
        p = first_non_gap(l)
        if p not in pos_to_list_seqs.keys():
            pos_to_list_seqs[p] = list()

        pos_to_list_seqs[p].append(l)

    # reappend into single list and sort per position
    output = list()
    output_counts = list()
    for p in sorted(pos_to_list_seqs.keys()):
        # get counts per item
        counter = Counter(pos_to_list_seqs[p])
        sorted_counter = sorted([(x, counter[x]) for x in counter], key=lambda v: v[0])

        output += [x[0] for x in sorted_counter]
        output_counts += [x[1] for x in sorted_counter]

    return output, output_counts


def print_reduced_msa(msa_t, sort_by_starting_position=False, n=None):
    # type: (MSAType, bool) -> str

    list_sequences = [x.seq._data for x in msa_t.list_alignment_sequences]

    if sort_by_starting_position:
        list_sequences, counts = sort_sequences_by_first_non_gap_and_consensus(list_sequences)

    out = ""
    counter = 0
    for s, c in zip(list_sequences, counts):
        print(f"{s}\t{c}")
        out += "{}    {}\n".format(s, c)

        if n is not None and counter >= n:
            break
        counter += 1

    return out


def create_extended_numpy_for_column_and_shifts(df, col, update_shifts, new_width):
    # type: (pd.DataFrame, str, List[int], int) -> np.ndarray
    df = df[~df[col].isna()]  # we only need non-NA
    example = df.at[df.index[0], col]

    n = len(df)  # number of examples
    w = new_width
    b = len(example)  # number of bases (letters)

    mat = np.zeros((n, w, b), dtype=float)

    # fill the array
    for n_pos, idx in enumerate(df.index):
        dict_arr = df.at[idx, col]

        # for each base
        for b_pos, letter in enumerate(sorted(dict_arr.keys())):
            for w_pos, value in enumerate(dict_arr[letter]):
                shifted_w_pos = w_pos + update_shifts[n_pos]
                mat[n_pos, shifted_w_pos, b_pos] = value
    return mat


def run_msa_on_sequence_file(pf_fasta, pf_msa, **kwargs):
    # type: (str, str, Dict[str, Any]) -> None

    num_processors = get_value(kwargs, "num_processors", None)
    output_order = get_value(kwargs, "outputorder", "input-order")

    log.debug("Number of processors for MSA: {}".format(num_processors))
    other_options = dict()
    if num_processors is not None:
        other_options["threads"] = num_processors

    clustalw_cline = ClustalOmegaCommandline(
        "clustalo", infile=pf_fasta, outfile=pf_msa,
        outputorder=output_order,
        force=True,
        outfmt="clustal",
        **other_options
    )

    clustalw_cline()


def run_msa_on_sequences(env, sequences, **kwargs):
    # type: (Environment, List[Seq], Dict[str, Any]) -> MSAType

    pd_work = env["pd-work"]
    fn_tmp_prefix = get_value(kwargs, "fn_tmp_prefix", "", default_if_none=True)

    # write sequences to file
    pf_fasta = os_join(pd_work, "{}tmp_sequences.fasta".format(fn_tmp_prefix))
    remove_p(pf_fasta)
    write_sequence_list_to_fasta_file(sequences, pf_fasta)

    # run msa
    pf_msa = os_join(pd_work, "{}tmp_msa.txt".format(fn_tmp_prefix))
    run_msa_on_sequence_file(pf_fasta, pf_msa, **kwargs)

    msa_t = MSAType.init_from_file(pf_msa)

    remove_p(pf_msa, pf_fasta)

    return msa_t


def gather_consensus_sequences(env, df, col):
    # type: (Environment, pd.DataFrame, str) -> List[str]

    sequences = list()

    for idx in df.index:
        d = df.at[idx, col]     # type: Dict[str, List[float]]

        num_positions = len(next(iter(d.values())))
        out = ""
        for n in range(num_positions):
            best_letter = None
            best_val = None

            for letter in d.keys():
                if best_letter is None:
                    best_letter = letter
                    best_val = d[letter][n]
                else:
                    if d[letter][n] > best_val:
                        best_letter = letter
                        best_val = d[letter][n]

            out += best_letter
        sequences.append(out)

    return sequences

def create_numpy_for_column_with_extended_motif(env, df, col, other=dict()):
    # type: (Environment, pd.DataFrame, str) -> np.ndarray

    example = df.at[df.index[0], col]

    # run alignment
    consensus_seqs = gather_consensus_sequences(env, df, col)
    msa_t = run_msa_on_sequences(env, consensus_seqs,  outputorder="input-order")
    n = len(df)  # number of examples
    w = msa_t.alignment_length()
    b = len(example)  # number of bases (letters)
    other["msa_t"] = msa_t

    # get position of shift
    shifts = list()
    for s in msa_t.list_alignment_sequences:
        p = 0
        for pos in range(len(s)):
            if s[pos] != "-":
                p = pos
                break
        shifts.append(p)

    msa_t = run_msa_on_sequences(env, consensus_seqs, outputorder="tree-order")

    print_reduced_msa(msa_t, sort_by_starting_position=True)

    return create_extended_numpy_for_column_and_shifts(df, col, shifts, w), shifts
