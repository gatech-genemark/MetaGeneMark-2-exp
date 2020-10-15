# Author: Karl Gemayel
# Created: 10/7/20, 3:45 PM

import logging
import argparse
import numpy as np
import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import next_name
from mg_io.general import load_obj
from mg_models.shelf import create_numpy_for_column_with_extended_motif, fix_genome_type, \
    get_consensus_sequence, helper_clusters_by_heuristic

parser = argparse.ArgumentParser("DRIVER DESCRIPTION.")

parser.add_argument('--pf-data', required=True)

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


def clusters_by_msa(env, df):
    # type: (Environment, pd.DataFrame) -> np.ndarray
    _, update_shifts = create_numpy_for_column_with_extended_motif(env, df, "RBS_MAT")

    return update_shifts


def clusters_by_heuristic(env, df):
    # type: (Environment, pd.DataFrame) -> np.ndarray
    return helper_clusters_by_heuristic(env, df)


def collect_sequences_in_groups(df, col):
    # type: (pd.DataFrame, str) -> Dict[int, List[str]]

    groups = {c: set() for c in df[col].unique()}       # type: Dict[int, Set]

    for idx in df.index:
        groups[df.loc[idx, col]].add(df.loc[idx, "CONSENSUS_RBS_MAT"])

    groups_sorted = dict()
    for g in groups:
        groups_sorted[g] = sorted(groups[g])

    return groups_sorted






def compare_clustering_algorithms_for_gc_and_group(env, df):
    # type: (Environment, pd.DataFrame) -> Dict[str, Any]

    msa_clusters = clusters_by_msa(env, df)
    heu_clusters = clusters_by_heuristic(env, df)

    df["MSA"] = msa_clusters
    df["Heuristic"] = heu_clusters

    groups_msa = collect_sequences_in_groups(df, "MSA")
    groups_heu = collect_sequences_in_groups(df, "Heuristic")

    print(groups_heu)
    print(groups_msa)

def compare_clustering_algorithms(env, df):
    # type: (Environment, pd.DataFrame) -> None

    low_gc = 0
    col = "RBS_MAT"

    for high_gc in range(30, 70, 5):

        df_gc_subset = df[(df["GC"] >= low_gc) & (df["GC"] < high_gc)]

        for group in list("ABCD"):

            df_subset = df_gc_subset[df_gc_subset["GENOME_TYPE"] == group]
            df_subset = df_subset[df_subset.groupby(f"CONSENSUS_{col}")[f"CONSENSUS_{col}"].transform(len) > 5]

            if len(df_subset) == 0:
                continue

            results_subset = compare_clustering_algorithms_for_gc_and_group(
                env, df_subset
            )

        low_gc = high_gc

def load_gms2_models_from_pickle(pf_mods):
    # type: (str) -> pd.DataFrame
    df = load_obj(pf_mods)  # type: pd.DataFrame

    df["Type"] = "Bacteria"

    df.reset_index(inplace=True)
    df[f"CONSENSUS_RBS_MAT"] = df.apply(lambda r: get_consensus_sequence(r["Mod"].items["RBS_MAT"]), axis=1)

    df["GENOME_TYPE"] = df.apply(lambda r: r["Mod"].items["GENOME_TYPE"], axis=1)
    # df["GC"] = df.apply(lambda r: r["Mod"].items["GC"], axis=1)
    fix_genome_type(df)

    return df


def merge_spacers_by_peak(df):

    shift_to_pos_dist = {0: [df.at[idx, "Mod"].items["RBS_POS_DISTR"] for idx in df.index]}
    # priors on peak positions per shift
    peak_prior_per_shift = dict()  # type: Dict[int, Dict[int, float]]

    position_distributions_by_shift_by_peak = dict()  # type: Dict[int, Dict[int, Dict[int, float]]]

    total_merged = dict()
    for s in sorted(shift_to_pos_dist.keys()):
        position_distributions_by_shift_by_peak[s] = dict()
        peak_prior_per_shift[s] = dict()

        # get list of position distributions for shift
        list_pos_dist = shift_to_pos_dist[s]

        # separate distributions per peak
        peak_to_list_pos_dist = dict()
        for l in list_pos_dist:
            peak = max(l, key=lambda key: l[key])  # get position of peak
            if peak not in peak_to_list_pos_dist:
                peak_to_list_pos_dist[peak] = list()
            peak_to_list_pos_dist[peak].append(l)

        # average positions (per peak)
        values = dict()
        peak_counter = 0
        for peak in peak_to_list_pos_dist.keys():
            values[peak] = dict()
            peak_counter = peak
            peak_prior_per_shift[s][peak_counter] = len(peak_to_list_pos_dist[peak])

            for l in peak_to_list_pos_dist[peak]:
                try:
                    for i in l.keys():
                        if i not in values[peak].keys():
                            values[peak][i] = list()
                        values[peak][i].append(l[i])

                        if i not in total_merged:
                            total_merged[i] = list()
                        total_merged[i].append(l[i])


                except Exception:
                    continue

            for i in values[peak].keys():
                values[peak][i] = np.mean(values[peak][i])

            total = sum(values[peak].values())
            for i in values[peak].keys():
                values[peak][i] /= total

            x = sorted(values[peak].keys())
            y = [values[peak][a] for a in x]

            position_distributions_by_shift_by_peak[s][peak_counter] = {
                a: b for a, b in zip(x, y)
            }

            peak_counter += 1

        # normalize prior for shift
        total = sum(peak_prior_per_shift[s].values())
        if total > 0:
            for pc in peak_prior_per_shift[s].keys():
                peak_prior_per_shift[s][pc] /= float(total)

    for i in total_merged.keys():
        total_merged[i] = np.mean(total_merged[i])

    total = sum(total_merged.values())
    for i in total_merged.keys():
        total_merged[i] /= total

    return position_distributions_by_shift_by_peak[0], peak_prior_per_shift[0], total_merged


def viz_stuff(env, df):
    # type: (Environment, pd.DataFrame) -> None
    low_gc = 0
    for high_gc in range(30, 70, 5):

        df_gc_subset = df[(df["GC"] >= low_gc) & (df["GC"] < high_gc)]

        for group in list("ABCD"):


            for mot in ["AGGAGG", "AAGGAG"]:
                df_subset = df_gc_subset[df_gc_subset["GENOME_TYPE"] == group]
                df_subset = df_subset[df_subset["CONSENSUS_RBS_MAT"] == mot]

                spacer_per_peak, peak_prior, cumul_spacer = merge_spacers_by_peak(df_subset)


                import matplotlib.pyplot as plt

                for j in range(2):

                    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

                    for peak in spacer_per_peak:
                        spacer = spacer_per_peak[peak]
                        spacer = [spacer[i] for i in range(len(spacer))]
                        axes[0].plot(range(len(spacer)), spacer, label=peak)

                    if j == 1:
                        x = sorted(cumul_spacer.keys())
                        y = [cumul_spacer[i] for i in x]
                        axes[0].plot(x, y, label="Merged", linestyle="dashed", color="black")

                    axes[0].set_xlabel("Distance from Gene Start")
                    axes[0].set_ylabel("Frequency")

                    x = sorted(peak_prior)
                    y = [peak_prior[i] for i in x]
                    axes[1].bar(x, y)
                    axes[1].set_xlabel("Peak")
                    axes[1].set_ylabel("Frequency of models with peak")
                    axes[1].set_ylim(0,1)
                    axes[1].set_xlim(0, 15)

                    fig.suptitle(f"{mot}, GC Range [{low_gc}, {high_gc}]")
                    fig.savefig(next_name(env["pd-work"]),  bbox_index="tight")

                    plt.show()

        low_gc = high_gc


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    df = load_gms2_models_from_pickle(args.pf_data)
    # compare_clustering_algorithms(env, df)

    # visualize stuff
    viz_stuff(env, df)





if __name__ == "__main__":
    main(my_env, parsed_args)
