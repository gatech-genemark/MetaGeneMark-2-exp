# Author: Karl Gemayel
# Created: 6/22/20, 3:41 PM

import logging
import argparse
import math
import operator

import numpy as np
import pandas as pd
from typing import *
from tqdm import tqdm
from statsmodels import api as sm
import matplotlib.pyplot as plt

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
import mg_viz.sns as sns
from mg_container.mgm_model import MGMModel, MGMModelGC
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import get_value, next_name, os_join
from mg_io.general import mkdir_p
from mg_models.building import build_mgm_motif_model_for_gc_v2
from mg_models.mgm_motif_model_all_gc import MGMMotifModelAllGC
from mg_models.mgm_motif_model_v2 import MGMMotifModelV2
from mg_models.shelf import read_archaea_bacteria_inputs, get_consensus_sequence, bin_by_gc

parser = argparse.ArgumentParser("DRIVER DESCRIPTION.")

parser.add_argument('--pf-bac', required=True, help="Collected GMS2 model files for bacteria")
parser.add_argument('--pf-arc', required=True, help="Collected GMS2 model files for archaea")
parser.add_argument('--pf-mgm', required=True, help="Base MGM model file")
parser.add_argument('--pf-output', required=True, help="Output MGM model file")
parser.add_argument('--components', nargs="+", choices=["Start Context", "RBS", "Promoter", "Start Codons"])
parser.add_argument('--genome-type', choices=["Archaea", "Bacteria"])
parser.add_argument('--plot', default=False, action="store_true")
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


def get_loess(local_x, local_y):
    loess = sm.nonparametric.lowess(local_y, local_x)
    return loess[:, 1]


def add_codon_probabilities(df, mgm, codons, **kwargs):
    # type: (pd.DataFrame, MGMModel, List[str], Dict[str, Any]) -> None

    genome_type = get_value(kwargs, "genome_type", required=True, choices=["Archaea", "Bacteria"])
    plot = get_value(kwargs, "plot", False, valid_type=bool)

    gc_step = 1
    gc_min = 30
    gc_max = 71

    if genome_type == "Archaea":
        gc_step = 5
        gc_max = 66

    df = df[df["Type"] == genome_type].copy()

    values_per_codon = dict()
    df.sort_values("GC", inplace=True)
    for c in codons:

        x = df["GC"].values
        y = df.apply(lambda r: r["Mod"].items[c], axis=1).astype(float).values
        y = get_loess(x, y)

        values_per_codon[c] = [x, y]

    df_tmp = pd.DataFrame(
        {
            "GC": df["GC"].values,
            **{c: values_per_codon[c][1] for c in values_per_codon}
        }
    )
    df_tmp = pd.melt(df_tmp[["GC"] + codons], ["GC"], var_name="Codon", value_name="Frequency")

    if plot:
        fig, ax = plt.subplots(1, 1)
        sns.scatterplot(df_tmp, "GC", "Frequency", hue="Codon", ax=ax, show=False,
                        sns_kwargs={"alpha": 0.4, "s": 2})

        for c in codons:
            ax.plot(values_per_codon[c][0], values_per_codon[c][1])

        plt.show()

    list_entries = list()

    # get average per GC
    for c in codons:
        current = 0
        x, y = values_per_codon[c]
        for gc_tag in range(gc_min, gc_max, gc_step):

            gc_left = gc_tag if gc_tag != gc_min else 0
            gc_right = gc_tag + gc_step if gc_tag != gc_max - gc_step else 100

            acc = 0
            total = 0

            while current < len(x) and gc_left <= x[current] < gc_right:
                acc += max(y[current], 0)
                total += 1
                current += 1

            avg = 0 if total == 0 else acc / float(total)
            list_entries.append({
                "Codon": c, "GC": gc_tag, "Probability": avg
            })

            # update MGM
            genome_tag = "B"        # FIXME
            mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[c] = avg

    df_tmp = pd.DataFrame(list_entries)

    if plot:
        sns.scatterplot(df_tmp, "GC", "Probability", hue="Codon")


def add_start_codon_probabilities(df, mgm, **kwargs):
    # type: (pd.DataFrame, MGMModel, Dict[str, Any]) -> None
    add_codon_probabilities(df, mgm, ["ATG", "GTG", "TTG"], **kwargs)


def add_stop_codon_probabilities(df, mgm, **kwargs):
    # type: (pd.DataFrame, MGMModel, Dict[str, Any]) -> None
    add_codon_probabilities(df, mgm, ["TAA", "TAG", "TGA"], **kwargs)


def compute_bin_averages(x, y, x_min, x_max, x_step):
    # type: (List[float], List[float], float, float, float) -> [List[float], List[float]]

    x_out = list()
    y_out = list()

    current = 0
    for x_tag in np.arange(x_min, x_max, x_step):


        gc_left = x_tag if x_tag != x_min else 0
        gc_right = x_tag + x_step if x_tag != x_max - x_step else 100

        acc = 0
        total = 0

        while current < len(x) and gc_left <= x[current] < gc_right:
            acc += max(y[current], 0)
            total += 1
            current += 1

        avg = y_out[-1] if total == 0 else acc / float(total)
        x_out.append(x_tag)
        y_out.append(avg)

    return [x_out, y_out]

def add_start_context_probabilities(df, mgm, tag, **kwargs):
    # type: (pd.DataFrame, MGMModel, str, Dict[str, Any]) -> None
    genome_type = get_value(kwargs, "genome_type", required=True, choices=["Archaea", "Bacteria"])
    plot = get_value(kwargs, "plot", False, valid_type=bool)

    gc_step = 1
    gc_min = 30
    gc_max = 71

    if genome_type == "Archaea":
        gc_step = 1
        gc_max = 71

    df = df[df["Type"] == genome_type].copy()
    df.sort_values("GC", inplace=True)

    example_sc = df.at[df.index[0], "Mod"].items[tag + "_MAT"]  # type: Dict[str, List[float]]
    words = sorted(set(example_sc.keys()))
    num_positions = len(next(iter(example_sc.values())))

    # create empty models for each gc bin
    sc_gc = dict()
    for gc_tag in np.arange(gc_min, gc_max, gc_step):
        sc_gc[gc_tag] = {
            w: [0] * num_positions for w in words
        }

    list_entries = dict()

    # get all words appearing in start contexts and all positions

    for p in range(num_positions):
        # if p > 1:
        #     break
        for w in tqdm(words, f"Building {p}", total=len(words)):
            x = [0.0] * len(df.index)
            y = [0.0] * len(df.index)
            for i, idx in enumerate(df.index):
                x[i] = float(df.at[idx, "GC"])
                y[i] = float(df.at[idx, "Mod"].items[tag + "_MAT"][w][p])

            list_entries["GC"] = x
            list_entries[f"{w}{p}"] = y

    df_tmp = pd.DataFrame(list_entries)
    for p in range(num_positions):
        num_words = len(words)
        num_rows = int(math.sqrt(num_words))
        num_cols = math.ceil(num_words / float(num_rows))

        if plot:
            fig, axes = plt.subplots(num_rows, num_cols, sharex="all", sharey="all", figsize=(20,20))

        for i, w in tqdm(enumerate(words), f"Plotting {p}", total=len(words)):

            if plot:
                ax = axes.ravel()[i]

            # sns.scatterplot(
            #     df_tmp, "GC", f"{w}{p}", sns_kwargs={"alpha": 0.3, "s": 2},
            #     ax=ax, show=False
            # )

            x = df_tmp["GC"].values
            y = df_tmp[f"{w}{p}"].values

            if plot:
                ax.scatter(x, y, alpha=0.3, s=2 if genome_type == "Bacteria" else 4)

            y = get_loess(x, y)

            x_bin, y_bin = compute_bin_averages(x, y, gc_min, gc_max, gc_step)

            # fill gc models
            for gc_tag, prob in zip(x_bin, y_bin):
                sc_gc[gc_tag][w][p] = prob

            if plot:
                ax.plot(x_bin, y_bin, "r")
                ax.set_title(w)
                ax.set_ylim(0, 1)

        if plot:
            fig.suptitle(f"Position {p}")
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            plt.savefig(next_name(my_env["pd-work"], ext="png"))
            plt.show()

    # add gc models to mgm
    genome_tag = "B"        # genome_type[0]    FIXME
    for gc_tag in sc_gc.keys():
        mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[tag + "_MAT"] = sc_gc[gc_tag]
        mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[f"{tag}"] = 1
        mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[f"{tag}_ORDER"] = 2
        mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[f"{tag}_WIDTH"] = 18
        mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[f"{tag}_MARGIN"] = -15

        if gc_step > 1:
            for curr in range(gc_tag, min(gc_max, gc_tag + gc_step)):
                mgm.items_by_species_and_gc[genome_tag][str(curr)].items[tag + "_MAT"] = sc_gc[gc_tag]
                mgm.items_by_species_and_gc[genome_tag][str(curr)].items[f"{tag}"] = 1
                mgm.items_by_species_and_gc[genome_tag][str(curr)].items[f"{tag}_ORDER"] = 2
                mgm.items_by_species_and_gc[genome_tag][str(curr)].items[f"{tag}_WIDTH"] = 18
                mgm.items_by_species_and_gc[genome_tag][str(curr)].items[f"{tag}_MARGIN"] = -15


    #
    #         values_per_word = dict()
    #         for w in words:
    #             x = [] * len(df.index)
    #         y = [] * len(df.index)
    #         for i, idx in enumerate(df.index):
    #             x[i] = df.at[idx, "GC"]
    #         y[i] = df.at[idx, tag][w]
    #
    #         values_per_codon = dict()
    #         for c in codons:
    #             df[c] = df[c].astype(float)
    #         x = df["GC"].values
    #         y = df[c].values
    #         y = get_loess(x, y)
    #
    #         values_per_codon[c] = [x, y]
    #
    #         df_tmp = pd.melt(df[["GC"] + codons], ["GC"], var_name="Codon", value_name="Frequency")
    #
    #         fig, ax = plt.subplots(1, 1)
    #         sns.scatterplot(df_tmp, "GC", "Frequency", hue="Codon", ax=ax, show=False,
    #                         sns_kwargs={"alpha": 0.4, "s": 2})
    #
    #         for c in codons:
    #             ax.plot(values_per_codon[c][0], values_per_codon[c][1])
    #
    #         plt.show()
    #
    #         list_entries = list()
    #
    #         # get average per GC
    #         for c in codons:
    #             current = 0
    #         x, y = values_per_codon[c]
    #         for gc_tag in range(gc_min, gc_max, gc_step):
    #             if
    #         gc_tag == gc_max - gc_step:
    #         print('hi')
    #
    #         gc_left = gc_tag if gc_tag != gc_min else 0
    #         gc_right = gc_tag + gc_step if gc_tag != gc_max - gc_step else 100
    #
    #         acc = 0
    #         total = 0
    #
    #         while current < len(x) and gc_left <= x[current] < gc_right:
    #             acc += max(y[current], 0)
    #             total += 1
    #             current += 1
    #
    #         avg = 0 if total == 0 else acc / float(total)
    #         list_entries.append({
    #             "Codon": c, "GC": gc_tag, "Probability": avg
    #         })
    #
    #         # update MGM
    #         mgm.items_by_species_and_gc[genome_type[0]][str(gc_tag)].items[c] = avg
    #
    # df_tmp = pd.DataFrame(list_entries)
    #
    # sns.scatterplot(df_tmp, "GC", "Probability", hue="Codon")


def build_mgm_motif_models_for_all_gc(env, df, name, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> MGMMotifModelAllGC
    has_name = df.apply(lambda r: name in r["Mod"].items, axis=1)
    df = df[has_name].copy()  # we only need non-NA

    bin_size = get_value(kwargs, "bin_size", 5, default_if_none=True)

    # get consensus sequences for all motifs
    df[f"CONSENSUS_{name}"] = df.apply(lambda r: get_consensus_sequence(r["Mod"].items[name]), axis=1)

    # bin dataframes by GC
    binned_dfs = bin_by_gc(df, step=bin_size)

    # for each binned dataframe, build specific model
    list_mgm_models = list()  # type: List[List[float, float, MGMMotifModelV2]]
    for info in binned_dfs:
        lower, upper, df_gc = info

        # if int(lower) != 40:
        #     continue

        mgm_mm = None
        if len(df_gc) > 1:
            mgm_mm = build_mgm_motif_model_for_gc_v2(env, df_gc, name, title=f"[{lower},{upper}]", **kwargs)

        if lower == 30 and upper == 35 and mgm_mm is None:
            print('hi')

        if mgm_mm is None:
            # use previous model
            if len(list_mgm_models) > 0:
                prev = list_mgm_models[-1][2]
                list_mgm_models.append([lower, upper, prev])
        else:
            list_mgm_models.append([lower, upper, mgm_mm])

    return MGMMotifModelAllGC(list_mgm_models)


def add_motif_probabilities(env, df, mgm, tag, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, str, Dict[str, Any]) -> None
    plot = get_value(kwargs, "plot", False)
    env = env.duplicate({"pd-work": os_join(env["pd-work"], tag)})
    mkdir_p(env["pd-work"])
    motif_by_gc = build_mgm_motif_models_for_all_gc(env, df, f"{tag}_MAT", plot=plot)


    # width = 6 if tag == "RBS" else 12
    # dur = 14 if tag == "RBS" else 28
    width = df.at[df.index[0], "Mod"].items[f"{tag}_WIDTH"]
    dur = df.at[df.index[0], "Mod"].items[f"{tag}_MAX_DUR"]

    # tag = "RBS"     # FIXME
    for gc in range(30, 71):
        motif = motif_by_gc.get_model_by_gc(gc)
        best_shift = max(motif._shift_prior.items(), key=operator.itemgetter(1))[0]
        mgm.items_by_species_and_gc["B"][str(gc)].items[f"{tag}_MAT"] = motif._motif[best_shift]
        mgm.items_by_species_and_gc["B"][str(gc)].items[f"{tag}_POS_DISTR"] = motif._spacer[best_shift]

        mgm.items_by_species_and_gc["B"][str(gc)].items[f"{tag}"] = 1
        mgm.items_by_species_and_gc["B"][str(gc)].items[f"{tag}_ORDER"] = 0
        mgm.items_by_species_and_gc["B"][str(gc)].items[f"{tag}_WIDTH"] = width
        mgm.items_by_species_and_gc["B"][str(gc)].items[f"{tag}_MARGIN"] = 0
        mgm.items_by_species_and_gc["B"][str(gc)].items[f"{tag}_MAX_DUR"] = dur


def build_mgm_models_from_gms2_models(env, df, mgm, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, Dict[str, Any]) -> None

    components = get_value(kwargs, "components", {"Start Codons", "Stop Codons", "Start Context", "RBS", "Promoter"},
                           valid_type=set)
    genome_type = get_value(kwargs, "genome_type", required=True, choices=["Archaea", "Bacteria"])
    plot = get_value(kwargs, "plot", False, valid_type=bool)

    # start/stop codons
    if "Start Codons" in components:
        add_start_codon_probabilities(df, mgm, genome_type=genome_type, plot=plot)
        # add_start_codon_probabilities(df, mgm, genome_type="Archaea", plot=plot)

    if "Stop Codons" in components:
        add_stop_codon_probabilities(df, mgm, genome_type=genome_type, plot=plot)
        # add_stop_codon_probabilities(df, mgm, genome_type="Archaea", plot=plot)

    # Start Context
    if "Start Context" in components:
        add_start_context_probabilities(df, mgm, "SC_RBS", genome_type=genome_type, plot=plot)
        add_start_context_probabilities(df, mgm, "SC_PROMOTER", genome_type=genome_type, plot=plot)
        # add_start_context_probabilities(df, mgm, "SC_RBS", genome_type="Archaea")

    # Motifs
    if "RBS" in components:
        if genome_type == "Archaea":
            add_motif_probabilities(
                env,
                df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin({"A", "D"}))],
                mgm,
                "RBS", plot=plot
            )
        else:
            add_motif_probabilities(
                env,
                df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin({"A", "C"}))],
                mgm,
                "RBS", plot=plot
            )

    if "Promoter" in components:
        if genome_type == "Archaea":
            add_motif_probabilities(
                env,
                df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin({"D"}))],
                mgm,
                "PROMOTER", plot=plot
            )
        else:
            add_motif_probabilities(
                env,
                df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin({"C"}))],
                mgm,
                "PROMOTER", plot=plot
            )

def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    df = read_archaea_bacteria_inputs(args.pf_arc, args.pf_bac)
    df = df.convert_dtypes().copy()
    # df = df.sample(100).copy()

    mgm = MGMModel.init_from_file(args.pf_mgm)
    build_mgm_models_from_gms2_models(env, df, mgm, components=args.components, genome_type=args.genome_type,
                                      plot=args.plot)
    mgm.to_file(args.pf_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
