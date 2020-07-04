# Author: Karl Gemayel
# Created: 6/22/20, 3:41 PM
import copy
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
from mg_general.general import get_value, next_name, os_join, fix_names
from mg_io.general import mkdir_p
from mg_models.building import build_mgm_motif_model_for_gc_v2
from mg_models.mgm_motif_model_all_gc import MGMMotifModelAllGC
from mg_models.mgm_motif_model_v2 import MGMMotifModelV2
from mg_models.shelf import read_archaea_bacteria_inputs, get_consensus_sequence, bin_by_gc
from mg_viz.general import FigureOptions

parser = argparse.ArgumentParser("DRIVER DESCRIPTION.")

parser.add_argument('--pf-bac', required=True, help="Collected GMS2 model files for bacteria")
parser.add_argument('--pf-arc', required=True, help="Collected GMS2 model files for archaea")
parser.add_argument('--pf-mgm', required=True, help="Base MGM model file")
parser.add_argument('--pf-output', required=True, help="Output MGM model file")
parser.add_argument('--components', nargs="+", choices=["Start Context", "RBS", "Promoter", "Start Codons"])
parser.add_argument('--genome-type', choices=["Archaea", "Bacteria"])
parser.add_argument('--gc-feature', default="GC")
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


def add_codon_probabilities(df, mgm, codons, gms2_group, **kwargs):
    # type: (pd.DataFrame, MGMModel, List[str], str, Dict[str, Any]) -> None

    genome_type = get_value(kwargs, "genome_type", required=True, choices=["Archaea", "Bacteria"])  # type: str
    plot = get_value(kwargs, "plot", False, valid_type=bool)
    gc_feature = get_value(kwargs, "gc_feature", "GC", valid_type=str)

    genome_tag = genome_type[0]

    gc_step = 1
    gc_min = 30
    gc_max = 71

    if genome_type == "Archaea":
        gc_step = 1
        gc_max = 71

    df = df[df["Type"] == genome_type].copy()

    values_per_codon = dict()
    df.sort_values(gc_feature, inplace=True)
    for c in codons:

        x = df[gc_feature].values
        y = df.apply(lambda r: r["Mod"].items[c], axis=1).astype(float).values
        y = get_loess(x, y)

        values_per_codon[c] = [x, y]



    # if plot:
    # df_tmp = pd.DataFrame(
    #     {
    #         "GC": df[gc_feature].values,
    #         **{c: values_per_codon[c][1] for c in values_per_codon}
    #     }
    # )
    #
    #     df_tmp = pd.melt(df_tmp[["GC"] + codons], ["GC"], var_name="Codon", value_name="Frequency")
    #     fig, ax = plt.subplots(1, 1)
    #     sns.scatterplot(df_tmp, "GC", "Frequency", hue="Codon", ax=ax, show=False,
    #                     sns_kwargs={"alpha": 0.4, "s": 2})
    #
    #     for c in codons:
    #         ax.plot(values_per_codon[c][0], values_per_codon[c][1])
    #
    #     plt.show()

    list_entries = list()

    # get average per GC
    for c in codons:
        current = 0
        x, y = values_per_codon[c]

        x_bin, y_bin = compute_bin_averages(x, y, gc_min, gc_max, gc_step)

        for gc_tag, prob in zip(x_bin, y_bin):
            mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[f"{c}_{gms2_group}"] = prob
            list_entries.append({
                "Codon": c, "GC": gc_tag, "Probability": prob
            })

        # for gc_tag in range(gc_min, gc_max, gc_step):
        #
        #     gc_left = gc_tag if gc_tag != gc_min else 0
        #     gc_right = gc_tag + gc_step if gc_tag != gc_max - gc_step else 100
        #
        #     acc = 0
        #     total = 0
        #
        #     while current < len(x) and gc_left <= x[current] < gc_right:
        #         acc += max(y[current], 0)
        #         total += 1
        #         current += 1
        #
        #     avg = 0 if total == 0 else acc / float(total)
        #     list_entries.append({
        #         "Codon": c, "GC": gc_tag, "Probability": avg
        #     })
        #
        #     if avg == 0:
        #         avg = mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[f"{c}"]
        #
        #     # update MGM
        #     mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[f"{c}_{gms2_group}"] = avg



    if plot:
        df_tmp = pd.DataFrame(list_entries)
        sns.scatterplot(df_tmp, "GC", "Probability", hue="Codon", figure_options=FigureOptions(title=gms2_group))


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

        if total == 0 and len(y_out) == 0:
            continue
        avg = y_out[-1] if total == 0 else acc / float(total)
        x_out.append(x_tag)
        y_out.append(avg)

    return [x_out, y_out]

def add_start_context_probabilities(df, mgm, input_tag, output_tag, **kwargs):
    # type: (pd.DataFrame, MGMModel, str, str, Dict[str, Any]) -> None
    genome_type = get_value(kwargs, "genome_type", required=True, choices=["Archaea", "Bacteria"])
    plot = get_value(kwargs, "plot", False, valid_type=bool)
    gc_feature = get_value(kwargs, "gc_feature", "GC", valid_type=str)


    gc_step = 1
    gc_min = 30
    gc_max = 71

    if genome_type == "Archaea":
        gc_step = 1
        gc_max = 71

    df = df[df["Type"] == genome_type].copy()
    df.sort_values(gc_feature, inplace=True)

    example_sc = df.at[df.index[0], "Mod"].items[input_tag + "_MAT"]  # type: Dict[str, List[float]]
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
        for w in tqdm(words, f"Building {p}", total=len(words)):
            x = [0.0] * len(df.index)
            y = [0.0] * len(df.index)
            for i, idx in enumerate(df.index):
                x[i] = float(df.at[idx, gc_feature])
                y[i] = float(df.at[idx, "Mod"].items[input_tag + "_MAT"][w][p])

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
    # for genome_tag in ["A", "B"]:        # genome_type[0]    FIXME
    genome_tag = genome_type[0]
    for gc_tag in sc_gc.keys():
        mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[output_tag + "_MAT"] = sc_gc[gc_tag]
        mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[f"{output_tag}"] = 1
        mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[f"{output_tag}_ORDER"] = 2
        mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[f"{output_tag}_WIDTH"] = 18
        mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[f"{output_tag}_MARGIN"] = -15

        if gc_step > 1:
            for curr in range(gc_tag, min(gc_max, gc_tag + gc_step)):
                mgm.items_by_species_and_gc[genome_tag][str(curr)].items[output_tag + "_MAT"] = sc_gc[gc_tag]
                mgm.items_by_species_and_gc[genome_tag][str(curr)].items[f"{output_tag}"] = 1
                mgm.items_by_species_and_gc[genome_tag][str(curr)].items[f"{output_tag}_ORDER"] = 2
                mgm.items_by_species_and_gc[genome_tag][str(curr)].items[f"{output_tag}_WIDTH"] = 18
                mgm.items_by_species_and_gc[genome_tag][str(curr)].items[f"{output_tag}_MARGIN"] = -15


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
    gc_feature = get_value(kwargs, "gc_feature", "GC", valid_type=str)


    bin_size = get_value(kwargs, "bin_size", 5, default_if_none=True)

    # get consensus sequences for all motifs
    df[f"CONSENSUS_{name}"] = df.apply(lambda r: get_consensus_sequence(r["Mod"].items[name]), axis=1)

    # bin dataframes by GC
    binned_dfs = bin_by_gc(df, step=bin_size, gc_feature=gc_feature)

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


def add_motif_probabilities(env, df, mgm, input_tag, output_tag, genome_type, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, str, str, str, Dict[str, Any]) -> None
    plot = get_value(kwargs, "plot", False)
    env = env.duplicate({"pd-work": os_join(env["pd-work"], input_tag)})
    mkdir_p(env["pd-work"])
    motif_by_gc = build_mgm_motif_models_for_all_gc(env, df, f"{input_tag}_MAT", plot=plot)


    # width = 6 if tag == "RBS" else 12
    # dur = 14 if tag == "RBS" else 28
    width = df.at[df.index[0], "Mod"].items[f"{input_tag}_WIDTH"]
    dur = df.at[df.index[0], "Mod"].items[f"{input_tag}_MAX_DUR"]

    # tag = "RBS"
    genome_tag = genome_type[0]

    for gc in range(30, 71):
        motif = motif_by_gc.get_model_by_gc(gc)
        best_shift = max(motif._shift_prior.items(), key=operator.itemgetter(1))[0]
        mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}_MAT"] = motif._motif[best_shift]
        mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}_POS_DISTR"] = motif._spacer[best_shift]

        mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}"] = 1
        mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}_ORDER"] = 0
        mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}_WIDTH"] = width
        mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}_MARGIN"] = 0
        mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}_MAX_DUR"] = dur


def build_mgm_models_from_gms2_models(env, df, mgm, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, Dict[str, Any]) -> None

    components = get_value(kwargs, "components", {"Start Codons", "Stop Codons", "Start Context", "RBS", "Promoter"},
                           valid_type=set)
    genome_type = get_value(kwargs, "genome_type", required=True, choices=["Archaea", "Bacteria"])
    plot = get_value(kwargs, "plot", False, valid_type=bool)
    gc_feature = get_value(kwargs, "gc_feature", "GC", valid_type=str)

    # start/stop codons
    if "Start Codons" in components:
        if genome_type == "Archaea":
            output_group = ["A", "D"]
            learn_from = [{"A"}, {"D"}]     # always learn RBS form group A

            for o, l in zip(output_group, learn_from):
                df_curr = df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin(l))]
                add_start_codon_probabilities(df_curr, mgm, genome_type=genome_type, plot=plot, gms2_group=o)
        if genome_type == "Bacteria":
            output_group = ["A", "B", "C", "X"]
            learn_from = [{"A"}, {"B"}, {"C"}, {"A"}]

            for o, l in zip(output_group, learn_from):
                df_curr = df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin(l))]
                add_start_codon_probabilities(df_curr, mgm, genome_type=genome_type, plot=plot, gms2_group=o)
        # add_start_codon_probabilities(df, mgm, genome_type="Archaea", plot=plot)

    # if "Stop Codons" in components:
        # add_stop_codon_probabilities(df, mgm, genome_type=genome_type, plot=plot)
        # add_stop_codon_probabilities(df, mgm, genome_type="Archaea", plot=plot)

    # Start Context
    if "Start Context" in components:
        if genome_type == "Archaea":
            output_group = ["A", "D"]
            learn_from = [{"A"}, {"D"}]     # always learn RBS form group A

            for o, l in zip(output_group, learn_from):
                df_curr = df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin(l))]
                add_start_context_probabilities(df_curr, mgm, "SC_RBS", f"SC_RBS_{o}", genome_type=genome_type, plot=plot)
        else:
            output_group = ["A", "B", "C", "X"]
            learn_from = [{"A"}, {"B"}, {"C"}, {"A"}]

            for o, l in zip(output_group, learn_from):
                df_curr = df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin(l))]
                add_start_context_probabilities(df_curr, mgm, "SC_RBS", f"SC_RBS_{o}", genome_type=genome_type, plot=plot)


        # promoter
        if genome_type == "Archaea":
            output_group = ["D"]
            learn_from = [{"D"}]     # always learn RBS form group A

            for o, l in zip(output_group, learn_from):
                df_curr = df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin(l))]
                add_start_context_probabilities(df_curr, mgm, "SC_PROMOTER", f"SC_PROMOTER_{o}", genome_type=genome_type, plot=plot)
        else:
            output_group = ["C"]
            learn_from = [{"C"}]

            for o, l in zip(output_group, learn_from):
                df_curr = df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin(l))]
                add_start_context_probabilities(df_curr, mgm, "SC_PROMOTER", f"SC_PROMOTER_{o}", genome_type=genome_type, plot=plot)



    # Motifs
    if "RBS" in components:
        if genome_type == "Archaea":
            output_group = ["A", "D"]
            learn_from = [{"A"}, {"D"}]     # always learn RBS form group A

            df_type = df[df["Type"] == genome_type]
            for o, l in zip(output_group, learn_from):
                add_motif_probabilities(
                    env,
                    df_type[(df_type["GENOME_TYPE"].isin(l))],
                    mgm,
                    "RBS", f"RBS_{o}", genome_type, plot=plot
                )
        else:
            output_group = ["A", "B", "C", "X"]
            learn_from = [{"A"}, {"B"}, {"C"}, {"A"}]

            df_type = df[df["Type"] == genome_type]
            for o, l in zip(output_group, learn_from):
                add_motif_probabilities(
                    env,
                    df_type[(df_type["GENOME_TYPE"].isin(l))],
                    mgm,
                    "RBS", f"RBS_{o}", genome_type, plot=plot
                )

    if "Promoter" in components:
        if genome_type == "Archaea":
            output_group = ["D"]
            learn_from = [{"D"}]  # always learn Promoter form group D

            df_type = df[df["Type"] == genome_type]
            for o, l in zip(output_group, learn_from):
                add_motif_probabilities(
                    env,
                    df_type[(df_type["GENOME_TYPE"].isin(l))],
                    mgm,
                    "PROMOTER", f"PROMOTER_{o}", genome_type, plot=plot
                )
        else:
            output_group = ["C"]
            learn_from = [{"C"}]  # always learn Promoter form group C

            df_type = df[df["Type"] == genome_type]
            for o, l in zip(output_group, learn_from):
                add_motif_probabilities(
                    env,
                    df_type[(df_type["GENOME_TYPE"].isin(l))],
                    mgm,
                    "PROMOTER", f"PROMOTER_{o}", genome_type, plot=plot
                )


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    df = read_archaea_bacteria_inputs(args.pf_arc, args.pf_bac)
    logger.debug(f"Removing genetic code 4: {(df['Genetic Code'] == 4).sum()}")
    df = df[df["Genetic Code"] != 4]
    df = df.convert_dtypes().copy()
    # df = df.sample(100).copy()

    # sns.scatterplot(df, "GC", "Average Gene GC")
    sns.scatterplot(df[df["Type"] == "Archaea"], "GC", "Average Gene GC")
    sns.scatterplot(df[df["Type"] == "Archaea"], "GC", "Median Gene GC")

    df_opg = df.groupby("Genome", as_index=False).first()
    for name, df_gcfid in df_opg.groupby("GENOME_TYPE", as_index=False):

        sns.distplot(df_gcfid, "GC", figure_options=FigureOptions(title=name))

    gc_feature = args.gc_feature

    pd_work = os_join(env["pd-work"], "_".join(args.components).lower())
    mkdir_p(pd_work)
    env = env.duplicate({'pd-work': pd_work})

    mgm = MGMModel.init_from_file(args.pf_mgm)
    for genome_type in {"Archaea", "Bacteria"}:
        build_mgm_models_from_gms2_models(env, df, mgm, components=args.components, genome_type=genome_type,
                                      plot=args.plot, gc_feature=gc_feature)
    mgm.to_file(args.pf_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
