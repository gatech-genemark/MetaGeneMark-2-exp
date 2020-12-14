# Author: Karl Gemayel
# Created: 6/22/20, 3:41 PM
import logging
import argparse
import math
import operator

import numpy as np
import pandas as pd
from typing import *

import seaborn
import yaml
from tqdm import tqdm
from statsmodels import api as sm
import matplotlib.pyplot as plt

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.mgm_model import MGMModel
from mg_general import Environment, add_env_args_to_parser
from mg_general.general import get_value, next_name, os_join, except_if_not_valid
from mg_io.general import mkdir_p
from mg_models.building import build_mgm_motif_model_for_gc_v2
from mg_models.mgm_motif_model_all_gc import MGMMotifModelAllGC
from mg_models.mgm_motif_model_v2 import MGMMotifModelV2
from mg_models.shelf import read_archaea_bacteria_inputs, get_consensus_sequence, bin_by_gc
from mg_options.learn_from import LearnFromOptions
from mg_viz.general import FigureOptions, save_figure

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #


parser = argparse.ArgumentParser("DRIVER DESCRIPTION.")

parser.add_argument('--pf-bac', required=True, help="Collected GMS2 model files for bacteria")
parser.add_argument('--pf-arc', required=True, help="Collected GMS2 model files for archaea")
parser.add_argument('--pf-mgm', required=True, help="Base MGM model file to update")
parser.add_argument('--pf-output', required=True, help="Output MGM model file")
parser.add_argument('--components', nargs="+",
                    choices=["Start Context", "RBS", "Promoter", "Start Codons", "Stop Codons"])
parser.add_argument('--genome-type', choices=["Archaea", "Bacteria"], default=None,
                    help="Set if only want to build for single set. Leave empty for both")

parser.add_argument('--cluster-by', default="msa", choices=["msa", "heuristic"], type=str.lower)

parser.add_argument('--gc-feature', default="GC")
parser.add_argument('--plot', default=False, action="store_true")
parser.add_argument('--pf-learn-from-options')

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


# class LearnFrom:
#     """Maps which group of genomes to learn from, for each feature
#     E.g. It can specify that Start codons for archaea group A should be learned only from
#     group A, or from group A and D.
#
#     For now, allowed groups are hard coded: Bacteria (A, B, C, X) and Archaea (A, D)
#     """
#
#     def __init__(self, info=None):
#         # type: (Dict[str, Dict[str, Dict[str, Set[str]]]]) -> None
#         self._values = LearnFrom._default_values()
#
#         if info is not None:
#             # update values based on input
#             for component, c_vals in self._values.items():
#                 if component in info:
#                     for gtype, g_vals in c_vals.items():
#                         if gtype in info[component]:
#                             for group, gr_vals in g_vals.items():
#                                 if group in info[component][gtype]:
#                                     self._values[component][gtype][group] = info[component][gtype][group]
#
#     def __getitem__(self, item):
#         # type: (str) -> Dict[str, Dict[str, Set[str]]]
#         return self._values[item]
#
#     @staticmethod
#     def _default_values():
#         # type: () -> Dict[str, Dict[str, Dict[str, Set[str]]]]
#         return {
#             component: {
#                 "Archaea": {"A": {"A"}, "D": {"D"}},
#                 "Bacteria": {"A": {"A"}, "B": {"B"}, "C": {"C"}, "X": {"X"}}
#             } for component in {"RBS", "PROMOTER", "Start Context", "Start Codons", "Stop Codons"}
#         }
#
#     @classmethod
#     def init_from_file(cls, pf_config):
#         # type: (str) -> LearnFrom
#         try:
#             f = open(pf_config, "r")
#             return LearnFrom(yaml.load(f, Loader=yaml.FullLoader))
#         except IOError:
#             logger.warning(f"Configuration File Not Found: {pf_config}. "
#                            f"Using defaults.")
#             return LearnFrom()


def get_loess(local_x, local_y):
    loess = sm.nonparametric.lowess(local_y, local_x)
    return loess[:, 1]


def visualize_start_codons(env, viz_collector):
    # type: (Environment, Dict[str, Dict[str, Dict[str, Any]]]) -> None

    list_entries = list()

    for genome_type in viz_collector:
        for group in viz_collector[genome_type]:
            if group == "X":
                continue
            for codon in viz_collector[genome_type][group]:
                vals = viz_collector[genome_type][group][codon]
                x = vals["x"]
                y = vals["y"]
                y_fit = vals["y_fit"]

                for i in range(len(x)):
                    list_entries.append({
                        "Genome Type": genome_type,
                        "Group": group if genome_type == "Bacteria" else f"A*,D*",
                        "Codon": codon,
                        "x": x[i],
                        "y": y[i],
                        "y_fit": y_fit[i]
                    })
            if genome_type == "Archaea":
                break

    df = pd.DataFrame(list_entries)
    g = seaborn.FacetGrid(df, col="Codon", hue="Group")
    g.map(plt.scatter, "x", "y", alpha=.3, s=2)
    g.map(plt.plot, "x", "y_fit")
    g.set_xlabels("GC")
    g.set_ylabels("Probability")
    g.add_legend()
    g.set_titles(col_template='{col_name}')
    leg = g._legend
    for lh in leg.legendHandles:
        lh.set_alpha(1)
        lh.set_sizes([14] * 3)
    #
    # g.fig.subplots_adjust(top=.8)
    # g.fig.suptitle(genome_type)


    g.fig.savefig(next_name(env["pd-work"]))
    plt.close()

    # plt.show()


def visualize_stop_codons(env, viz_collector):
    visualize_start_codons(env, viz_collector)


def add_codon_probabilities(env, df, mgm, codons, gms2_group, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, List[str], str, Dict[str, Any]) -> None

    genome_type = get_value(kwargs, "genome_type", required=True, choices=["Archaea", "Bacteria"])  # type: str
    plot = get_value(kwargs, "plot", False, valid_type=bool)
    gc_feature = get_value(kwargs, "gc_feature", "GC", valid_type=str)
    viz_collector = get_value(kwargs, "viz_collector", None)

    genome_tag = genome_type[0]

    gc_step = 1
    gc_min = 30
    gc_max = 71

    if genome_type == "Archaea":
        gc_step = 1
        gc_max = 71

    df = df[df["Type"] == genome_type].copy()

    list_entries = list()

    fig, ax = plt.subplots()
    # values_per_codon = dict()
    df.sort_values(gc_feature, inplace=True)
    for c in codons:

        x = df[gc_feature].values
        y_original = df.apply(lambda r: r["Mod"].items[c], axis=1).astype(float).values
        y = get_loess(x, y_original)

        # values_per_codon[c] = [x, y]

        # x, y = values_per_codon[c]

        x_bin, y_bin = compute_bin_averages(x, y, gc_min, gc_max, gc_step)

        for gc_tag, prob in zip(x_bin, y_bin):
            mgm.items_by_species_and_gc[genome_tag][str(gc_tag)].items[f"{c}_{gms2_group}"] = prob
            list_entries.append({
                "Codon": c, "GC": gc_tag, "Probability": prob
            })

        ax.scatter(x, y_original, alpha=0.3, s=2 if genome_type == "Bacteria" else 4, label=c)
        ax.plot(x, y)

        if viz_collector is not None:
            viz_collector[c] = dict()
            viz_collector[c]["x"] = x
            viz_collector[c]["y"] = y_original
            viz_collector[c]["y_fit"] = y

    plt.title(f"{genome_type}: {gms2_group}")
    leg = ax.legend()
    for lh in leg.legendHandles:
        lh.set_alpha(1)
        lh.set_sizes([4] * 3)
    plt.xlabel("GC")
    plt.ylabel("Probability")
    plt.xlim([20, 80])
    plt.ylim([0, 1])
    plt.savefig(next_name(env["pd-work"]))
    plt.close()
    # plt.show()

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

    # if plot:
    #     df_tmp = pd.DataFrame(list_entries)
    #     sns.scatterplot(df_tmp, "GC", "Probability", hue="Codon", figure_options=FigureOptions(title=gms2_group))


def add_start_codon_probabilities(env, df, mgm, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, Dict[str, Any]) -> None
    add_codon_probabilities(env, df, mgm, ["ATG", "GTG", "TTG"], **kwargs)


def add_stop_codon_probabilities(env, df, mgm, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, Dict[str, Any]) -> None
    add_codon_probabilities(env, df, mgm, ["TAA", "TAG", "TGA"], **kwargs)


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


def add_start_context_probabilities(env, df, mgm, input_tag, output_tag, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, str, str, Dict[str, Any]) -> None
    genome_type = get_value(kwargs, "genome_type", required=True, choices=["Archaea", "Bacteria"])
    plot = get_value(kwargs, "plot", False, valid_type=bool)
    gc_feature = get_value(kwargs, "gc_feature", "GC", valid_type=str)
    pd_figures = get_value(kwargs, "pd_figures", env["pd-work"])
    pd_figures = os_join(pd_figures, input_tag)
    mkdir_p(pd_figures)

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
        for w in words:
            x = [0.0] * len(df.index)
            y = [0.0] * len(df.index)
            for i, idx in enumerate(df.index):
                x[i] = float(df.at[idx, gc_feature])
                y[i] = float(df.at[idx, "Mod"].items[input_tag + "_MAT"][w][p])

            list_entries["GC"] = x
            list_entries[f"{w}{p}"] = y

    df_tmp = pd.DataFrame(list_entries)
    for p in tqdm(range(num_positions), f"Building {input_tag}", total=num_positions):
        num_words = len(words)
        num_rows = int(math.sqrt(num_words))
        num_cols = math.ceil(num_words / float(num_rows))

        if plot:
            fig, axes = plt.subplots(num_rows, num_cols, sharex="all", sharey="all", figsize=(20, 20))

        for i, w in enumerate(words):

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
            plt.savefig(next_name(pd_figures, ext="png"))
            plt.close()
            # plt.show()

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
        #
        # if int(lower) != 45:
        #     continue

        mgm_mm = None
        if len(df_gc) > 1:
            mgm_mm = build_mgm_motif_model_for_gc_v2(env, df_gc, name, title=f"[{lower},{upper}]", **kwargs)

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
    pd_figures = get_value(kwargs, "pd_figures", env["pd-work"])
    pd_figures = os_join(pd_figures, input_tag)
    cluster_by = get_value(kwargs, "cluster_by", "msa")
    # env = env.duplicate({"pd-work": os_join(env["pd-work"], input_tag)})
    mkdir_p(pd_figures)
    motif_by_gc = build_mgm_motif_models_for_all_gc(env, df, f"{input_tag}_MAT", plot=plot, pd_figures=pd_figures,
                                                    cluster_by=cluster_by)

    # width = 6 if tag == "RBS" else 12
    # dur = 14 if tag == "RBS" else 28
    width = df.at[df.index[0], "Mod"].items[f"{input_tag}_WIDTH"]
    dur = df.at[df.index[0], "Mod"].items[f"{input_tag}_MAX_DUR"]

    # tag = "RBS"
    genome_tag = genome_type[0]

    for gc in range(30, 71):

        motif = motif_by_gc.get_model_by_gc(gc)

        if True or "RBS" in output_tag:
            # create a label for each shift
            for shift, prob in motif._shift_prior.items():
                prob /= 100.0
                output_tag_ws = f"{output_tag}_{int(shift)}"
                try:
                    mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag_ws}_MAT"] = motif._motif[shift]
                    mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag_ws}_POS_DISTR"] = \
                        motif._spacer[
                            shift]
                except KeyError:
                    pass

                mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag_ws}"] = 1
                mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag_ws}_ORDER"] = 0
                mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag_ws}_WIDTH"] = width
                mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag_ws}_MARGIN"] = 0
                mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag_ws}_MAX_DUR"] = dur
                mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag_ws}_SHIFT"] = prob
        else:
            # promoter aren't shifted (for now)
            best_shift = max(motif._shift_prior.items(), key=operator.itemgetter(1))[0]
            mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}_MAT"] = motif._motif[best_shift]
            mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}_POS_DISTR"] = motif._spacer[
                best_shift]

            mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}"] = 1
            mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}_ORDER"] = 0
            mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}_WIDTH"] = width
            mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}_MARGIN"] = 0
            mgm.items_by_species_and_gc[genome_tag][str(gc)].items[f"{output_tag}_MAX_DUR"] = dur


def _build_start_or_stop_codons(env, df, mgm, genome_type, codons, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, str, List[str], Dict[str, Any]) -> None
    plot = get_value(kwargs, "plot", False, valid_type=bool)
    learn_from = get_value(kwargs, "learn_from", default_value_callable=LearnFromOptions)  # type: LearnFromOptions

    learn_from_component = learn_from["Start Codons" if "ATG" in codons else "Stop Codons"]  # get for component

    viz_collector = dict()
    if genome_type == "Archaea":

        viz_collector[genome_type] = dict()
        for o, l in learn_from_component[genome_type].items():
            viz_collector[genome_type][o] = dict()
            df_curr = df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin(l))]
            add_codon_probabilities(env, df_curr, mgm, codons, genome_type=genome_type, plot=plot, gms2_group=o,
                                    viz_collector=viz_collector[genome_type][o])

    if genome_type == "Bacteria":

        viz_collector[genome_type] = dict()

        for o, l in learn_from_component[genome_type].items():
            viz_collector[genome_type][o] = dict()
            df_curr = df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin(l))]
            add_codon_probabilities(env, df_curr, mgm, codons, genome_type=genome_type, plot=plot,
                                    gms2_group=o,
                                    viz_collector=viz_collector[genome_type][o])

    if plot:
        pd_figures = get_value(kwargs, "pd_figures", env["pd-work"])
        visualize_start_codons(env.duplicate({"pd-work": pd_figures}), viz_collector)


def _build_start_context(env, df, mgm, genome_type, tag, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, str, str, Dict[str, Any]) -> None

    learn_from = get_value(kwargs, "learn_from", default_value_callable=LearnFromOptions)  # type: LearnFromOptions
    except_if_not_valid(tag, {"SC_RBS", "SC_PROMOTER"})

    learn_from_component = learn_from["Start Context"]  # get for component

    if genome_type == "Archaea":
        for o, l in learn_from_component[genome_type].items():
            if "PROMOTER" in tag and o != "D":
                continue  # promoters are only in group D
            df_curr = df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin(l))]

            # NOTE: SC_PROMOTER is intentionally learned from SC_RBS. This is not a bug
            # GMS2 has equal values for SC_RBS and SC_PROMOTER. Training from SC_RBS allows us
            # to learn from group A genomes as well (if needed).
            add_start_context_probabilities(env, df_curr, mgm, "SC_RBS", f"{tag}_{o}", genome_type=genome_type,
                                            **kwargs)
    else:
        # Bacteria
        for o, l in learn_from_component[genome_type].items():
            if "PROMOTER" in tag and o != "C":
                continue  # promoters are only in group C
            df_curr = df[(df["Type"] == genome_type) & (df["GENOME_TYPE"].isin(l))]

            # NOTE: SC_PROMOTER is intentionally learned from SC_RBS. This is not a bug
            # GMS2 has equal values for SC_RBS and SC_PROMOTER. Training from SC_RBS therefore allows us
            # to learn from group A genomes as well.
            add_start_context_probabilities(env, df_curr, mgm, "SC_RBS", f"{tag}_{o}", genome_type=genome_type,
                                            **kwargs)


def _build_motifs(env, df, mgm, genome_type, tag, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, str, str, Dict[str, any]) -> None
    learn_from = get_value(kwargs, "learn_from", default_value_callable=LearnFromOptions)  # type: LearnFromOptions
    except_if_not_valid(tag, {"RBS", "PROMOTER"})

    learn_from_component = learn_from[tag]  # get for component

    if genome_type == "Archaea":

        df_type = df[df["Type"] == genome_type]
        for o, l in learn_from_component[genome_type].items():
            if "PROMOTER" in tag and o != "D":
                continue  # promoters are only in group D

            add_motif_probabilities(
                env,
                df_type[(df_type["GENOME_TYPE"].isin(l))],
                mgm,
                f"{tag}", f"{tag}_{o}", genome_type, **kwargs
            )
    else:

        df_type = df[df["Type"] == genome_type]
        for o, l in learn_from_component[genome_type].items():
            if "PROMOTER" in tag and o != "C":
                continue  # promoters are only in group C

            add_motif_probabilities(
                env,
                df_type[(df_type["GENOME_TYPE"].isin(l))],
                mgm,
                f"{tag}", f"{tag}_{o}", genome_type, **kwargs
            )


def _build_start_codons(env, df, mgm, genome_type, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, str, Dict[str, Any]) -> None
    _build_start_or_stop_codons(env, df, mgm, genome_type, ["ATG", "GTG", "TTG"], **kwargs)


def _build_stop_codons(env, df, mgm, genome_type, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, str, Dict[str, Any]) -> None
    _build_start_or_stop_codons(env, df, mgm, genome_type, ["TAA", "TAG", "TGA"], **kwargs)


def build_mgm_models_from_gms2_models(env, df, mgm, **kwargs):
    # type: (Environment, pd.DataFrame, MGMModel, Dict[str, Any]) -> None

    components = get_value(
        kwargs, "components",
        {"Start Codons", "Stop Codons", "Start Context", "RBS", "Promoter"},
        valid_type=set
    )

    genome_type = get_value(kwargs, "genome_type", required=True, choices=["Archaea", "Bacteria"])

    # remove genome type from kwargs (to avoid duplicates when passing **kwargs)
    kwargs = kwargs.copy()
    kwargs.pop("genome_type")

    # start/stop codons
    if "Start Codons" in components:
        logger.info("Building start codon models")
        _build_start_codons(env, df, mgm, genome_type, **kwargs)

    if "Stop Codons" in components:
        logger.info("Building stop codon models")
        _build_stop_codons(env, df, mgm, genome_type, **kwargs)

    # Motifs
    if "RBS" in components:
        logger.info("Building RBS models")
        _build_motifs(env, df, mgm, genome_type, "RBS", **kwargs)

    if "Promoter" in components:
        logger.info("Building promoter models")
        _build_motifs(env, df, mgm, genome_type, "PROMOTER", **kwargs)

    # Start Context
    if "Start Context" in components:
        # for RBS Start Context
        logger.info("Building RBS start context models")
        _build_start_context(env, df, mgm, genome_type, "SC_RBS", **kwargs)

        logger.info("Building promoter start context models")
        _build_start_context(env, df, mgm, genome_type, "SC_PROMOTER", **kwargs)


def plot_gc_distributions(env, df):
    # type: (Environment, pd.DataFrame) -> None

    pd_figures = os_join(env["pd-work"], "figures/gc_distributions")
    mkdir_p(pd_figures)

    df = df.sort_values(["GENOME_TYPE", "Type"])
    g = seaborn.FacetGrid(df, col='GENOME_TYPE', hue="Type", col_wrap=3)
    g.map(seaborn.distplot, "GC", hist=False)
    g.set_titles("{col_name}")
    g.add_legend(loc="lower right", bbox_to_anchor=[0.833, 0.25])
    save_figure(FigureOptions(save_fig=next_name(pd_figures)))
    # plt.show()

    df_num_genomes = df.rename(columns={"GENOME_TYPE": "Group"}).groupby(
        ["Type", "Group"]
    ).size().to_frame("Number of genomes")

    # write to file
    df_num_genomes.to_csv(next_name(pd_figures, ext="csv"))


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    # Read data
    df = read_archaea_bacteria_inputs(args.pf_arc, args.pf_bac)

    # Clean up data
    logger.debug(f"Removing genetic code 4: {(df['Genetic Code'] == 4).sum()}")
    df = df[df["Genetic Code"] != 4]
    df = df[~((df["Type"] == "Archaea") & (df["GENOME_TYPE"].isin({"B", "C", "X"})))]
    df = df.convert_dtypes().copy()

    # Plot GC/Group information
    if args.plot:
        plot_gc_distributions(env, df)

    gc_feature = args.gc_feature

    # Read base MGM model
    mgm = MGMModel.init_from_file(args.pf_mgm)

    learn_from = LearnFromOptions.init_from_dict(env, args.pf_learn_from_options, vars(args))

    # Build models for archaea and bacteria
    for genome_type in ["Archaea", "Bacteria"]:

        if args.genome_type and args.genome_type != genome_type:
            # skip if specified by command line argument
            continue

        logger.info(f"Building models for {genome_type}")

        pd_figures = os_join(env["pd-work"], f"figures_{genome_type.lower()}")
        mkdir_p(pd_figures)

        build_mgm_models_from_gms2_models(
            env, df, mgm, components=args.components, genome_type=genome_type,
            plot=args.plot, gc_feature=gc_feature,
            pd_figures=pd_figures,
            learn_from=learn_from,
            cluster_by=args.cluster_by
        )

    # write new MGM file to output
    mgm.to_file(args.pf_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
