# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import pandas as pd
from typing import *
import matplotlib.pyplot as plt

# noinspection All
import seaborn

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_general import Environment
from mg_general.general import get_value
from mg_io.general import save_obj
from mg_models.building import build_mgm_motif_model_for_gc_v2
from mg_viz.colormap import ColorMap as CM

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_models.mgm_motif_model import MGMMotifModel
from mg_models.mgm_motif_model_all_gc import MGMMotifModelAllGC
from mg_models.shelf import bin_by_gc, get_consensus_sequence, read_archaea_bacteria_inputs

parser = argparse.ArgumentParser("Build MGM start models.")

parser.add_argument('--pf-input-arc', required=True, help="Input file")
parser.add_argument('--pf-input-bac', required=True, help="Input file")

parser.add_argument('--pf-output', required=True)

parser.add_argument('--pd-work', required=False, default=None, help="Path to working directory")
parser.add_argument('--pd-data', required=False, default=None, help="Path to data directory")
parser.add_argument('--pd-results', required=False, default=None, help="Path to results directory")
parser.add_argument("-l", "--log", dest="loglevel", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    help="Set the logging level", default='WARNING')

parsed_args = parser.parse_args()

# ------------------------------ #
#           Main Code            #
# ------------------------------ #

# Load environment variables
my_env = Environment(pd_data=parsed_args.pd_data,
                     pd_work=parsed_args.pd_work,
                     pd_results=parsed_args.pd_results)

# Setup logger
logging.basicConfig(level=parsed_args.loglevel)
logger = logging.getLogger("logger")  # type: logging.Logger


def plot_candidate_codons(env, df, codons, cmap=None):
    # type: (Environment, pd.DataFrame, List[str]) -> None

    fig, ax = plt.subplots()

    for c in sorted(codons):
        seaborn.regplot(df["GC"].astype(float).values, df[c].astype(float).values, label=c,
                        lowess=True, scatter_kws={"s": 5, "alpha": 0.1},
                        color=cmap[c]
                        )

    ax.set_ylim([-0.05, 1.05])
    ax.set_ylabel("Probability")
    ax.set_xlabel("GC")
    leg = ax.legend()
    for lh in leg.legendHandles:
        lh.set_alpha(1)

    plt.show()

    # bacteria vs archaea
    fig, axes = plt.subplots(1, 2, sharex="all", sharey="all")

    for t, ax in zip(["Bacteria", "Archaea"], axes.ravel()):
        df_tmp = df[df["Type"] == t]
        for c in sorted(codons):
            seaborn.regplot(df_tmp["GC"].astype(float).values, df_tmp[c].astype(float).values, label=c,
                            lowess=True, scatter_kws={"s": 5, "alpha": 0.1}, ax=ax, color=cmap[c])

        ax.set_ylim([-0.05, 1.05])
        ax.set_ylabel("Probability")
        ax.set_xlabel("GC")
        ax.set_title(t)
        leg = ax.legend()
        for lh in leg.legendHandles:
            lh.set_alpha(1)

    plt.show()

    # group
    fig, axes = plt.subplots(2, 2, sharex="all", sharey="all")

    for t, ax in zip(list("ABCD"), axes.ravel()):
        df_tmp = df[df["GENOME_TYPE"] == t]
        for c in sorted(codons):
            seaborn.regplot(df_tmp["GC"].astype(float).values, df_tmp[c].astype(float).values, label=c,
                            lowess=True, scatter_kws={"s": 5, "alpha": 0.1}, ax=ax, color=cmap[c])

        ax.set_ylim([-0.05, 1.05])
        ax.set_ylabel("Probability")
        ax.set_xlabel("GC")
        ax.set_title(t)
        leg = ax.legend()
        for lh in leg.legendHandles:
            lh.set_alpha(1)

    plt.show()


def plot_candidate_starts(env, df):
    # type: (Environment, pd.DataFrame) -> None
    plot_candidate_codons(env, df, ["ATG", "GTG", "TTG"],
                          CM.get_map("starts"))


def plot_candidate_stops(env, df):
    # type: (Environment, pd.DataFrame) -> None
    from mg_viz.colormap import ColorMap as CM
    plot_candidate_codons(env, df, ["TAA", "TAG", "TGA"],
                          CM.get_map("stops"))


def build_mgm_motif_models_for_all_gc(env, df, name, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> MGMMotifModelAllGC
    df = df[~df[name].isna()].copy()  # we only need non-NA

    bin_size = get_value(kwargs, "bin_size", 5, default_if_none=True)

    # get consensus sequences for all motifs
    df[f"CONSENSUS_{name}"] = df.apply(lambda r: get_consensus_sequence(r[name]), axis=1)

    # bin dataframes by GC
    binned_dfs = bin_by_gc(df, step=bin_size)

    # for each binned dataframe, build specific model
    list_mgm_models = list()  # type: List[Tuple[float, float, MGMMotifModel]]
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


def build_mgm_models(env, df, pf_output):
    # type: (Environment, pd.DataFrame, str) -> None

    type_model_group = {
        "Bacteria": {
            "RBS": {
                "AC", "B"
            },
            "PROMOTER": {
                "C"
            }
        },
        "Archaea": {
            "RBS": {
                "AD"
            },
            "PROMOTER": {
                "D"
            }
        }
    }

    name_to_models = dict()  # type: Dict[str, Dict[str, Dict[str, MGMMotifModelAllGC]]]
    for species_type in type_model_group.keys():
        name_to_models[species_type] = dict()  # type: Dict[str, Dict[str, MGMMotifModelAllGC]]
        for name in type_model_group[species_type].keys():
            name_to_models[species_type][name] = dict()
            for group in type_model_group[species_type][name]:
                # if group != "A" and species_type != "Archaea":
                #     continue
                name_to_models[species_type][name][group] = build_mgm_motif_models_for_all_gc(
                    env, df[(df["Type"] == species_type) & (df["GENOME_TYPE"].isin(set(group)))], name + "_MAT"
                )

    save_obj(name_to_models, pf_output)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = read_archaea_bacteria_inputs(args.pf_input_arc, args.pf_input_bac)

    # build_mgm_models(env, df, args.pf_output)
    df = df[(df["GENOME_TYPE"] != "C") | (df["GENOME_TYPE"] == "C") & (df["GC"] > 40)].copy()

    plot_candidate_starts(env, df)
    plot_candidate_stops(env, df)


if __name__ == "__main__":
    main(my_env, parsed_args)
