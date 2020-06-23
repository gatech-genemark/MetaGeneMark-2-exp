# Author: Karl Gemayel
# Created: 6/22/20, 3:41 PM

import logging
import argparse
import pandas as pd
from typing import *
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
from mg_general.general import get_value
from mg_models.shelf import read_archaea_bacteria_inputs

parser = argparse.ArgumentParser("DRIVER DESCRIPTION.")

parser.add_argument('--pf-bac', required=True, help="Collected GMS2 model files for bacteria")
parser.add_argument('--pf-arc', required=True, help="Collected GMS2 model files for archaea")
parser.add_argument('--pf-mgm', required=True, help="Base MGM model file")
parser.add_argument('--pf-output', required=True, help="Output MGM model file")

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
        df[c] = df[c].astype(float)
        x = df["GC"].values
        y = df[c].values
        y = get_loess(x, y)

        values_per_codon[c] = [x, y]

    df_tmp = pd.melt(df[["GC"] + codons], ["GC"], var_name="Codon", value_name="Frequency")

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
            if gc_tag == gc_max - gc_step:
                print('hi')

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
            mgm.items_by_species_and_gc[genome_type[0]][str(gc_tag)].items[c] = avg

    df_tmp = pd.DataFrame(list_entries)

    sns.scatterplot(df_tmp, "GC", "Probability", hue="Codon")


def add_start_codon_probabilities(df, mgm, **kwargs):
    # type: (pd.DataFrame, MGMModel, Dict[str, Any]) -> None
    add_codon_probabilities(df, mgm, ["ATG", "GTG", "TTG"], **kwargs)


def add_stop_codon_probabilities(df, mgm, **kwargs):
    # type: (pd.DataFrame, MGMModel, Dict[str, Any]) -> None
    add_codon_probabilities(df, mgm, ["TAA", "TAG", "TGA"], **kwargs)


def add_start_context_probabilities(df, mgm, tag, **kwargs):
    # type: (pd.DataFrame, MGMModel, str, Dict[str, Any]) -> None
    genome_type = get_value(kwargs, "genome_type", required=True, choices=["Archaea", "Bacteria"])

    gc_step = 1
    gc_min = 30
    gc_max = 71

    if genome_type == "Archaea":
        gc_step = 5
        gc_max = 66

    df = df[df["Type"] == genome_type].copy()

    example_sc = df.at[df.index[0], tag]  # type: Dict[str, List[float]]
    words = set(example_sc.keys())
    num_positions = len(next(iter(example_sc.values())))

    def empty_sc(l_words, l_num_positions):
        # type: (Iterable[str], int) -> Dict[str, List[float]]
        l_result = dict()
        for w in l_words:
            l_result[w] = [0] * l_num_positions
        return l_result

    sc_per_gc = dict()      # type: Dict[float, Dict[str, List[float]]]
    for gc_tag in range(gc_min, gc_max, gc_step):
        sc_per_gc[gc_tag] = empty_sc(words, num_positions)



    # get all words appearing in start contexts and all positions


    for p in range(num_positions):



    df.sort_values("GC", inplace=True)
    values_per_word = dict()
    for w in words:
        x = [] * len(df.index)
        y = [] * len(df.index)
        for i, idx in enumerate(df.index):
            x[i] = df.at[idx, "GC"]
            y[i] = df.at[idx, tag][w]

    values_per_codon = dict()
    for c in codons:
        df[c] = df[c].astype(float)
        x = df["GC"].values
        y = df[c].values
        y = get_loess(x, y)

        values_per_codon[c] = [x, y]

    df_tmp = pd.melt(df[["GC"] + codons], ["GC"], var_name="Codon", value_name="Frequency")

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
            if gc_tag == gc_max - gc_step:
                print('hi')

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
            mgm.items_by_species_and_gc[genome_type[0]][str(gc_tag)].items[c] = avg

    df_tmp = pd.DataFrame(list_entries)

    sns.scatterplot(df_tmp, "GC", "Probability", hue="Codon")


def build_mgm_models_from_gms2_models(df, mgm, **kwargs):
    # type: (pd.DataFrame, MGMModel, Dict[str, Any]) -> None

    components = get_value(kwargs, "components", {"Start Codons", "Stop Codons", "Start Context", "RBS", "Promoter"},
                           valid_type=set)

    # start/stop codons
    # if "Start Codons" in components:
    #     add_start_codon_probabilities(df, mgm, genome_type="Bacteria")
    #     add_start_codon_probabilities(df, mgm, genome_type="Archaea")
    #
    # if "Stop Codons" in components:
    #     add_stop_codon_probabilities(df, mgm, genome_type="Bacteria")
    #     add_stop_codon_probabilities(df, mgm, genome_type="Archaea")

    # Start Context
    if "Start Context" in components:
        add_start_context_probabilities(df, mgm, "SC_RBS_MAT", genome_type="Bacteria")
        add_start_context_probabilities(df, mgm, "SC_RBS_MAT", genome_type="Archaea")


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    df = read_archaea_bacteria_inputs(args.pf_arc, args.pf_bac)
    df = df.convert_dtypes().copy()

    mgm = MGMModel.init_from_file(args.pf_mgm)
    build_mgm_models_from_gms2_models(df, mgm)
    mgm.to_file(args.pf_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
