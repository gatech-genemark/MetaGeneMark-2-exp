# Karl Gemayel
# Georgia Institute of Technology
#
# Created: 3/17/20
import logging
import argparse
import numpy as np
import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from mg_general import Environment
from mg_general.general import get_value
from mg_io.general import load_obj, save_obj

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_models.gms2_noncoding import GMS2Noncoding
from mg_models.mgm_motif_model import MGMMotifModel
from mg_models.mgm_motif_model_all_gc import MGMMotifModelAllGC
from mg_models.mgm_motif_model_v2 import MGMMotifModelV2
from mg_models.shelf import bin_by_gc, get_consensus_sequence, get_position_distributions_by_shift, \
    create_numpy_for_column_with_extended_motif

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


def fix_genome_type(df):
    # type: (pd.DataFrame) -> None
    df["GENOME_TYPE"] = df["GENOME_TYPE"].apply(lambda x: x.strip().split("-")[1].upper())
    df.loc[df["GENOME_TYPE"] == "D2", "GENOME_TYPE"] = "D"


def read_archaea_bacteria_inputs(pf_arc, pf_bac):
    # type: (str, str) -> pd.DataFrame
    df_bac = load_obj(pf_bac)  # type: pd.DataFrame
    df_arc = load_obj(pf_arc)  # type: pd.DataFrame
    df_bac["Type"] = "Bacteria"
    df_arc["Type"] = "Archaea"

    df = pd.concat([df_bac, df_arc], sort=False)
    df.reset_index(inplace=True)
    fix_genome_type(df)

    return df


def mat_to_dict(mat):
    # type: (np.ndarray) -> Dict[str, List[float]]

    index_to_letter = {
        i: x for i, x in enumerate(list("ACGT"))
    }

    result = dict()

    P, L = mat.shape
    for l in range(L):
        l_str = index_to_letter[l]
        result[l_str] = list(mat[:, l])

    return result


def get_average_zero_order_noncoding(df):
    # type: (pd.DataFrame) -> np.ndarray
    list_arr = list()
    for idx in df.index:
        mod = GMS2Noncoding(df.at[idx, "NON_MAT"])
        list_arr.append(mod.pwm_to_array(0))

    avg = np.mean(list_arr, 0)

    return avg


def build_mgm_motif_model_for_gc(env, df, col, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> Union[MGMMotifModel, None]

    min_consensus_occurence = get_value(kwargs, "min_consensus_occurence", 5)
    title = get_value(kwargs, "title", "")

    # filter out consensus sequences that don't appear very frequently
    df = df[df.groupby(f"CONSENSUS_{col}")[f"CONSENSUS_{col}"].transform(len) > min_consensus_occurence]

    if len(df) <= 1:
        return

    # run alignment of consensus sequences and get back shifts
    collect = dict()
    array, update_shifts = create_numpy_for_column_with_extended_motif(env, df, col,
                                                                       collect)

    # two widths:
    original_width = len(df.iloc[0][f"CONSENSUS_{col}"])
    extended_width = array.shape[1]

    # array has shape: N x P x L
    # N: Number of motifs
    # P: Width of extended motif
    # L: number of letters
    N, P, L = array.shape

    extended_motif = np.sum(array, 0)
    extended_motif = np.divide(
        extended_motif,
        extended_motif.sum(1).reshape(P, 1)
    )

    extended_motif_dict = mat_to_dict(extended_motif)

    # get prior probabilities on shift position
    counter = Counter(update_shifts)
    total = sum(counter.values())
    to_add = sorted(set(range(max(counter.keys()))).difference(counter.keys()))
    normalized = [[x, 100 * counter[x] / total] for x in counter] + [[x, 0] for x in to_add]
    normalized = np.array(normalized)
    shifts_dict = {normalized[x, 0]: normalized[x, 1] for x in range(normalized.shape[0])}

    # get position distributions
    col_pos = col.replace("_MAT", "_POS_DISTR")
    shift_to_pos_dist = get_position_distributions_by_shift(df, col_pos, update_shifts)

    position_distributions_by_shift = dict()  # type: Dict[int, Dict[int, float]]
    for s in sorted(shift_to_pos_dist.keys()):
        list_pos_dist = shift_to_pos_dist[s]

        # average positions
        values = dict()
        for l in list_pos_dist:
            try:
                for i in l.keys():
                    if i not in values.keys():
                        values[i] = list()
                    values[i].append(l[i])
            except Exception:
                continue
        for i in values.keys():
            values[i] = np.mean(values[i])

        total = sum(values.values())
        for i in values.keys():
            values[i] /= total

        x = sorted(values.keys())
        y = [values[a] for a in x]

        position_distributions_by_shift[s] = {
            a: b for a, b in zip(x, y)
        }

    # compile into single model

    avg_bgd = get_average_zero_order_noncoding(df)

    mgm_mm = MGMMotifModel(shifts_dict, extended_motif_dict, original_width, position_distributions_by_shift,
                           avg_gc=df["GC"].mean(),
                           avg_bgd=avg_bgd)

    MGMMotifModelVisualizer.visualize(mgm_mm, title=title, msa_t=collect["msa_t"],
                                      raw_motif_data=[array, update_shifts])
    return mgm_mm


def build_mgm_motif_model_for_gc_v2(env, df, col, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> Union[MGMMotifModelV2, None]

    min_consensus_occurence = get_value(kwargs, "min_consensus_occurence", 5)
    title = get_value(kwargs, "title", "")

    # filter out consensus sequences that don't appear very frequently
    df = df[df.groupby(f"CONSENSUS_{col}")[f"CONSENSUS_{col}"].transform(len) > min_consensus_occurence]

    if len(df) <= 1:
        return

    original_width = len(df.iloc[0][f"CONSENSUS_{col}"])

    # run alignment of consensus sequences and get back shifts
    collect = dict()
    array, update_shifts = create_numpy_for_column_with_extended_motif(env, df, col,
                                                                       collect)

    # Separate motifs per shift
    unique_shifts = sorted(set(update_shifts))
    array_per_shift = {
        s: list() for s in unique_shifts
    }

    for i in range(len(update_shifts)):
        shift = update_shifts[i]
        array_per_shift[shift].append(array[i, shift:shift + original_width, :])

    raw_array_per_shift = {
        x: np.array(array_per_shift[x]) for x in array_per_shift.keys()
    }

    for s in unique_shifts:
        array_per_shift[s] = np.sum(array_per_shift[s], axis=0)
        array_per_shift[s] = np.divide(
            array_per_shift[s],
            array_per_shift[s].sum(1).reshape(original_width, 1)
        )
        array_per_shift[s] = mat_to_dict(array_per_shift[s])

    # get prior probabilities on shift position
    counter = Counter(update_shifts)
    total = sum(counter.values())
    to_add = sorted(set(range(max(counter.keys()))).difference(counter.keys()))
    normalized = [[x, 100 * counter[x] / total] for x in counter] + [[x, 0] for x in to_add]
    normalized = np.array(normalized)
    shifts_dict = {normalized[x, 0]: normalized[x, 1] for x in range(normalized.shape[0])}

    # get position distributions
    col_pos = col.replace("_MAT", "_POS_DISTR")
    shift_to_pos_dist = get_position_distributions_by_shift(df, col_pos, update_shifts)

    position_distributions_by_shift = dict()  # type: Dict[int, Dict[int, float]]
    for s in sorted(shift_to_pos_dist.keys()):
        list_pos_dist = shift_to_pos_dist[s]

        # average positions
        values = dict()
        for l in list_pos_dist:
            try:
                for i in l.keys():
                    if i not in values.keys():
                        values[i] = list()
                    values[i].append(l[i])
            except Exception:
                continue
        for i in values.keys():
            values[i] = np.mean(values[i])

        total = sum(values.values())
        for i in values.keys():
            values[i] /= total

        x = sorted(values.keys())
        y = [values[a] for a in x]

        position_distributions_by_shift[s] = {
            a: b for a, b in zip(x, y)
        }

    # compile into single model

    avg_bgd = get_average_zero_order_noncoding(df)

    mgm_mm = MGMMotifModelV2(shifts_dict, array_per_shift, original_width, position_distributions_by_shift,
                             avg_gc=df["GC"].mean(),
                             avg_bgd=avg_bgd)

    MGMMotifModelVisualizerV2.visualize(mgm_mm, title=title, msa_t=collect["msa_t"],
                                        raw_motif_data=raw_array_per_shift)
    return mgm_mm



def plot_candidate_codons(env, df, codons, cmap=None):
    # type: (Environment, pd.DataFrame, List[str]) -> None

    fig, ax = plt.subplots()
    from sbsp_viz.colormap import ColorMap as CM

    for c in sorted(codons):
        seaborn.regplot(df["GC"].astype(float).values, df[c].astype(float).values, label=c,
                        lowess=True, scatter_kws={"s": 5, "alpha": 0.1},
                        color=cmap[c]
                        )

    ax.set_ylim([-0.05,1.05])
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
    from sbsp_viz.colormap import ColorMap as CM
    plot_candidate_codons(env, df, ["ATG", "GTG", "TTG"],
                          CM.get_map("starts"))


def plot_candidate_stops(env, df):
    # type: (Environment, pd.DataFrame) -> None
    from sbsp_viz.colormap import ColorMap as CM
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
