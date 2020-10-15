# Author: Karl Gemayel
# Created: 6/23/20, 8:16 PM

import logging
from typing import *
from collections import Counter

import pandas as pd
import numpy as np

from mg_general import Environment
from mg_general.general import get_value
from mg_models.gms2_noncoding import GMS2Noncoding
from mg_models.mgm_motif_model import MGMMotifModel
from mg_models.mgm_motif_model_v2 import MGMMotifModelV2
from mg_models.shelf import create_numpy_for_column_with_extended_motif, get_position_distributions_by_shift, \
    create_numpy_for_column_by_clustering
from mg_viz.mgm_motif_model import MGMMotifModelVisualizer
from mg_viz.mgm_motif_model_v2 import MGMMotifModelVisualizerV2

log = logging.getLogger(__name__)


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
        mod = GMS2Noncoding(df.at[idx, "Mod"].items["NON_MAT"])
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
    plot = get_value(kwargs, "plot", False, valid_type=bool)
    gc_feature = get_value(kwargs, "gc_feature", "GC", valid_type=str)
    cluster_by = get_value(kwargs, "cluster_by", "msa")


    # filter out consensus sequences that don't appear very frequently
    df = df[df.groupby(f"CONSENSUS_{col}")[f"CONSENSUS_{col}"].transform(len) > min_consensus_occurence]

    if len(df) <= 1:
        return

    original_width = len(df.iloc[0][f"CONSENSUS_{col}"])

    # run alignment of consensus sequences and get back shifts
    collect = dict()
    if cluster_by == "msa":
        array, update_shifts = create_numpy_for_column_with_extended_motif(env, df, col,
                                                                       collect)
    else:
        print("Heuristic")
        array, update_shifts = create_numpy_for_column_by_clustering(env, df, col,
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
                             avg_gc=df[gc_feature].mean(),
                             avg_bgd=avg_bgd)

    if plot:
        pd_figures = get_value(kwargs, "pd_figures", env["pd-work"])
        MGMMotifModelVisualizerV2.visualize(mgm_mm, title=title, msa_t=collect["msa_t"],
                                        raw_motif_data=raw_array_per_shift, cluster_by=cluster_by,
                                            pd_figures=pd_figures)
    return mgm_mm