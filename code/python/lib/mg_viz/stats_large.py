# Author: Karl Gemayel
# Created: 8/5/20, 8:25 AM

import logging
import os

import pandas as pd
from typing import *

import seaborn
import matplotlib.pyplot as plt

from mg_general import Environment
from mg_io.general import load_obj, save_obj
from mg_viz.colormap import ColorMap as CM
from mg_general.general import next_name, get_value
from mg_stats.small import _helper_join_reference_and_tidy_data, prl_join_reference_and_tidy_data

logger = logging.getLogger(__name__)


def case_insensitive_match(df, col, value):
    # type: (pd.DataFrame, str, str) -> pd.Series
    return df[col].apply(lambda x: x.lower()) == value.lower()


def plot_gc_stats_side_by_side(env, df_tidy, columns, tool_order, reference, **kwargs):
    col_to_ylim = get_value(kwargs, "col_to_ylim", dict())

    fig, axes = plt.subplots(1, len(columns), sharex="all", figsize=(8 * len(columns), 6))
    reg_kws = {"lowess": True, "scatter_kws": {"s": 2}}

    ax = None
    for ax, col in zip(axes, columns):
        for t in tool_order:
            if t.lower() == reference.lower():
                continue
            df_curr = df_tidy[case_insensitive_match(df_tidy, "Tool", t)]

            seaborn.regplot(
                df_curr["Genome GC"], df_curr[col], label=t, color=CM.get_map("tools")[t.lower()],
                **reg_kws, ax=ax
            )

            if col in col_to_ylim:
                ax.set_ylim(*col_to_ylim[col])

    if ax is not None:
        fig.subplots_adjust(bottom=0.3)

        handles, labels = ax.get_legend_handles_labels()
        for lh in handles:
            lh.set_alpha(1)
            lh.set_sizes([8] * (len(tool_order)))
        leg = fig.legend(handles, labels, bbox_to_anchor=(0.5, 0.2), loc='upper center',
                         ncol=len(tool_order), bbox_transform=fig.transFigure)

        fig.savefig(next_name(env["pd-work"]), bbox_inches="tight")

    plt.show()


def reorder_pivot_by_tool(df_pivoted, tool_order):
    # type: (pd.DataFrame, List[str]) -> pd.DataFrame
    return df_pivoted.reorder_levels([1, 0], 1)[
        [x.upper() for x in tool_order]].reorder_levels(
        [1, 0], 1
    ).sort_index(1, 0, sort_remaining=False)


def viz_stats_large_3p_sn_sp(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))

    plot_gc_stats_side_by_side(
        env, df_tidy, ["Number of Found", "Sensitivity"], tool_order, reference,
        col_to_ylim={"Specificity": (0.5, 1), "Sensitivity": (0.5, 1)}
    )


def stats_large_3p_predictions_vs_found(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))

    plot_gc_stats_side_by_side(
        env, df_tidy, ["Number of Predictions", "Number of Found", "Specificity"], tool_order, reference,
        col_to_ylim={"Specificity": (0.5, 1), "Sensitivity": (0.5, 1)}
    )


def viz_stats_large_5p_error_vs_sensitivity(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))

    df_tidy["Error Rate"] = df_tidy["Number of Error"] / df_tidy["Number of Found"]  # FIXME: compute before
    plot_gc_stats_side_by_side(
        env, df_tidy, ["Error Rate", "Sensitivity"], tool_order, reference,
        col_to_ylim={"Specificity": (0.5, 1), "Sensitivity": (0.5, 1), "Error Rate": (0, 0.5)}
    )


def viz_stats_large_3p(env, df_per_gene, tools, list_ref, **kwargs):
    pf_checkpoint = get_value(kwargs, "pf_checkpoint", None)
    if not pf_checkpoint or not os.path.isfile(pf_checkpoint):
        reference, df_tidy = prl_join_reference_and_tidy_data(env, df_per_gene, tools, list_ref)
        if pf_checkpoint:
            save_obj([reference, df_tidy], pf_checkpoint)
    else:
        reference, df_tidy = load_obj(pf_checkpoint)

    # Number of Predictions versus number of found
    stats_large_3p_predictions_vs_found(env, df_tidy, reference, tool_order=tools)

    # Number of Sensitivity and specificity
    viz_stats_large_3p_sn_sp(env, df_tidy, reference, tool_order=tools)


def viz_stats_large_5p(env, df_per_gene, tools, list_ref, **kwargs):
    pf_checkpoint = get_value(kwargs, "pf_checkpoint", None)
    if not pf_checkpoint or not os.path.isfile(pf_checkpoint):
        reference, df_tidy = prl_join_reference_and_tidy_data(env, df_per_gene, tools, list_ref)
        if pf_checkpoint:
            save_obj([reference, df_tidy], pf_checkpoint)
    else:
        reference, df_tidy = load_obj(pf_checkpoint)

    # Number of found vs number of 5' error
    viz_stats_large_5p_error_vs_sensitivity(env, df_tidy, reference, tool_order=tools)
