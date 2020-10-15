# Author: Karl Gemayel
# Created: 8/5/20, 8:25 AM

import logging
import math
import os
from textwrap import wrap

import pandas as pd
from typing import *

import seaborn
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from mg_general import Environment
from mg_io.general import load_obj, save_obj
from mg_viz.colormap import ColorMap as CM
from mg_general.general import next_name, get_value
from mg_stats.small import _helper_join_reference_and_tidy_data, prl_join_reference_and_tidy_data
from mg_viz.general import set_size
from mg_viz.shelf import number_formatter, update_tool_names_to_full

logger = logging.getLogger(__name__)


def case_insensitive_match(df, col, value):
    # type: (pd.DataFrame, str, str) -> pd.Series
    return df[col].apply(lambda x: x.lower()) == value.lower()


def plot_gc_stats_side_by_side(env, df_tidy, columns, tool_order, reference, **kwargs):
    col_to_ylim = get_value(kwargs, "col_to_ylim", dict())
    col_wrap = get_value(kwargs, "col_wrap", len(columns))
    num_rows = math.ceil(len(columns) / float(col_wrap))
    wrap_val = get_value(kwargs, "wrap_val", None)
    figsize = get_value(kwargs, "figsize", (8 * col_wrap, 6 * num_rows))
    col_x = get_value(kwargs, "col_x", "Genome GC")
    col_x_text = get_value(kwargs, "col_x", "GC")
    legend_cols = get_value(kwargs, "legend_cols", len(tool_order))
    legend_pos = get_value(kwargs, "legend_pos", "bottom")
    fig, axes = plt.subplots(num_rows, col_wrap, figsize=figsize)
    reg_kws = {"lowess": True, "scatter_kws": {"s": 0.1, "alpha": 0.3},
               "line_kws": {"linewidth": 1}}
    from collections import abc

    axes_unr = axes
    if not isinstance(axes, abc.Iterable):
        axes = [axes]
    else:
        axes = axes.ravel()

    ax = None
    i = j = 0
    fontsize="small"
    for ax, col in zip(axes, columns):
        for t in tool_order:
            if t.lower() == reference.lower():
                continue
            df_curr = df_tidy[case_insensitive_match(df_tidy, "Tool", t)]

            seaborn.regplot(
                df_curr[col_x], df_curr[col], label=t, color=CM.get_map("tools")[t.lower()],
                **reg_kws, ax=ax
            )

            if col in col_to_ylim:
                ax.set_ylim(*col_to_ylim[col])

            if max(df_curr[col]) > 2000:
                ax.yaxis.set_major_formatter(FuncFormatter(number_formatter))

            if i != num_rows - 1:
                ax.set_xlabel("")
            else:
                ax.set_xlabel(col_x_text, fontsize=fontsize)

            if wrap_val:
                col_text = "\n".join(wrap(col, wrap_val, break_long_words=False))
            else:
                col_text = col
            ax.set_ylabel(col_text, wrap=True, fontsize=fontsize)
            ax.tick_params(labelsize=fontsize, length=2)

        j += 1
        if j == col_wrap:
            i += 1
            j = 0

    if ax is not None:
        if legend_pos == "bottom":
            fig.subplots_adjust(bottom=0.2)
        else:
            fig.subplots_adjust(right=0.8)
        handles, labels = ax.get_legend_handles_labels()

        # labels = [{
        #     "mgm": "MGM",
        #     "mgm2": "MGM2",
        #     "mga": "MGA",
        #     "mprodigal": "MProdigal",
        #     "fgs": "FGS",
        #     "gms2": "GMS2",
        #     "prodigal": "Prodigal"
        # }[l.lower()] for l in labels]
        labels = update_tool_names_to_full(labels)

        if legend_pos == "bottom" or True:
            leg = fig.legend(handles, labels, bbox_to_anchor=(0.5, 0.1), loc='upper center', ncol=legend_cols,
                             bbox_transform=fig.transFigure, frameon=False,
                             fontsize="xx-small")
        else:
            leg = fig.legend(handles, labels, bbox_to_anchor=(1.05, 0.5), loc='center left',
                             frameon=False,
                             fontsize=18)
        for lh in leg.legendHandles:
            lh.set_alpha(1)
            lh.set_sizes([18] * (len(tool_order)))

        if num_rows > 1:
            for i in range(col_wrap):
                fig.align_ylabels(axes_unr[:,i])

        if legend_pos == "bottom" or True:
            if num_rows == 1:
                fig.tight_layout(rect=[0,0.05,1,1])
            else:
                fig.tight_layout(rect=[0,0.1,1,1])
        # else:
        #     fig.tight_layout(rect=[0, 0, 1, 1])
        fig.savefig(next_name(env["pd-work"]), bbox_extra_artists=(leg,)) #bbox_inches='tight'

    plt.show()


def reorder_pivot_by_tool(df_pivoted, tool_order):
    # type: (pd.DataFrame, List[str]) -> pd.DataFrame
    return df_pivoted.reorder_levels([1, 0], 1)[
        [x.upper() for x in tool_order]].reorder_levels(
        [1, 0], 1
    ).sort_index(1, 0, sort_remaining=False)


def stats_large_3p_reference(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))

    # if reference not in tool_order:
    #     tool_order = [reference] + tool_order

    # number of genes per clade
    df_grouped = df_tidy.groupby(["Clade", "Tool"], as_index=False).sum()
    df_grouped["Sensitivity"] = df_grouped["Number of Found"] / df_grouped["Number in Reference"]
    df_grouped["Specificity"] = df_grouped["Number of Found"] / df_grouped["Number of Predictions"]

    # df_pivoted = reorder_pivot_by_tool(
    #     df_grouped.pivot(index=["Clade", "Number in Reference"], columns="Tool", values=["Sensitivity", "Specificity"]), tool_order
    # )

    df_pivoted = reorder_pivot_by_tool(df_grouped.pivot_table(
        index=["Clade", "Number in Reference"], columns="Tool",
        values=["Sensitivity", "Specificity"]).reset_index(
        level=1),
        tool_order)

    df_pivoted.to_csv(
        next_name(env["pd-work"], ext="csv")
    )


def stats_large_5p_overall(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))

    # if reference not in tool_order:
    #     tool_order = [reference] + tool_order

    # number of genes per clade
    df_grouped = df_tidy.groupby(["Clade", "Tool"], as_index=False).sum()
    df_grouped["Error Rate"] = df_grouped["Number of Error"] / df_grouped["Number of Found"]

    df_pivoted = reorder_pivot_by_tool(df_grouped.pivot_table(
        index=["Clade", "Number in Reference"], columns="Tool",
        values=["Error Rate"]).reset_index(
        level=1),
        tool_order)

    df_pivoted.to_csv(
        next_name(env["pd-work"], ext="csv")
    )


def viz_stats_large_3p_sn_sp(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))

    df_tidy["3' FN Error Rate"] = 1- df_tidy["Sensitivity"]
    df_tidy["3' FP Error Rate"] = 1 - df_tidy["Specificity"]

    plot_gc_stats_side_by_side(
        env, df_tidy, ["Sensitivity", "Specificity", "Number of Found", "Number of Predictions"],
        tool_order, reference, col_wrap=2, wrap_val=10, figsize=set_size(433.62001, subplots=(2,2), legend="bottom"),
        col_to_ylim={"Specificity": (0.5, 1), "Sensitivity": (0.5, 1)}
    )

    plot_gc_stats_side_by_side(
        env, df_tidy, ["3' FN Error Rate", "3' FP Error Rate"],
        tool_order, reference, col_wrap=2, wrap_val=10, figsize=set_size(433.62001, subplots=(1, 2), legend="bottom"),
        col_to_ylim={"3' FN Error Rate": (0, 0.2), "3' FP Error Rate": (0, 0.2)},
        legend_cols = math.ceil(len(tool_order)), legend_pos="right"
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

    df_tidy["Gene Start Error Rate"] = df_tidy["Number of Error"] / df_tidy["Number of Found"]  # FIXME: compute before
    df_tidy["3' FN Error Rate"] = 1 - df_tidy["Sensitivity"]
    plot_gc_stats_side_by_side(
        env, df_tidy, ["Gene Start Error Rate", "3' FN Error Rate"], tool_order, reference,
        col_wrap=2, wrap_val=10, figsize=set_size("thesis", subplots=(1, 2), legend="bottom"),
        col_to_ylim={"Specificity": (0.5, 1), "Gene Start Error Rate": (0, 0.3), "3' FN Error Rate": (0, 0.3)}
    )

    df_tidy["Gene 5' Error Rate"] = df_tidy["Gene Start Error Rate"]
    plot_gc_stats_side_by_side(
        env, df_tidy, ["Gene 5' Error Rate", "3' FN Error Rate"], tool_order, reference,
        col_wrap=2, wrap_val=10, figsize=set_size("thesis", subplots=(1, 2), legend="bottom"),
        col_to_ylim={"Specificity": (0.5, 1), "Gene 5' Error Rate": (0, 0.3), "3' FN Error Rate": (0, 0.15)}
    )

    print(df_tidy.groupby("Tool", as_index=False).mean().to_csv(index=False))
    print(df_tidy.groupby("Tool", as_index=False).sum().to_csv(index=False))

def viz_stats_large_5p_error_vs_gc_by_clade(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))

    df_tidy["Error Rate"] = df_tidy["Number of Error"] / df_tidy["Number of Found"]

    clades_sorted = sorted(df_tidy["Clade"].unique())
    num_clades = len(clades_sorted)
    num_rows = 2
    subplots=(num_rows, math.ceil(num_clades/ float(num_rows)))
    figsize = set_size("thesis", subplots=subplots,legend="bottom", titles=True)
    col_x = "Genome GC"
    col_x_text = "GC"

    fig, axes = plt.subplots(subplots[0], subplots[1], figsize=figsize, sharex="all", sharey="all")
    reg_kws = {"lowess": True, "scatter_kws": {"s": 0.1, "alpha": 0.3},
               "line_kws": {"linewidth": 1}}
    from collections import abc

    axes_unr = axes
    if not isinstance(axes, abc.Iterable):
        axes = [axes]
    else:
        axes = axes.ravel()

    ax = None
    fontsize = "xx-small"
    counter = 0
    for ax, col in zip(axes, clades_sorted):
        for t in tool_order:
            if t.lower() == reference.lower():
                continue
            df_curr = df_tidy[case_insensitive_match(df_tidy, "Tool", t)]
            df_curr = df_curr[df_curr["Clade"] == col]

            seaborn.regplot(
                df_curr[col_x], df_curr["Error Rate"], label=t, color=CM.get_map("tools")[t.lower()],
                **reg_kws, ax=ax
            )

            # if col in col_to_ylim:
            #     ax.set_ylim(*col_to_ylim[col])

        if max(df_curr["Error Rate"]) > 2000:
            ax.yaxis.set_major_formatter(FuncFormatter(number_formatter))

        ax.set_xlabel(col_x_text, fontsize=fontsize)
        ax.set_title(col, fontsize=fontsize)

        ax.set_ylabel("Error Rate", wrap=True, fontsize=fontsize)
        ax.tick_params(labelsize=fontsize, length=2)
        if counter == 0:
            ax.set_ylabel("Error Rate", wrap=True, fontsize=fontsize)
        else:
            ax.set_ylabel("")


    if ax is not None:
        fig.subplots_adjust(bottom=0.2)
        handles, labels = ax.get_legend_handles_labels()

        # labels = [{
        #               "mgm": "MGM",
        #               "mgm2": "MGM2",
        #               "mga": "MGA",
        #               "mprodigal": "MProdigal",
        #               "fgs": "FGS",
        #               "gms2": "GMS2",
        #               "prodigal": "Prodigal"
        #           }[l.lower()] for l in labels]
        labels = update_tool_names_to_full(labels)

        leg = fig.legend(handles, labels, bbox_to_anchor=(0.5, 0.1), loc='upper center', ncol=len(tool_order),
                         bbox_transform=fig.transFigure, frameon=False,
                         fontsize=fontsize)
        for lh in leg.legendHandles:
            lh.set_alpha(1)
            lh.set_sizes([18] * (len(tool_order)))

        # if num_rows > 1:
        #     for i in range():
        #         fig.align_ylabels(axes_unr[:, i])

        if num_rows == 1:
            fig.tight_layout(rect=[0, 0.05, 1, ])
        else:
            fig.tight_layout(rect=[0, 0.1, 1, 1])
        fig.savefig(next_name(env["pd-work"]), bbox_extra_artists=(leg,))  # bbox_inches='tight'

    plt.show()

def viz_stats_large_3p(env, df_per_gene, tools, list_ref, **kwargs):
    pf_checkpoint = get_value(kwargs, "pf_checkpoint", None)
    if not pf_checkpoint or not os.path.isfile(pf_checkpoint):
        reference, df_tidy = prl_join_reference_and_tidy_data(env, df_per_gene, tools, list_ref)
        if pf_checkpoint:
            save_obj([reference, df_tidy], pf_checkpoint)
    else:
        reference, df_tidy = load_obj(pf_checkpoint)

    # Reference stats
    df_tidy.loc[df_tidy["Tool"] == "MGM2_AUTO", "Tool"] = "MGM2"
    reference = reference.replace("MGM2_AUTO", "MGM2")
    tools = tools.copy()
    for i in range(len(tools)):
        if tools[i].upper() == "MGM2_AUTO":
            tools[i] = "MGM2"

    stats_large_3p_reference(env, df_tidy, reference, tool_order=tools)

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

    df_tidy.loc[df_tidy["Tool"] == "MGM2_AUTO", "Tool"] = "MGM2"
    reference = reference.replace("MGM2_AUTO", "MGM2")
    tools = tools.copy()
    for i in range(len(tools)):
        if tools[i].upper() == "MGM2_AUTO":
            tools[i] = "MGM2"

    stats_large_5p_overall(env, df_tidy, reference, tool_order=tools)

    # Number of found vs number of 5' error
    viz_stats_large_5p_error_vs_sensitivity(env, df_tidy, reference, tool_order=tools)

    viz_stats_large_5p_error_vs_gc_by_clade(env, df_tidy, reference, tool_order=tools)
