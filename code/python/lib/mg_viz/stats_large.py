# Author: Karl Gemayel
# Created: 8/5/20, 8:25 AM

import logging
import pandas as pd
from typing import *

import seaborn
import matplotlib.pyplot as plt

from mg_general import Environment
from mg_viz.colormap import ColorMap as CM
from mg_general.general import next_name, get_value
from mg_stats.small import _helper_join_reference_and_tidy_data

logger = logging.getLogger(__name__)


def case_insensitive_match(df, col, value):
    # type: (pd.DataFrame, str, str) -> pd.Series
    return df[col].apply(lambda x: x.lower()) == value.lower()


def reorder_pivot_by_tool(df_pivoted, tool_order):
    # type: (pd.DataFrame, List[str]) -> pd.DataFrame
    return df_pivoted.reorder_levels([1, 0], 1)[
        [x.upper() for x in tool_order]].reorder_levels(
        [1, 0], 1
    ).sort_index(1, 0, sort_remaining=False)


def viz_stats_large_3p_sn_sp(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))

    reorder_pivot_by_tool(
        df_tidy.pivot(index="Genome", columns="Tool", values=["Sensitivity", "Specificity"]),
        tool_order
    ).to_csv(
        next_name(env["pd-work"], ext="csv")
    )


def stats_large_3p_predictions_vs_found(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))

    reorder_pivot_by_tool(
        df_tidy.pivot(index="Genome", columns="Tool", values=["Number of Predictions", "Number of Found"]),
        tool_order
    ).to_csv(
        next_name(env["pd-work"], ext="csv")
    )


def viz_stats_large_5p_error_vs_found(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))

    # plot against GC
    fig, axes = plt.subplots(1, 2, sharex="all", figsize=(16, 6))
    reg_kws = {"lowess": True, "scatter_kws": {"s": 2}}

    ax = axes[0]
    df_tidy["Error Rate"] = df_tidy["Number of Error"] / df_tidy["Number of Found"]     # FIXME: compute before
    col = "Error Rate"
    for t in tool_order:
        if t.lower() == reference.lower():
            continue
        df_curr = df_tidy[case_insensitive_match(df_tidy, "Tool", t)]

        seaborn.regplot(
            df_curr["Genome GC"], df_curr[col], label=t, color=CM.get_map("tools")[t.lower()],
            **reg_kws, ax=ax
        )

    ax = axes[1]
    col = "Sensitivity"
    for t in tool_order:
        df_curr = df_tidy[case_insensitive_match(df_tidy, "Tool", t)]

        if t.lower() == reference.lower():
            continue
        else:
            seaborn.regplot(
                df_curr["Genome GC"], df_curr[col], label=t, color=CM.get_map("tools")[t.lower()],
                **reg_kws, ax=ax
            )

    handles, labels = ax.get_legend_handles_labels()
    for lh in handles:
        lh.set_alpha(1)
        lh.set_sizes([8] * (len(tool_order)-1))

    leg = fig.legend(handles, labels, bbox_to_anchor=(0.5, -0.05), loc='lower center',
                     ncol=len(tool_order)-1)
    fig.tight_layout()
    fig.savefig(next_name(env["pd-work"]), bbox_extra_artists=(leg,), bbox_inches='tight')

    plt.show()


def viz_stats_large_3p(env, df_per_gene, tools, list_ref):
    reference, df_tidy = _helper_join_reference_and_tidy_data(env, df_per_gene, tools, list_ref)

    # Number of Predictions versus number of found
    stats_large_3p_predictions_vs_found(env, df_tidy, reference, tool_order=tools)

    # Number of Sensitivity and specificity
    viz_stats_large_3p_sn_sp(env, df_tidy, reference, tool_order=tools)


def viz_stats_large_5p(env, df_per_gene, tools, list_ref):
    reference, df_tidy = _helper_join_reference_and_tidy_data(env, df_per_gene, tools, list_ref)

    # Number of found vs number of 5' error
    viz_stats_large_5p_error_vs_found(env, df_tidy, reference, tool_order=tools)
