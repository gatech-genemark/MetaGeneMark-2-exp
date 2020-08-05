# Author: Karl Gemayel
# Created: 8/5/20, 8:25 AM

import logging
import pandas as pd
from typing import *

from mg_general import Environment
from mg_general.general import next_name, get_value
from mg_stats.small import _helper_join_reference_and_tidy_data

logger = logging.getLogger(__name__)


def reorder_pivot_by_tool(df_pivoted, tool_order):
    # type: (pd.DataFrame, List[str]) -> pd.DataFrame
    return df_pivoted.reorder_levels([1, 0], 1)[
        [x.upper() for x in tool_order]].reorder_levels(
        [1, 0], 1
    ).sort_index(1, 0, sort_remaining=False)


def viz_stats_small_3p_sn_sp(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))

    reorder_pivot_by_tool(
        df_tidy.pivot(index="Genome", columns="Tool", values=["Sensitivity", "Specificity"]),
        tool_order
    ).to_csv(
        next_name(env["pd-work"], ext="csv")
    )


def stats_small_3p_predictions_vs_found(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))

    reorder_pivot_by_tool(
        df_tidy.pivot(index="Genome", columns="Tool", values=["Number of Predictions", "Number of Found"]),
        tool_order
    ).to_csv(
        next_name(env["pd-work"], ext="csv")
    )


def viz_stats_small_5p_error_vs_found(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))

    reorder_pivot_by_tool(
        df_tidy.pivot(index="Genome", columns="Tool", values=["Number of Found", "Number of Error"]),
        tool_order
    ).to_csv(next_name(env["pd-work"], ext="csv"))


def viz_stats_small_3p(env, df_per_gene, tools, list_ref):
    reference, df_tidy = _helper_join_reference_and_tidy_data(env, df_per_gene, tools, list_ref)

    # Number of Predictions versus number of found
    stats_small_3p_predictions_vs_found(env, df_tidy, reference, tool_order=tools)

    # Number of Sensitivity and specificity
    viz_stats_small_3p_sn_sp(env, df_tidy, reference, tool_order=tools)


def viz_stats_small_5p(env, df_per_gene, tools, list_ref):
    reference, df_tidy = _helper_join_reference_and_tidy_data(env, df_per_gene, tools, list_ref)

    # Number of found vs number of 5' error
    viz_stats_small_5p_error_vs_found(env, df_tidy, reference, tool_order=tools)
