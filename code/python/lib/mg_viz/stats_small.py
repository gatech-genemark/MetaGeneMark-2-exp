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

def stats_small_3p_missed_vs_length(env, df_per_gene, tools, reference):
    # type: (Environment, pd.DataFrame, List[str], str) -> None

    df_with_reference = df_per_gene[~df_per_gene[f"5p-{reference}"].isna()].copy()

    # sort by reference length
    df_with_reference.sort_values(f"Length({reference})", inplace=True)
    df_with_reference.reset_index(inplace=True)


    # collect in bins
    bins = [[0, 150], [150, 300], [300, 600], [600, 900], [900, float('inf')]]

    list_entries = list()
    for t in tools + [reference]:
        df_tool = df_per_gene[~df_per_gene[f"5p-{t}"].isna()]
        for b in bins:

            # get entries where length of reference is in bin
            len_series = df_with_reference[f"Length({reference})"]
            mask = (len_series >= b[0]) & (len_series < b[1])
            df_curr = df_with_reference[mask]

            # count number of genes that are found by tool
            number_found = df_curr[f"5p-{t}"].count()
            number_reference = df_curr[f"5p-{reference}"].count()
            number_predicted = ((df_tool[f"Length({t})"] >= b[0]) & (df_tool[f"Length({t})"] < b[1])).sum()

            if t == reference:
                number_found = 0

            list_entries.append({
                "Tool": t, "Found": number_found, "Reference": number_reference,
                "Missed": (number_reference-number_found),
                "Predictions": number_predicted,
                "Exceed": max(number_predicted - number_found, 0),
                "Bin Lower": b[0],
                "Bin Upper": b[1]
            })


    df = pd.DataFrame(list_entries)
    df["Rate"] = 100 * df["Found"] / df["Reference"]

    def reorder(l_df, l_order):
        # type: (pd.DataFrame, List[str]) -> pd.DataFrame
        return l_df.reindex(l_order)

    # print(df.pivot(index="Tool", columns="Bin Upper", values="Rate").to_csv())
    reorder(df.pivot(index="Tool", columns="Bin Upper", values="Missed"), [reference] + tools).to_csv(
        next_name(env["pd-work"], ext="csv")
    )
    #
    # print(df.pivot(index="Tool", columns="Bin Upper", values="Predictions").to_csv())
    df.pivot(index="Tool", columns="Bin Upper", values="Exceed").reindex([reference] + tools).to_csv(
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

    # Number of missed per gene length bin
    stats_small_3p_missed_vs_length(env, df_per_gene, tools, reference)


def viz_stats_small_5p(env, df_per_gene, tools, list_ref):
    reference, df_tidy = _helper_join_reference_and_tidy_data(env, df_per_gene, tools, list_ref)

    # Number of found vs number of 5' error
    viz_stats_small_5p_error_vs_found(env, df_tidy, reference, tool_order=tools)
