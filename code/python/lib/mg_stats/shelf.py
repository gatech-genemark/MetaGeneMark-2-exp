# Author: Karl Gemayel
# Created: 7/23/20, 8:43 PM

import logging
from functools import reduce

import pandas as pd
from typing import *

from mg_general import Environment
from mg_general.general import fix_names

log = logging.getLogger(__name__)


def all_columns_equal(df, columns=None):
    # type: (pd.DataFrame, List[str]) -> pd.Series
    """Return True/False series shown rows where all columns have the same value"""

    if columns is None:
        columns = df.columns.values

    # create condition list
    conditions = list()
    for i in range(1, len(columns)):
        conditions.append(f"(df[{columns[i-1]}] == df[{columns[i]}])")

    return eval(" & ".join(conditions))


def create_joint_reference_from_list(df, list_reference):
    # type: (pd.DataFrame, List[str]) -> str

    reference = "=".join(list_reference)

    reference_rows = all_columns_equal(df, [f"'5p-{r}'" for r in list_reference])
    reference_values = df.loc[reference_rows, f"5p-{list_reference[0]}"]
    df.loc[reference_rows, f"5p-{reference}"] = reference_values

    reference_rows = all_columns_equal(df, [f"'3p-{r}'" for r in list_reference])
    reference_values = df.loc[reference_rows, f"3p-{list_reference[0]}"]
    df.loc[reference_rows, f"3p-{reference}"] = reference_values

    list_partial = [f"'Partial3p-{r}'" for r in list_reference if f"Partial3p-{r}" in df.columns.values]
    if len(list_partial) > 0:

        reference_rows = all_columns_equal(df, list_partial)
        reference_values = df.loc[reference_rows, f"Partial3p-{list_reference[0]}"]
        df.loc[reference_rows, f"Partial3p-{reference}"] = reference_values



    list_partial = [f"'Partial5p-{r}'" for r in list_reference if f"Partial5p-{r}" in df.columns.values]
    if len(list_partial) > 0:

        reference_rows = all_columns_equal(df, list_partial)
        reference_values = df.loc[reference_rows, f"Partial5p-{list_reference[0]}"]
        df.loc[reference_rows, f"Partial5p-{reference}"] = reference_values


    return reference


def update_dataframe_with_stats(df, tools, reference):
    # type: (pd.DataFrame, List[str], List[str]) -> pd.DataFrame
    for t in tools:
        tag_5p = f"5p:Match({t}={reference})"
        tag_3p = f"3p:Match({t}={reference})"

        # match by 5prime end
        df[tag_5p] = df[f"5p-{t}"] == df[f"5p-{reference}"]

        # all tools have a prediction
        df[tag_3p] = df[[f"5p-{t}", f"5p-{reference}"]].notnull().all(axis=1)

        df[f"Length({t})"] = df.apply(
            lambda r: abs(r[f"3p-{t}"] - r[f"5p-{t}"]) + 1,
            axis=1
        )

    df[f"Length({reference})"] = df.apply(
        lambda r: abs(r[f"3p-{reference}"] - r[f"5p-{reference}"]) + 1,
        axis=1
    )

    # remove short

    if f"Partial5p-{reference}" in df.columns:
        before = len(df)
        df = df[~((df[f"Length({reference})"] < 90) & (
                (df[f"Partial5p-{reference}"]) | (df[f"Partial3p-{reference}"]))
                  )]
        after = len(df)
        log.info(f"Filtered {before-after} short partial genes")
    return df


def tidy_genome_level(env, df):
    # type: (Environment, pd.DataFrame) -> pd.DataFrame
    """Creates a tidy dataframe for all metrics"""
    values_to_melt = ["Match", "Number of Error", "Number of Found", "Number of Match", "Number of Predictions",
                      "Number of IC5p Match", "Number of IC5p Found", "Number of IC3p Match", "Number of IC3p Found",
                      "Number of Comp Match", "Number of Comp Found", "Precision", "Recall", "WR", "Number of Missed",
                      "Sensitivity", "Specificity", "Error Rate",
                      "IC3p Match", "IC5p Match", "Comp Match"]
    df_total = list()

    list_index = [x for x in ["Genome", "Clade", "Chunk Size", "Genome GC", "Number in Reference"] if x in df.columns]
    for v in values_to_melt:
        value_vars = [x for x in df.columns if v == x.split("(")[0].strip()]
        if len(value_vars) == 0:
            continue
        df_curr = pd.melt(df, id_vars=list_index,
                          value_vars=value_vars,
                          var_name="Combination", value_name=v)
        df_curr["Tool"] = df_curr["Combination"].apply(lambda x: x.split("(")[1].split(",")[0].upper())
        df_total.append(df_curr)

    return reduce(lambda df1, df2: pd.merge(df1, df2, on=list_index + ["Tool"],
                                            how="outer"), df_total)


def _helper_df_joint_reference(df, reference):
    # create joint references when possible
    if len(reference) > 1:
        reference = create_joint_reference_from_list(df, reference)
    else:
        reference = reference[0]

    return reference


def check_tools_and_reference_lists(df, tools, ref_5p, ref_3p):
    # type: (pd.DataFrame, Union[List[str], None], List[str], List[str]) -> List[str]
    """Verifies references exist in dataframe, and returns cleaned up tools list (without references in it)"""
    # get tools list
    # If not provided, extract from df
    # Make sure it doesn't contain any references
    all_tools = sorted(set([x.split("-")[1] for x in df.columns if "5p-" in x]))

    # check that references exist
    for list_ref in [ref_5p, ref_3p]:
        for ref in list_ref:
            if ref not in all_tools:
                raise ValueError(f"Unknown reference {ref}")

    if tools is None:
        tools = all_tools
        tools = sorted(set(tools).difference({*ref_5p}).difference({*ref_3p}))

    return tools


def read_small_stats_per_gene(pf_data, parse_names=False):
    # type: (pd.DataFrame, bool) -> pd.DataFrame
    df = pd.read_csv(pf_data)
    if parse_names:
        df["Genome"] = df[["Genome"]].apply(fix_names, axis=1)
    return df
