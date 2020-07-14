# Author: Karl Gemayel
# Created: 7/10/20, 12:03 PM

import logging
import argparse
import pandas as pd
from typing import *
import matplotlib.pyplot as plt

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_general import Environment, add_env_args_to_parser
from mg_general.general import fix_names
from mg_viz import sns

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_viz.general import FigureOptions

parser = argparse.ArgumentParser("Visualize statistics collected per gene by comparing to a reference set.")

parser.add_argument('--pf-data', required=True)
parser.add_argument('--reference', required=True, nargs='+', help="List of tools to be used as reference "
                                                                  "(ground-truth). Note: If more than one provided, "
                                                                  "their intersection (in terms of 5' end is taken as"
                                                                  " reference.")

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


def get_stats_at_gcfid_level_with_reference(df, tools, reference):
    # type: (pd.DataFrame, List[str], str) -> pd.DataFrame

    list_entries = list()

    for gcfid, df_group in df.groupby("Genome", as_index=False):


        result = dict()
        for t in tools:

            tag = ",".join([t, reference])
            tag_eq = "=".join([t, reference])

            result[f"Match({tag})"] = 100 * df_group[f"5p:Match({tag_eq})"].sum()/ float(df_group[f"3p:Match({tag_eq})"].sum())

        result["Genome"] = gcfid
        list_entries.append(result)

    return pd.DataFrame(list_entries)


def viz_stats_per_gene_with_reference(env, df, tools, reference):
    # type: (Environment, pd.DataFrame, List[str], str) -> None

    df_gcfid = get_stats_at_gcfid_level_with_reference(df, tools, reference)

    # values on same plot
    df_tidy = pd.melt(df_gcfid, id_vars=["Genome"],
                      value_vars=[x for x in df_gcfid.columns if "Match(" in x],
                      var_name="Combination", value_name="Match")

    fig, ax = plt.subplots(figsize=(12,4))
    sns.barplot(df_tidy, "Genome", "Match", hue="Combination", ax=ax,
                figure_options=FigureOptions(ylim=[60, 100]))

    fig.show()

def update_dataframe_with_stats(df, tools, reference):
    # type: (pd.DataFrame, List[str], str) -> None



    for t in tools:
        tag_5p = f"5p:Match({t}={reference})"
        tag_3p = f"3p:Match({t}={reference})"

        # match by 5prime end
        df[tag_5p] = df[f"5p-{t}"] == df[f"5p-{reference}"]

        # all tools have a prediction
        df[tag_3p] = df[[f"5p-{t}", f"5p-{reference}"]].notnull().all(axis=1)


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



    return reference


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = pd.read_csv(args.pf_data)
    df["Genome"] = df[["Genome"]].apply(fix_names, axis=1)

    reference = args.reference
    tools = sorted(
        set([x.split("-")[1] for x in df.columns if "5p-" in x])
    )

    # check that reference is one of the tools
    if isinstance(reference, list):
        # check that all are part of tools
        for r in reference:
            if r not in tools:
                raise ValueError(f"Unknown reference {r}")

        tools = sorted(set(tools).difference(
            {*reference}
        ))
        reference = create_joint_reference_from_list(df, reference)


    else:
        if reference not in tools:
            raise ValueError(f"Unknown reference {reference}")

        tools = sorted(set(tools).difference(
            {reference, "prodigal", "gms2"}
        ))

    update_dataframe_with_stats(df, tools, reference)

    # tools = sorted([x.split("-")[1] for x in df.columns if "5p-" in x])
    viz_stats_per_gene_with_reference(env, df, tools, reference)


if __name__ == "__main__":
    main(my_env, parsed_args)
