# Author: Karl Gemayel
# Created: 7/10/20, 12:03 PM

import logging
import argparse
import math
import numpy as np
import pandas as pd
from typing import *
import matplotlib.pyplot as plt
import seaborn

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_general import Environment, add_env_args_to_parser
from mg_general.general import next_name, fix_names, get_value
from mg_stats.shelf import create_joint_reference_from_list, update_dataframe_with_stats, tidy_genome_level, \
    _helper_df_joint_reference
from mg_viz import sns

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_viz.general import FigureOptions
from mg_viz.shelf import get_order_by_rank

parser = argparse.ArgumentParser("Visualize statistics collected per gene by comparing to a reference set.")

parser.add_argument('--pf-data', required=True)
parser.add_argument('--ref-3p', required=True, nargs='+', help="List of tools to be used as 3prime reference "
                                                               "(ground-truth). Note: If more than one provided, "
                                                               "their intersection (in terms of 5' end is taken as"
                                                               " reference.")

parser.add_argument('--ref-5p', required=True, nargs='+', help="List of tools to be used as 5prime reference "
                                                               "(ground-truth). Note: If more than one provided, "
                                                               "their intersection (in terms of 5' end is taken as"
                                                               " reference.")

parser.add_argument('--tools', nargs="+", help="If set, only compare these tools. Otherwise all tools are chosen")
parser.add_argument('--parse-names', action='store_true', help="If set, try to shorten genome names. Useful only "
                                                               "genome ID's in the data are actually names")

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

            if df_group[f"3p:Match({tag_eq})"].sum() == 0:
                result[f"Match({tag})"] = np.nan
                result[f"Number of Error({tag})"] = np.nan
                result[f"Number of Found({tag})"] = np.nan
                result[f"Number of Predictions({tag})"] = np.nan
            else:
                result[f"Match({tag})"] = 100 * df_group[f"5p:Match({tag_eq})"].sum() / float(
                    df_group[f"3p:Match({tag_eq})"].sum())
                result[f"Number of Error({tag})"] = df_group[f"3p:Match({tag_eq})"].sum() - df_group[
                    f"5p:Match({tag_eq})"].sum()
                result[f"Number of Match({tag})"] = df_group[f"5p:Match({tag_eq})"].sum()

                result[f"Number of Found({tag})"] = df_group[f"3p:Match({tag_eq})"].sum()
                result[f"Number of Predictions({tag})"] = df[f"3p-{t}"].count()

                result[f"Sensitivity({t},{reference})"] = result[f"Number of Found({t},{reference})"] / df_group[
                    f"5p-{reference}"].count()
                result[f"Specificity({t},{reference})"] = result[f"Number of Found({t},{reference})"] / df_group[
                    f"5p-{t}"].count()

        result["Genome"] = gcfid
        result["Genome GC"] = df_group.at[df_group.index[0], "Genome GC"]
        result["Number in Reference"] = df_group[f"5p-{reference}"].count()

        list_entries.append(result)

    return pd.DataFrame(list_entries)


def viz_stats_as_function_of_reference_length(env, df_per_gene, tools, reference):
    # type: (Environment, pd.DataFrame, List[str], str) -> None

    list_entries = list()

    max_length = np.nanmax(df_per_gene[f"Length({reference})"])

    for genome, df_genome in df_per_gene.groupby("Genome", as_index=False):

        df_genome = df_genome.sort_values(f"Length({reference})")
        max_length = np.nanmax(df_genome[f"Length({reference})"])
        if math.isnan(max_length):
            continue

        for l in range(100, int(max_length) + 100, 100):
            curr_df = df_genome[df_genome[f"Length({reference})"] <= l]

            for t in tools:
                tag = ",".join([t, reference])
                tag_eq = "=".join([t, reference])

                list_entries.append({
                    "Genome": genome,
                    "Length": l,
                    "Tool": t,
                    "Number of Found": curr_df[f"3p:Match({tag_eq})"].sum(),
                    "Number of Matches": curr_df[f"5p:Match({tag_eq})"].sum(),
                    "Match Rate": 100 * curr_df[
                        f"5p:Match({tag_eq})"].sum() / float(curr_df[f"3p:Match({tag_eq})"].sum())
                })

    df = pd.DataFrame(list_entries)
    df["Genome"] = df.apply(fix_names, axis=1)

    if len(df["Genome"]) < 20:

        seaborn.lmplot("Length", "Number of Matches", data=df, hue="Genome")

        g = seaborn.FacetGrid(df, col="Genome", col_wrap=4, hue="Tool", sharey=False)
        # g.map(plt.plot, "Length", "Number of Found", linestyle="dashed")
        g.map(plt.plot, "Length", "Match Rate")
        # g.map(plt.plot, "x", "y_fit")
        g.set_xlabels("Length")
        g.set_titles("{col_name}")
        g.set(ylim=(0, 100))
        # g.set(xlim=(0, 5100))
        g.set_ylabels("Number of predictions")
        g.add_legend()
        plt.show()

        print("Hi")
    else:
        seaborn.lmplot("Length", "Number of Matches", data=df, hue="Tool")
        plt.show()


def viz_stats_per_gene_with_reference(env, df, tools, reference):
    # type: (Environment, pd.DataFrame, List[str], str) -> None

    # viz_stats_as_function_of_reference_length(env, df, tools, reference)
    df_gcfid = get_stats_at_gcfid_level_with_reference(df, tools, reference)
    # df_gcfid = df_gcfid.dropna().copy()

    if len(df_gcfid) < 10:

        # values on same plot
        df_tidy = pd.melt(df_gcfid, id_vars=["Genome"],
                          value_vars=[x for x in df_gcfid.columns if "Match(" in x],
                          var_name="Combination", value_name="Match")

        df_tidy["Combination"] = df_tidy["Combination"].apply(lambda x: x.split("(")[1].split(",")[0])
        combination_order = get_order_by_rank(df_tidy, "Genome", "Match", "Combination")

        df_tidy["Error"] = 100 - df_tidy["Match"]

        fig, ax = plt.subplots(figsize=(12, 4))

        sns.barplot(df_tidy, "Genome", "Error", hue="Combination", ax=ax,
                    sns_kwargs={"hue_order": combination_order}, show=False,
                    figure_options=FigureOptions(ylim=[0, 40]))

        fig.show()

        print(df_tidy.pivot("Genome", columns="Combination", values="Error").to_csv(float_format="%.2f"))

        # Number of errors
        df_tidy = pd.melt(df_gcfid, id_vars=["Genome"],
                          value_vars=[x for x in df_gcfid.columns if "Number of Error(" in x],
                          var_name="Combination", value_name="Error")
        df_tidy["Combination"] = df_tidy["Combination"].apply(lambda x: x.split("(")[1].split(",")[0])
        print(df_tidy.pivot("Genome", columns="Combination", values="Error").to_csv(float_format="%.2f"))

        fig, ax = plt.subplots(figsize=(12, 4))

        sns.barplot(df_tidy, "Genome", "Error", hue="Combination", ax=ax,
                    sns_kwargs={"hue_order": combination_order}, show=False,
                    figure_options=FigureOptions())

        fig.show()

        # errors as a function of gene length
    else:
        df_tidy = pd.melt(df_gcfid, id_vars=["Genome", "Genome GC"],
                          value_vars=[x for x in df_gcfid.columns if "Match(" in x],
                          var_name="Combination", value_name="Match")

        df_tidy["Combination"] = df_tidy["Combination"].apply(lambda x: x.split("(")[1].split(",")[0])
        combination_order = get_order_by_rank(df_tidy, "Genome", "Match", "Combination")

        df_tidy["Error"] = 100 - df_tidy["Match"]

        # fig, ax = plt.subplots(figsize=(8,4))
        g = seaborn.lmplot(
            "Genome GC", "Error", data=df_tidy, hue="Combination",
            lowess=True, scatter_kws={"alpha": 0.3, "s": 2},
            legend=False,
        )
        plt.legend(loc='upper left')
        for lh in plt.legend().legendHandles:
            # for lh in g._legend.legendHandles:
            lh.set_alpha(1)
            lh.set_sizes([4] * len(combination_order))

        plt.show()

        df_tidy_rate = df_tidy

        # Number of errors
        df_tidy = pd.melt(df_gcfid, id_vars=["Genome", "Genome GC"],
                          value_vars=[x for x in df_gcfid.columns if "Number of Error(" in x],
                          var_name="Combination", value_name="Error")
        df_tidy["Combination"] = df_tidy["Combination"].apply(lambda x: x.split("(")[1].split(",")[0])

        # fig, ax = plt.subplots(figsize=(8,4))
        g = seaborn.lmplot(
            "Genome GC", "Error", data=df_tidy, hue="Combination",
            lowess=True, scatter_kws={"alpha": 0.3, "s": 2},
            legend=False,
        )
        plt.legend(loc='upper left')
        for lh in plt.legend().legendHandles:
            # for lh in g._legend.legendHandles:
            lh.set_alpha(1)
            lh.set_sizes([4] * len(combination_order))

        plt.show()

        df_tidy_rate["Metric"] = "Percentage"
        df_tidy["Metric"] = "Number"

        df_total = pd.concat([df_tidy, df_tidy_rate], sort=False, ignore_index=True)
        # g = seaborn.FacetGrid(df_total, col="Metric", hue="Combination", sharey=False)
        # g.map(seaborn.lmplot, "Genome GC", "Error")
        # plt.show()
        #
        seaborn.lmplot("Genome GC", "Error", df_total, col="Metric", hue="Combination", sharey=False,
                       scatter_kws={"alpha": 0.3, "s": 2}, lowess=True,
                       legend=False,
                       )
        plt.legend(loc='upper left')
        for lh in plt.legend().legendHandles:
            # for lh in g._legend.legendHandles:
            lh.set_alpha(1)
            lh.set_sizes([4] * len(combination_order))

        plt.savefig(next_name(env["pd-work"]))
        plt.show()


def viz_stats_3p_number_of_predictions_number_of_found(env, df_tidy, reference):
    # type: (Environment, pd.DataFrame, str) -> None
    print(df_tidy.pivot(index="Genome", columns="Tool", values=["Number of Predictions", "Number of Found"]).to_csv())


def viz_stats_3p_sensitivity_specificity(env, df_tidy, reference):
    # type: (Environment, pd.DataFrame, str) -> None
    print(df_tidy.pivot(index="Genome", columns="Tool", values=["Sensitivity", "Specificity"]).to_csv())

def viz_stats_3p_missed_vs_length(env, df_per_gene, tools, reference):
    # type: (Environment, pd.DataFrame, List[str], str) -> None

    df_with_reference = df_per_gene[~df_per_gene[f"5p-{reference}"].isna()].copy()

    # sort by reference length
    df_with_reference.sort_values(f"Length({reference})", inplace=True)
    df_with_reference.reset_index(inplace=True)

    min_length = df_with_reference.iloc[0][f"Length({reference})"]


    list_entries = list()
    for t in tools:
        curr_length = min_length
        position = 0
        number_found = 0

        while True:

            while position < len(df_with_reference) and curr_length >= df_with_reference.at[
                df_with_reference.index[position], f"Length({reference})"
            ]:

                if not pd.isnull(df_with_reference.at[df_with_reference.index[position], f"5p-{t}"]):
                    number_found += 1
                position += 1

            list_entries.append({
                "Tool": t, "Found": number_found, "Reference": position, "Length": curr_length
            })

            if position >= len(df_with_reference):
                break
            curr_length = df_with_reference.at[df_with_reference.index[position], f"Length({reference})"]



    df = pd.DataFrame(list_entries)
    df["Rate"] = df["Found"] / df["Reference"]
    seaborn.lineplot("Length", "Rate", data=df[(~df["Tool"].isin({"gms2", "prodigal"})) & (df["Length"] < 500)],
                     hue="Tool")
    plt.show()
    print("hi")



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
    # df.set_index("Tool", drop=True, inplace=True)
    # df.reindex(tools).reset_index(inplace=True)

    def reorder(l_df, l_order):
        # type: (pd.DataFrame, List[str]) -> pd.DataFrame
        return l_df.reindex(l_order)

    print(df.pivot(index="Tool", columns="Bin Upper", values="Rate").to_csv())
    print(reorder(df.pivot(index="Tool", columns="Bin Upper", values="Missed"), [reference] + tools).to_csv())

    print(df.pivot(index="Tool", columns="Bin Upper", values="Predictions").to_csv())
    print(df.pivot(index="Tool", columns="Bin Upper", values="Exceed").reindex([reference] + tools).to_csv())












def viz_stats_3p(env, df_per_gene, tools, list_ref):
    # type: (Environment, pd.DataFrame, List[str], List[str]) -> None
    """Visualize statistics at 3prime level"""

    reference = _helper_df_joint_reference(df_per_gene, list_ref)
    df_per_gene = update_dataframe_with_stats(df_per_gene, tools, reference).copy()

    #### Genome Level
    # compute stats per genome
    df_stats_gcfid = list()
    for _, df_group in df_per_gene.groupby("Genome", as_index=False):
        df_stats_gcfid.append(get_stats_at_gcfid_level_with_reference(df_group, tools, reference))
    df_per_genome = pd.concat(df_stats_gcfid, ignore_index=True, sort=False)

    df_tidy = tidy_genome_level(env, df_per_genome)
    df_tidy = df_tidy[df_tidy["Tool"].apply(lambda x: x.lower()).isin(tools + [reference])]
    # Number of Predictions, number of found
    viz_stats_3p_number_of_predictions_number_of_found(env, df_tidy, reference)

    # Number of Predictions, Precision
    # viz_stats_3p_number_of_predictions_precision(env, df_tidy, reference)
    viz_stats_3p_sensitivity_specificity(env, df_tidy, reference)

    #### Gene Level

    viz_stats_3p_missed_vs_length(env, df_per_gene, tools, reference)


def viz_stats_5p_number_of_found_number_of_5prime_match(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    print(df_tidy.pivot(index="Genome", columns="Tool", values=["Number of Found", "Number of Match"]).to_csv())


def viz_stats_5p_number_of_found_number_of_5prime_error(env, df_tidy, reference, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, Any]) -> None
    tool_order = get_value(kwargs, "tool_order", sorted(df_tidy["Tool"].unique()))
    print(df_tidy.pivot(
        index="Genome", columns="Tool", values=["Number of Found", "Number of Error"]
    ).reorder_levels([1,0], 1)[[x.upper() for x in tool_order]].reorder_levels([1,0],1).sort_index(1,0, sort_remaining=False).to_csv(
    ))


def viz_stats_5p(env, df_per_gene, tools, list_ref):
    # type: (Environment, pd.DataFrame, List[str], List[str]) -> None
    """Visualize statistics at 5prime level"""
    reference = _helper_df_joint_reference(df_per_gene, list_ref)
    df_per_gene = update_dataframe_with_stats(df_per_gene, tools, reference).copy()

    #### Genome Level
    # compute stats per genome
    df_stats_gcfid = list()
    for _, df_group in df_per_gene.groupby("Genome", as_index=False):
        df_stats_gcfid.append(get_stats_at_gcfid_level_with_reference(df_group, tools, reference))
    df_per_genome = pd.concat(df_stats_gcfid, ignore_index=True, sort=False)

    df_tidy = tidy_genome_level(env, df_per_genome)
    df_tidy = df_tidy[df_tidy["Tool"].apply(lambda x: x.lower()).isin(tools + [reference])]
    # Number of Found, number of 5prime match
    viz_stats_5p_number_of_found_number_of_5prime_match(env, df_tidy, reference)
    viz_stats_5p_number_of_found_number_of_5prime_error(env, df_tidy, reference, tool_order=tools)
    # Number of Predictions, Precision
    # viz_stats_3p_number_of_predictions_precision(env, df_tidy, reference)
    # viz_stats_5p_sensitivity_specificity(env, df_tidy, reference)

    #### Gene Level


def viz_stats_per_gene(env, df_per_gene, tools, list_ref_5p, list_ref_3p):
    # type: (Environment, pd.DataFrame, List[str], List[str], List[str]) -> None

    viz_stats_3p(env, df_per_gene, tools, list_ref_3p)
    viz_stats_5p(env, df_per_gene, tools, list_ref_5p)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = pd.read_csv(args.pf_data)
    if args.parse_names:
        df["Genome"] = df[["Genome"]].apply(fix_names, axis=1)

    # get tools list
    # If not provided, extract from df
    # Make sure it doesn't contain any references
    all_tools = sorted(set([x.split("-")[1] for x in df.columns if "5p-" in x]))

    # check that references exist
    for list_ref in [args.ref_5p, args.ref_3p]:
        for ref in list_ref:
            if ref not in all_tools:
                raise ValueError(f"Unknown reference {ref}")

    if args.tools is not None:
        tools = args.tools
    else:
        tools = all_tools
        tools = sorted(set(tools).difference({*args.ref_5p}).difference({*args.ref_3p}))

    # update_dataframe_with_stats(df, tools, args.ref_5p, args.ref_3p)

    viz_stats_per_gene(env, df, tools, args.ref_5p, args.ref_3p)


if __name__ == "__main__":
    main(my_env, parsed_args)
