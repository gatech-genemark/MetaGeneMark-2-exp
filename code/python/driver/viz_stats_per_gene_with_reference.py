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
from mg_general.general import next_name, fix_names
from mg_stats.shelf import create_joint_reference_from_list
from mg_viz import sns

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_viz.general import FigureOptions
from mg_viz.shelf import get_order_by_rank

parser = argparse.ArgumentParser("Visualize statistics collected per gene by comparing to a reference set.")

parser.add_argument('--pf-data', required=True)
parser.add_argument('--reference', required=True, nargs='+', help="List of tools to be used as reference "
                                                                  "(ground-truth). Note: If more than one provided, "
                                                                  "their intersection (in terms of 5' end is taken as"
                                                                  " reference.")
parser.add_argument('--tools', nargs="+", help="If set, only compare these tools. Otherwise all tools are chosen")

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

        if gcfid == 'GCF_000661085.1_Myco_sp_TKK-01-0051_V1':
            print("yo")

        result = dict()
        for t in tools:

            tag = ",".join([t, reference])
            tag_eq = "=".join([t, reference])
            
            if df_group[f"3p:Match({tag_eq})"].sum() == 0:
                result[f"Match({tag})"] = np.nan
                result[f"Number of Error({tag})"] = np.nan
            else:
                result[f"Match({tag})"] = 100 * df_group[f"5p:Match({tag_eq})"].sum() / float(df_group[f"3p:Match({tag_eq})"].sum())
                result[f"Number of Error({tag})"] = df_group[f"3p:Match({tag_eq})"].sum() - df_group[f"5p:Match({tag_eq})"].sum()


        result["Genome"] = gcfid
        result["Genome GC"] = df_group.at[df_group.index[0], "Genome GC"]
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

        for l in range(100, int(max_length)+100, 100):
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
    df_gcfid = df_gcfid.dropna().copy()

    if len(df_gcfid) < 10:

        # values on same plot
        df_tidy = pd.melt(df_gcfid, id_vars=["Genome"],
                          value_vars=[x for x in df_gcfid.columns if "Match(" in x],
                          var_name="Combination", value_name="Match")


        df_tidy["Combination"] = df_tidy["Combination"].apply(lambda x: x.split("(")[1].split(",")[0])
        combination_order = get_order_by_rank(df_tidy, "Genome", "Match", "Combination")

        df_tidy["Error"] = 100 - df_tidy["Match"]

        fig, ax = plt.subplots(figsize=(12,4))

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


def update_dataframe_with_stats(df, tools, reference):
    # type: (pd.DataFrame, List[str], str) -> None

    for t in tools:
        tag_5p = f"5p:Match({t}={reference})"
        tag_3p = f"3p:Match({t}={reference})"

        # match by 5prime end
        df[tag_5p] = df[f"5p-{t}"] == df[f"5p-{reference}"]

        # all tools have a prediction
        df[tag_3p] = df[[f"5p-{t}", f"5p-{reference}"]].notnull().all(axis=1)

    df[f"Length({reference})"] = df.apply(
        lambda r: abs(r[f"3p-{reference}"] - r[f"5p-{reference}"]) + 1,
        axis=1
    )



def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = pd.read_csv(args.pf_data)
    # df["Genome"] = df[["Genome"]].apply(fix_names, axis=1)
    # df = df.sample(500)

    reference = args.reference
    tools = sorted(
        set([x.split("-")[1] for x in df.columns if "5p-" in x])
    )
    if args.tools is not None:
        tools = args.tools + args.reference

    # check that reference is one of the tools
    if len(reference) > 1:
        # check that all are part of tools
        for r in reference:
            if r not in tools:
                raise ValueError(f"Unknown reference {r}")

        tools = sorted(set(tools).difference(
            {*reference}
        ))
        reference = create_joint_reference_from_list(df, reference)


    else:
        reference = reference[0]
        if reference not in tools:
            raise ValueError(f"Unknown reference {reference}")

        tools = sorted(set(tools).difference(
            {reference}
        ))

    update_dataframe_with_stats(df, tools, reference)

    # tools = sorted([x.split("-")[1] for x in df.columns if "5p-" in x])
    viz_stats_per_gene_with_reference(env, df, tools, reference)


if __name__ == "__main__":
    main(my_env, parsed_args)
