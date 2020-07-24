# Author: Karl Gemayel
# Created: 6/29/20, 3:41 PM

import logging
import argparse
import numpy as np
import pandas as pd
from typing import *
import matplotlib.pyplot as plt
import seaborn

# noinspection All
from seaborn import FacetGrid

import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_general import Environment, add_env_args_to_parser
from mg_general.general import all_elements_equal, fix_names
from mg_general.shelf import powerset
from mg_stats.shelf import create_joint_reference_from_list
from mg_viz import sns
from mg_viz.general import FigureOptions, square_subplots

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Visualize statistics collected per gene.")

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

        result = dict()
        for t in tools:

            tag = ",".join([t, reference])
            tag_eq = "=".join([t, reference])

            if df_group[f"3p:Match({tag_eq})"].sum() == 0:
                result[f"Match({tag})"] = np.nan
                result[f"Number of Error({tag})"] = np.nan
                result[f"Error Rate({tag})"] = np.nan
            else:
                result[f"Match({tag})"] = 100 * df_group[f"5p:Match({tag_eq})"].sum() / float(
                    df_group[f"3p:Match({tag_eq})"].sum())
                result[f"Error Rate({tag})"] = 100 - result[f"Match({tag})"]
                result[f"Number of Error({tag})"] = df_group[f"3p:Match({tag_eq})"].sum() - df_group[
                    f"5p:Match({tag_eq})"].sum()

        result["Genome"] = gcfid
        result["Genome GC"] = df_group.at[df_group.index[0], "Genome GC"]
        result["Chunk Size"] = df_group.at[df_group.index[0], "Chunk Size"]

        list_entries.append(result)

    return pd.DataFrame(list_entries)


# def get_stats_at_gcfid_level(df, tools, reference):
#     # type: (pd.DataFrame, List[str], str) -> pd.DataFrame
#
#     list_entries = list()
#
#     ps = powerset(tools, min_len=2)
#
#
#
#     for gcfid, df_group in df.groupby("Genome", as_index=False):
#         result = dict()
#
#         for comb in ps:
#             tag = ",".join(comb)
#             tag_eq = "=".join(comb)
#
#             result[f"Match({tag})"] = 100 * df_group[f"5p:{tag_eq}"].sum()/ float(df_group[f"3p:{tag_eq}"].sum())
#
#         result["Genome"] = gcfid
#         result["Chunk Size"] = df_group.at[df_group.index[0], "Chunk Size"]
#         list_entries.append(result)
#
# return pd.DataFrame(list_entries)


def viz_stats_at_gcfid_level(df, tools):
    pass


def ridgeplot(df):
    # Create the data
    names = sorted(set(df["Genome"]))
    x = df["GC Diff"].values

    g = df.apply(lambda r: f"{r['Genome']} ({r['Genome GC']:.2f})", axis=1)
    df = pd.DataFrame(dict(x=x, g=g))

    hue_order = sorted(set(g), key=lambda x: float(x.split("(")[1].split(")")[0]))
    # Initialize the FacetGrid object
    pal = seaborn.cubehelix_palette(10, rot=-.25, light=.7)
    g = seaborn.FacetGrid(df, row="g", hue="g",
                          hue_order=hue_order,
                          row_order=hue_order,
                          aspect=15, height=.5, palette=pal)

    # Draw the densities in a few steps
    g.map(seaborn.kdeplot, "x", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
    g.map(seaborn.kdeplot, "x", clip_on=False, color="w", lw=2, bw=.2)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)

    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight="bold", color=color,
                ha="left", va="center", transform=ax.transAxes)

    g.map(label, "x")

    # Set the subplots to overlap
    g.fig.subplots_adjust(hspace=-.25)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)
    plt.show()


def viz_stats_per_gene(env, df, tools, reference):
    # type: (Environment, pd.DataFrame, List[str], str) -> None

    df_stats_gcfid = list()
    for _, df_group in df.groupby("Chunk Size", as_index=False):
        df_stats_gcfid.append(get_stats_at_gcfid_level_with_reference(df_group, tools, reference))
    df_all_chunks = pd.concat(df_stats_gcfid, ignore_index=True, sort=False)

    for chunk_size, df_stats_gcfid in df_all_chunks.groupby("Chunk Size", as_index=False):
        # ps = powerset(tools, min_len=2)

        # for t in tools:
        #     tag = ",".join([t, reference])
        #     tag_eq = "=".join([t, reference])
        #     sns.barplot(df_stats_gcfid, "Genome", f"Match({tag})", figure_options=FigureOptions(ylim=[50, 100]))

        # values on same plot
        df_tidy = pd.melt(df_stats_gcfid, id_vars=["Genome"],
                          value_vars=[x for x in df_stats_gcfid.columns if "Match(" in x],
                          var_name="Combination", value_name="Match")

        fig, ax = plt.subplots(figsize=(12, 4))
        sns.barplot(df_tidy, "Genome", "Match", hue="Combination", ax=ax, show=False, figure_options=FigureOptions(
            ylim=[60, 100]
        ))
        ax.set_title(f"Chunk Size: {chunk_size}")
        plt.show()

        # values on same plot
        df_tidy = pd.melt(df_stats_gcfid, id_vars=["Genome"],
                          value_vars=[x for x in df_stats_gcfid.columns if "Error Rate(" in x],
                          var_name="Combination", value_name="Error Rate")

        fig, ax = plt.subplots(figsize=(12, 4))
        sns.barplot(df_tidy, "Genome", "Error Rate", hue="Combination", ax=ax, show=False,
                    figure_options=FigureOptions(
                        # ylim=[60, 100]
                    ))
        ax.set_title(f"Chunk Size: {chunk_size}")
        plt.show()

        # values on same plot
        df_tidy = pd.melt(df_stats_gcfid, id_vars=["Genome", "Chunk Size"],
                          value_vars=[x for x in df_stats_gcfid.columns if "Number of Error(" in x],
                          var_name="Combination", value_name="Number of Error")

        fig, ax = plt.subplots(figsize=(12, 4))
        sns.barplot(df_tidy, "Genome", "Number of Error", hue="Combination", ax=ax, show=False,
                    figure_options=FigureOptions(
                        # ylim=[60, 100]
                    ))
        ax.set_title(f"Chunk Size: {chunk_size}")
        plt.show()

        # used
        g = seaborn.FacetGrid(df_tidy, col="Genome", col_wrap=4)
        g.map(seaborn.barplot, "Combination", "Number of Error", order=sorted(df_tidy["Combination"].unique()))
        # g.map(plt.plot, "x", "y_fit")
        g.set_xlabels("GC")
        g.set_ylabels("Probability")
        g.add_legend()
        plt.show()

        # # GC
        # sns.lmplot(df, "Genome GC", "Gene GC", hue="Genome")
        #
        # df["GC Diff"] = df["Genome GC"] - df["Gene GC"]
        #
        # sorted_genomes = sorted(set(df["Genome"]))
        # num_genomes = len(sorted_genomes)
        # num_rows, num_cols = square_subplots(num_genomes)
        #
        # fig, axes = plt.subplots(num_rows, num_cols, sharex="all", sharey="all")
        #
        # for i in range(num_genomes):
        #     genome = sorted_genomes[i]
        #     ax = axes.ravel()[i]
        #
        #     sns.distplot(df[df["Genome"] == genome], "GC Diff", ax=ax, show=False)
        #
        # plt.show()

        # sns.lmplot(df, "Genome GC", "GC Diff")
        #
        # ridgeplot(df)


def tools_match_for_dataframe_row(r, tools):
    # type: (pd.Series, Iterable[str]) -> bool

    # check all tools make a prediction for current gene
    list_5ps = list()

    for t in tools:
        if r[f"5p-{t}"] is None:
            return False
        list_5ps.append(r[f"5p-{t}"])

    return all_elements_equal(list_5ps)


# def update_dataframe_with_stats(df):
#     # type: (pd.DataFrame) -> None
#
#     tools = sorted([x.split("-")[1] for x in df.columns if "5p-" in x])
#
#     ps = powerset(tools, min_len=2)
#
#     # get all combinations of tools, and set default to False (i.e. not matching by 5' end)
#     for comb in ps:
#         tag_5p = "5p:" + "=".join(comb)
#         tag_3p = "3p:" + "=".join(comb)
#
#         # match by 5prime end
#         df[tag_5p] = df.apply(tools_match_for_dataframe_row, axis=1, tools=comb)
#
#         # all tools have a prediction
#         df[tag_3p] = df[[f"5p-{t}" for t in comb]].notnull().all(axis=1)

def update_dataframe_with_stats(df, tools, reference):
    # type: (pd.DataFrame, List[str], str) -> None

    for t in tools:
        tag_5p = f"5p:Match({t}={reference})"
        tag_3p = f"3p:Match({t}={reference})"

        # match by 5prime end
        df[tag_5p] = df[f"5p-{t}"] == df[f"5p-{reference}"]

        # all tools have a prediction
        df[tag_3p] = df[[f"5p-{t}", f"5p-{reference}"]].notnull().all(axis=1)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = pd.read_csv(args.pf_data)
    df["Genome"] = df[["Genome"]].apply(fix_names, axis=1)

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
    viz_stats_per_gene(env, df, tools, reference)


if __name__ == "__main__":
    main(my_env, parsed_args)
