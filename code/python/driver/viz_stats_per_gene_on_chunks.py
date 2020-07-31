# Author: Karl Gemayel
# Created: 6/29/20, 3:41 PM

import logging
import seaborn
import argparse
import numpy as np
import pandas as pd
from typing import *
from functools import reduce
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_viz.general import square_subplots
from mg_general import Environment, add_env_args_to_parser
from mg_stats.shelf import update_dataframe_with_stats, tidy_genome_level, \
    _helper_df_joint_reference
from mg_general.general import all_elements_equal, fix_names, next_name
from mg_viz.colormap import ColorMap as CM

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_viz.shelf import number_formatter

parser = argparse.ArgumentParser("Visualize statistics collected per gene.")

parser.add_argument('--pf-data', required=True)

parser.add_argument('--ref-5p', required=False, nargs="+", help="Reference(s) on which to compare 5' predictions")
parser.add_argument('--ref-3p', required=False, nargs="+", help="Reference(s) on which to compare 3' predictions")

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
                result[f"Error Rate({tag})"] = np.nan
                result[f"Number of Error({tag})"] = np.nan
                result[f"Number of Match({tag})"] = np.nan
                # result[f"Number of Predictions({t},{t})"] = np.nan
            else:
                result[f"Match({tag})"] = 100 * df_group[f"5p:Match({tag_eq})"].sum() / float(
                    df_group[f"3p:Match({tag_eq})"].sum())
                result[f"Error Rate({tag})"] = 100 - result[f"Match({tag})"]
                result[f"Number of Error({tag})"] = df_group[f"3p:Match({tag_eq})"].sum() - df_group[
                    f"5p:Match({tag_eq})"].sum()

                result[f"Number of Found({tag})"] = df_group[f"3p:Match({tag_eq})"].sum()
                result[f"Number of Missed({tag})"] = df_group[f"5p-{reference}"].count() - df_group[
                    f"3p:Match({tag_eq})"].sum()
                result[f"Number of Match({tag})"] = df_group[f"5p:Match({tag_eq})"].sum()
                # result[f"Number of Predictions({t},{t})"] = df_group[f"5p-{t}"].count()

                result[f"Number of IC5p Match({tag})"] = (
                        df_group[f"5p:Match({tag_eq})"] & df_group[f"Partial5p-{reference}"]).sum()
                result[f"Number of IC5p Found({tag})"] = (
                        df_group[f"3p:Match({tag_eq})"] & df_group[f"Partial5p-{reference}"]).sum()
                result[f"IC5p Match({tag})"] = 100 * result[f"Number of IC5p Match({tag})"] / result[
                    f"Number of IC5p Found({tag})"]

                result[f"Number of IC3p Match({tag})"] = (
                        df_group[f"5p:Match({tag_eq})"] & df_group[f"Partial3p-{reference}"]).sum()
                result[f"Number of IC3p Found({tag})"] = (
                        df_group[f"3p:Match({tag_eq})"] & df_group[f"Partial3p-{reference}"]).sum()
                result[f"IC3p Match({tag})"] = 100 * result[f"Number of IC3p Match({tag})"] / result[
                    f"Number of IC3p Found({tag})"]

                result[f"Number of Comp Match({tag})"] = (
                        df_group[f"5p:Match({tag_eq})"] & ~(
                        df_group[f"Partial5p-{reference}"] | df_group[f"Partial3p-{reference}"])).sum()
                result[f"Number of Comp Found({tag})"] = (
                        df_group[f"3p:Match({tag_eq})"] & ~(
                        df_group[f"Partial5p-{reference}"] | df_group[f"Partial3p-{reference}"])).sum()
                result[f"Comp Match({tag})"] = 100 * result[f"Number of Comp Match({tag})"] / result[
                    f"Number of Comp Found({tag})"]

        for t in tools + [reference]:
            result[f"Number of Predictions({t},{t})"] = df_group[f"5p-{t}"].count()
            result[f"Runtime({t},{t})"] = df_group[f"Runtime"].mean()
            if t != reference:
                result[f"Precision({t},{reference})"] = result[f"Number of Found({t},{reference})"] / result[
                    f"Number of Predictions({t},{t})"]
                result[f"Recall({t},{reference})"] = result[f"Number of Found({t},{reference})"] / df_group[
                    f"5p-{reference}"].count()
                result[f"WR({t},{reference})"] = (result[f"Number of Predictions({t},{t})"] - result[
                    f"Number of Found({t},{reference})"]) / result[f"Number of Predictions({t},{t})"]

                result[f"Sensitivity({t},{reference})"] = result[f"Number of Found({t},{reference})"] / df_group[
                    f"5p-{reference}"].count()
                result[f"Specificity({t},{reference})"] = result[f"Number of Found({t},{reference})"] / result[
                    f"Number of Predictions({t},{t})"]

            # result[f"Runtime({t, t})"] = df_group[f"Runtime"].mean()

        result["Genome"] = gcfid
        result["Genome GC"] = df_group.at[df_group.index[0], "Genome GC"]
        result["Chunk Size"] = df_group.at[df_group.index[0], "Chunk Size"]
        result["Number in Reference"] = result[f"Number of Predictions({reference},{reference})"]

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


def number_and_match(env, df_total, hue_order, col_number, col_perc, sup_title):
    g = seaborn.FacetGrid(df_total, col="Genome", col_wrap=4, hue="Tool", sharey=False)
    g.map(plt.plot, "Chunk Size", col_number, linestyle="dashed")
    # g.map(plt.plot, "x", "y_fit")
    g.set_xlabels("Chunk Size")
    g.set_titles("{col_name}")
    g.set(ylim=(0, None))
    g.set(xlim=(0, 5100))
    g.set_ylabels("Number of predictions")
    g.add_legend()

    for ax, (_, subdata) in zip(g.axes, df_total.groupby('Genome')):
        ax2 = ax.twinx()
        subdata = subdata.sort_values("Chunk Size")
        for hue in hue_order:
            subdata_hue = subdata[subdata["Tool"] == hue]
            ax2.plot(subdata_hue["Chunk Size"], subdata_hue[col_perc], label=hue)
            ax2.set_ylim(40, 100)
        # subdata.plot(x='data_sondage', y='impossible', ax=ax2, legend=False, color='r')

    plt.tight_layout()
    plt.savefig(next_name(env["pd-work"]))
    plt.suptitle(sup_title)
    plt.show()


def viz_number_of_predictions_for_short(env, df):
    # type: (Environment, pd.DataFrame) -> None

    df = df[df["Tool"] != "VERIFIED"]
    hue_order = sorted(df["Tool"].unique())
    df["Found%"] = 100 * df["Number of Found"] / df["Number of Predictions"]

    g = seaborn.FacetGrid(df, col="Genome", col_wrap=4, hue="Tool", sharey=True)

    xlim = (0, 5100)

    g.map(plt.plot, "Chunk Size", "Number of Found", linestyle="dashed")
    g.map(plt.plot, "Chunk Size", "Number of Predictions")

    for ax in g.axes:
        ax.yaxis.set_major_formatter(FuncFormatter(number_formatter))

    # g.map(plt.plot, "x", "y_fit")
    g.set_xlabels("Chunk Size (nt)")
    g.set_titles("{col_name}", style="italic")
    g.set(ylim=(0, None))
    g.set(xlim=xlim)
    g.set_ylabels("Number of predictions")

    # for ax, (_, subdata) in zip(g.axes, df.groupby('Genome')):
    #     # ax2 = ax.twinx()
    #     ax2 = inset_axes(ax, width="50%", height="50%", loc=1, borderpad=1)
    #
    #     subdata = subdata.sort_values("Chunk Size")
    #     for hue in hue_order:
    #         subdata_hue = subdata[subdata["Tool"] == hue]
    #         ax2.plot(subdata_hue["Chunk Size"], subdata_hue["Found%"], label=hue)
    #         ax2.set_ylim(40,100)
    #         ax2.set_ylabel("TPR")
    #         ax2.set_xlim(*xlim)
    #         ax2.set_xticks([])
    #         ax2.set_yticks([])
    #
    #     # subdata.plot(x='data_sondage', y='impossible', ax=ax2, legend=False, color='r')

    # plt.tight_layout()
    g.add_legend()

    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    g = seaborn.FacetGrid(df, col="Genome", col_wrap=4, hue="Tool", sharey=False)

    g.map(plt.plot, "Chunk Size", "Found%")
    # g.map(plt.plot, "Chunk Size", "Number of Found", linestyle="dashed")
    # g.map(plt.plot, "x", "y_fit")
    g.set_xlabels("Chunk Size (nt)")
    g.set_titles("{col_name}", style="italic")
    g.set(ylim=(0, None))
    g.set(xlim=(0, 5100))
    g.set_ylabels("Number of predictions")
    g.add_legend()

    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    df = df[df["Tool"].isin({"MGM2", "MPRODIGAL", "FGS", "MGA"})]
    g = seaborn.FacetGrid(df, col="Genome", col_wrap=4, hue="Tool", sharey=True)

    g.map(plt.plot, "Recall", "Precision")
    # g.map(plt.plot, "Chunk Size", "Number of Found", linestyle="dashed")
    # g.map(plt.plot, "x", "y_fit")
    # g.set_xlabels("Chunk Size (nt)")
    g.set_titles("{col_name}", style="italic")
    g.set(ylim=(0, 1))
    g.set(xlim=(0, 1))
    # g.set_ylabels("Number of predictions")
    g.add_legend()

    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    g = seaborn.FacetGrid(df, col="Genome", col_wrap=4, hue="Tool", sharey=True, palette=CM.get_map("tools"))

    g.map(plt.plot, "Chunk Size", "WR", linestyle="dashed")
    g.map(plt.plot, "Chunk Size", "Precision")

    # g.map(plt.plot, "Chunk Size", "Number of Found", linestyle="dashed")
    # g.map(plt.plot, "x", "y_fit")
    g.set_titles("{col_name}", style="italic")
    g.set(ylim=(0, 1))
    g.set(xlim=(0, 5100))
    g.set_xlabels("Chunk Size (nt)")
    g.set_ylabels("Score")
    g.add_legend()

    # inset_ylim=(0, max(df["Number of Predictions"])+100)
    # for ax, (_, subdata) in zip(g.axes, df.groupby('Genome')):
    #     # ax2 = ax.twinx()
    #     ax2 = inset_axes(ax, width="40%", height="40%", loc=7, borderpad=1)
    #
    #     subdata = subdata.sort_values("Chunk Size")
    #     for hue in hue_order:
    #         subdata_hue = subdata[subdata["Tool"] == hue]
    #         ax2.plot(subdata_hue["Chunk Size"], subdata_hue["Number of Predictions"], label=hue,
    #                  color=CM.get_map("tools")[hue])
    #         # ax2.set_ylim(40,100)
    #         ax2.set_ylabel("Total Predictions")
    #         ax2.set_xlim(*xlim)
    #         ax2.set_ylim(*inset_ylim)
    #
    #         ax2.yaxis.set_major_formatter(FuncFormatter(number_formatter))
    #         ax2.set_xticks([])
    #         # ax2.set_yticks([])
    #
    #     # subdata.plot(x='data_sondage', y='impossible', ax=ax2, legend=False, color='r')

    plt.savefig(next_name(env["pd-work"]))
    plt.show()


def viz_plot_per_genome_y_error_x_chunk(env, df):
    genomes = sorted(df["Genome"].unique())
    nrows, ncols = square_subplots(len(genomes))

    values_to_melt = ["Match", "Number of Error", "Number of Found", "Number of Match", "Number of Predictions",
                      "Number of IC5p Match", "Number of IC5p Found", "Number of IC3p Match", "Number of IC3p Found",
                      "Number of Comp Match", "Number of Comp Found", "Precision", "Recall", "WR", "Number of Missed",
                      "IC3p Match", "IC5p Match", "Comp Match"]
    df_total = list()
    for v in values_to_melt:
        if v == "Precision":
            print('hi')
        df_curr = pd.melt(df, id_vars=["Genome", "Chunk Size", "Genome GC"],
                          value_vars=[x for x in df.columns if v == x.split("(")[0].strip()],
                          var_name="Combination", value_name=v)
        df_curr["Tool"] = df_curr["Combination"].apply(lambda x: x.split("(")[1].split(",")[0].upper())
        df_total.append(df_curr)

    df_total = reduce(lambda df1, df2: pd.merge(df1, df2, on=["Genome", "Chunk Size", "Genome GC", "Tool"],
                                                how="outer"), df_total)
    viz_number_of_predictions_for_short(env, df_total)
    # return
    # df_total = pd.melt(
    #     df_total,
    #     id_vars=["Genome", "Chunk Size", "Genome GC", "Combination"],
    #     value_vars=values_to_melt,
    #     var_name="Metric", value_name="Score")
    # dfs = [df_tmp.set_index(["Genome", "Chunk Size", "Genome GC"]) for df_tmp in df_total]
    # dfs = pd.concat(dfs, ignore_index=True, sort=False, axis=1)
    hue_order = sorted(df_total["Tool"].unique())
    g = seaborn.FacetGrid(df_total, col="Genome", col_wrap=4, hue="Tool", sharey=True, hue_order=hue_order)
    # g.map(plt.plot, "Chunk Size", "Match", marker="o")
    g.map(plt.plot, "Chunk Size", "Number of Found", linestyle="--")
    # g.map(plt.plot, "x", "y_fit")
    g.set_xlabels("Chunk Size")
    g.set_ylabels("Metric")
    g.set(ylim=(0, None))
    g.set(xlim=(None, None))
    g.add_legend()

    for ax, (_, subdata) in zip(g.axes, df_total.groupby('Genome')):
        ax2 = ax.twinx()
        subdata = subdata.sort_values("Chunk Size")
        for hue in hue_order:
            subdata_hue = subdata[subdata["Tool"] == hue]
            ax2.plot(subdata_hue["Chunk Size"], subdata_hue["Match"], label=hue)
            ax2.set_ylim(40, 100)
        # subdata.plot(x='data_sondage', y='impossible', ax=ax2, legend=False, color='r')

    plt.tight_layout()
    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    g = seaborn.FacetGrid(df_total, col="Genome", col_wrap=4, hue="Tool", sharey=False)

    g.map(plt.plot, "Chunk Size", "Number of Predictions")
    # g.map(plt.plot, "x", "y_fit")
    g.set_xlabels("Chunk Size")
    g.set_titles("{col_name}")
    g.set(ylim=(0, None))
    # g.set(xlim=(None,5100))
    g.set_ylabels("Number of predictions")
    g.add_legend()

    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    # Incomplete
    g = seaborn.FacetGrid(df_total, col="Genome", col_wrap=4, hue="Tool", sharey=False)
    g.map(plt.plot, "Chunk Size", "Number of IC5p Found", linestyle="dashed")
    # g.map(plt.plot, "x", "y_fit")
    g.set_xlabels("Chunk Size")
    g.set_titles("{col_name}")
    g.set(ylim=(0, None))
    # g.set(xlim=(0, 5100))
    g.set_ylabels("Number of predictions")
    g.add_legend()

    for ax, (_, subdata) in zip(g.axes, df_total.groupby('Genome')):
        ax2 = ax.twinx()
        subdata = subdata.sort_values("Chunk Size")
        for hue in hue_order:
            subdata_hue = subdata[subdata["Tool"] == hue]
            ax2.plot(subdata_hue["Chunk Size"], subdata_hue["IC5p Match"], label=hue)
            ax2.set_ylim(40, 100)
        # subdata.plot(x='data_sondage', y='impossible', ax=ax2, legend=False, color='r')

    plt.tight_layout()
    plt.savefig(next_name(env["pd-work"]))
    plt.suptitle("IC5p")
    plt.show()

    number_and_match(env, df_total, hue_order, "Number of IC5p Match", "IC5p Match", "IC5p")
    number_and_match(env, df_total, hue_order, "Number of IC3p Match", "IC3p Match", "IC3p")
    number_and_match(env, df_total, hue_order, "Number of Comp Match", "Comp Match", "Comp")

    df_comprehensive = df_total.groupby(["Chunk Size", "Tool"], as_index=False).sum()
    df_comprehensive = pd.melt(df_comprehensive, id_vars=["Chunk Size", "Tool"],
                               value_vars=[f"Number of {x} Match" for x in ["IC3p", "IC5p", "Comp"]] + [
                                   "Number of Match"],
                               var_name="Partial", value_name="Value")

    df_comprehensive_2 = pd.melt(df_total.groupby(["Chunk Size", "Tool"], as_index=False).sum(),
                                 id_vars=["Chunk Size", "Tool"],
                                 value_vars=[f"Number of {x} Found" for x in ["IC3p", "IC5p", "Comp"]] + [
                                     "Number of Found"],
                                 var_name="Partial", value_name="Value")

    df_comprehensive["Match"] = 100 * df_comprehensive["Value"] / df_comprehensive_2["Value"]

    g = seaborn.lmplot("Chunk Size", "Match", data=df_comprehensive, hue="Tool", col="Partial", lowess=True)
    g.set(xlim=(0, 5010), ylim=(0, 100))
    plt.show()

    print(df_comprehensive.to_csv())

    return

    fig, axes = plt.subplots(2, 4, sharey="all", sharex="all")
    axes = axes.ravel()
    for i, g in enumerate(genomes):
        ax = axes[i]  # type: plt.Axes

        df_curr = df[df["Genome"] == g]
        df_curr = pd.melt(df_curr, id_vars=["Genome", "Chunk Size"],
                          value_vars=[x for x in df_curr.columns if "Number of Error(" in x],
                          var_name="Combination", value_name="Number of Error")

        seaborn.lineplot("Chunk Size", "Number of Error", data=df_curr, hue="Combination", ax=ax, legend=False)

    plt.show()
    fig, axes = plt.subplots(2, 4, sharey="all", sharex="all")
    axes = axes.ravel()

    for i, g in enumerate(genomes):
        ax = axes[i]  # type: plt.Axes
        df_curr = df[df["Genome"] == g]
        df_curr = pd.melt(df_curr, id_vars=["Genome", "Chunk Size"],
                          value_vars=[x for x in df_curr.columns if "Number of Found(" in x],
                          var_name="Combination", value_name="Number of Found")
        seaborn.lineplot("Chunk Size", "Number of Found", data=df_curr, hue="Combination", ax=ax, legend=False)

    plt.show()


def viz_plot_per_genome_5p(env, df_gcfid):
    # type: (Environment, pd.DataFrame) -> None
    pass


def viz_stats_genome_level(env, df_gcfid, tools, reference, **kwargs):
    # type: (Environment, pd.DataFrame, List[str], str, Dict[str, Any]) -> None

    # 3' analysis
    viz_plot_per_genome_y_error_x_chunk(env, df_gcfid)

    # 5' analysis
    viz_plot_per_genome_5p(env, df_gcfid)


def viz_stats_3p_number_of_predictions_number_of_found(env, df_tidy, reference):
    # type: (Environment, pd.DataFrame, str) -> None

    g = seaborn.FacetGrid(df_tidy, col="Genome", col_wrap=4, hue="Tool", sharey=True, palette=CM.get_map("tools"))

    g.map(plt.plot, "Chunk Size", "Number of Predictions", linestyle="dashed")
    g.map(plt.plot, "Chunk Size", "Number of Found")

    g.set_titles("{col_name}", style="italic")
    # g.set(ylim=(0, 1))
    g.set(xlim=(0, 5100))
    g.set_xlabels("Chunk Size (nt)")
    g.set_ylabels("Score")
    g.add_legend()
    plt.savefig(next_name(env["pd-work"]))
    plt.show()


def viz_stats_3p_sensitivity_specificity(env, df_tidy, reference):
    # type: (Environment, pd.DataFrame, str) -> None

    g = seaborn.FacetGrid(df_tidy, col="Genome", col_wrap=4, hue="Tool", sharey=True, palette=CM.get_map("tools"))

    g.map(plt.plot, "Chunk Size", "Sensitivity")

    g.set_titles("{col_name}", style="italic")
    # g.set(ylim=(0, 1))
    g.set(xlim=(0, 5100))
    g.set_xlabels("Chunk Size (nt)")
    g.set_ylabels("Sensitivity")
    g.add_legend()
    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    g = seaborn.FacetGrid(df_tidy, col="Genome", col_wrap=4, hue="Tool", sharey=True, palette=CM.get_map("tools"))

    g.map(plt.plot, "Chunk Size", "Specificity")

    g.set_titles("{col_name}", style="italic")
    # g.set(ylim=(0, 1))
    g.set(xlim=(0, 5100))
    g.set_xlabels("Chunk Size (nt)")
    g.set_ylabels("Specificity")
    g.add_legend()
    plt.savefig(next_name(env["pd-work"]))
    plt.show()


def viz_stats_3p_number_of_predictions_precision(env, df_tidy, reference):
    # type: (Environment, pd.DataFrame, str) -> None
    df_tidy = df_tidy[df_tidy["Tool"].apply(lambda x: x.lower()) != reference.lower()]
    df_tidy = df_tidy[df_tidy["Tool"].apply(lambda x: x.lower()) != "mgm"]

    hue_order = sorted(df_tidy["Tool"].unique())
    cw = 4
    g = seaborn.FacetGrid(df_tidy, col="Genome", col_wrap=cw, hue="Tool", hue_order=hue_order,
                          sharey=True, palette=CM.get_map("tools"))

    g.map(plt.plot, "Chunk Size", "Number of Predictions", linestyle="dashed")
    g.map(plt.plot, "Chunk Size", "Precision")

    g.set_titles("{col_name}", style="italic")
    # g.set(ylim=(0, 1))
    g.set(xlim=(0, 5100))
    g.set_xlabels("Chunk Size (nt)")
    g.set_ylabels("Number of Predictions")
    g.add_legend()

    counter = 0
    for ax, (_, subdata) in zip(g.axes, df_tidy.groupby('Genome')):
        ax.yaxis.set_major_formatter(FuncFormatter(number_formatter))
        if counter == 0:
            yticklabels = ax.get_yticklabels()
        ax2 = ax.twinx()
        subdata = subdata.sort_values("Chunk Size")
        for hue in hue_order:
            subdata_hue = subdata[subdata["Tool"] == hue]
            ax2.plot(subdata_hue["Chunk Size"], subdata_hue["Precision"], label=hue,
                     color=CM.get_map("tools")[hue.lower()])
            ax2.set_ylim(0, 1)

        counter += 1
        if counter % cw == 0:
            ax2.set_ylabel("Precision")
        else:
            ax2.set_yticks([])

        # if counter % cw != 1:
        #     ax.set_yticklabels([])
        # else:
        #     ax.set_yticklabels(yticklabels)

    plt.tight_layout()
    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    # tabulate
    # df_piv = df_tidy.pivot(index="Genome", columns="Tool", values=["Precision"])
    # print(df_piv.to_csv())

    def f_mi(x):
        d = []
        d.append(x['a'].sum())
        d.append(x['a'].max())
        d.append(x['b'].mean())
        d.append((x['c'] * x['d']).sum())
        return pd.Series(d, index=[['a', 'a', 'b', 'c_d'],
                                   ['sum', 'max', 'mean', 'prodsum']])

    df1 = df_tidy.groupby(["Chunk Size", "Tool"], as_index=False).sum()
    # df_tidy.groupby(["Chunk Size", "Tool"], as_index=False).agg(
    #     {**{x: ['sum'] for x in df_tidy.columns if x not in  {"Chunk Size", "Tool", "Number in Reference"}},
    #      'Number in Reference': ['sum']})

    df1["Precision"] = df1["Number of Found"] / df1["Number of Predictions"]
    df1["WR"] = (df1["Number of Predictions"] - df1["Number of Found"]) / df1["Number of Predictions"]

    df1["Sensitivity"] = df1["Number of Found"] / df1["Number in Reference"]
    df1["Specificity"] = df1["Number of Found"] / df1["Number of Predictions"]
    print(df1.pivot(index="Chunk Size", columns="Tool", values=["Precision", "Number of Found"]))

    print(df1.pivot(index="Chunk Size", columns="Tool", values=["Precision"]))
    print(df1.pivot(index="Chunk Size", columns="Tool",
                    values=["Precision", "Number of Missed", "Number of Predictions"]))

    print(df1.pivot(index="Chunk Size", columns="Tool",
                    values=["Sensitivity", "Specificity", "Number of Found", "Number in Reference",
                            "Number of Predictions"]).to_csv())
    print("hi")

    df1 = df_tidy.groupby(["Chunk Size", "Tool"], as_index=False).mean()
    print(df1.pivot(index="Chunk Size", columns="Tool",
                    values=["Sensitivity", "Specificity"]).to_csv())


def viz_stats_5p_number_of_errors_number_of_found(env, df_tidy, reference):
    # type: (Environment, pd.DataFrame, str) -> None
    g = seaborn.FacetGrid(df_tidy, col="Genome", col_wrap=4, hue="Tool", sharey=True, palette=CM.get_map("tools"))

    g.map(plt.plot, "Chunk Size", "Number of Found", linestyle="dashed")
    g.map(plt.plot, "Chunk Size", "Number of Error")

    g.set_titles("{col_name}", style="italic")
    # g.set(ylim=(0, 1))
    g.set(xlim=(0, 5100))
    g.set_xlabels("Chunk Size (nt)")
    g.set_ylabels("Sensitivity")
    g.add_legend()
    plt.savefig(next_name(env["pd-work"]))
    plt.show()


def viz_stats_5p_error_rate(env, df_tidy, reference):
    # type: (Environment, pd.DataFrame, str) -> None
    g = seaborn.FacetGrid(df_tidy, col="Genome", col_wrap=4, hue="Tool", sharey=True, palette=CM.get_map("tools"))

    g.map(plt.plot, "Chunk Size", "Error Rate")

    g.set_titles("{col_name}", style="italic")
    # g.set(ylim=(0, 1))
    g.set(xlim=(0, 5100))
    g.set_xlabels("Chunk Size (nt)")
    g.set_ylabels("Error Rate")
    g.add_legend()
    plt.savefig(next_name(env["pd-work"]))
    plt.show()


def viz_stats_5p_error_rate_partial(env, df_tidy, reference):
    # type: (Environment, pd.DataFrame, str) -> None
    df_tidy = df_tidy[df_tidy["Tool"].apply(lambda x: x.lower()) != "verified"].copy()

    for cond in ["IC5p", "IC3p", "Comp"]:
        df_tidy[f"Error Rate {cond}"] = (df_tidy[f"Number of {cond} Found"] - df_tidy[f"Number of {cond} Match"]) / \
                                        df_tidy[f"Number of {cond} Found"]

        g = seaborn.FacetGrid(df_tidy, col="Genome", col_wrap=4, hue="Tool", sharey=True, palette=CM.get_map("tools"))

        g.map(plt.plot, "Chunk Size", f"Error Rate {cond}")

        g.set_titles("{col_name}", style="italic")
        # g.set(ylim=(0, 1))
        g.set(xlim=(0, 5100))
        g.set_xlabels("Chunk Size (nt)")
        g.set_ylabels("Error Rate")
        g.add_legend()
        plt.suptitle({
                         "IC5p": "Incomplete at 5' end",
                         "IC3p": "Incomplete at 3' end",
                         "Comp": "Complete genes"
                     }[cond])
        plt.savefig(next_name(env["pd-work"]))
        plt.show()

    # show 5p error by condition (combine all tools)
    df2 = df_tidy.groupby(["Chunk Size", "Tool"], as_index=False).sum()
    df2_tidy = pd.melt(
        df2, id_vars=["Chunk Size", "Tool"],
        value_vars=[f"Number of {cond} Found" for cond in ["IC5p", "IC3p", "Comp"]],
        var_name="Condition", value_name="Found"
    )
    df2_tidy["Condition"] = df2_tidy["Condition"].apply(lambda x: x.split()[2])
    df_tmp = pd.melt(
        df2, id_vars=["Chunk Size", "Tool"],
        value_vars=[f"Number of {cond} Match" for cond in ["IC5p", "IC3p", "Comp"]],
        var_name="Condition", value_name="Match"
    )
    df_tmp["Condition"] = df_tmp["Condition"].apply(lambda x: x.split()[2])

    df2_tidy = reduce(lambda df1, df2: pd.merge(df1, df2, on=["Chunk Size", "Condition", "Tool"],
                                                how="outer"), [df2_tidy, df_tmp])

    df2_tidy[f"Error Rate"] = (df2_tidy[f"Found"] - df2_tidy[f"Match"]) / df2_tidy[f"Found"]

    df2_tidy["Condition"].replace({
        "IC5p": "Incomplete at Gene Start",
        "IC3p": "Incomplete at Gene End",
        "Comp": "Complete genes"
    }, inplace=True)

    hue_order = sorted(df_tidy["Tool"].unique())
    g = seaborn.FacetGrid(df2_tidy, col="Condition", hue="Tool", sharey=True, palette=CM.get_map("tools"),
                          hue_order=hue_order)

    g.map(plt.plot, "Chunk Size", f"Error Rate")
    g.set_titles("{col_name}", style="italic")
    # g.set(ylim=(0, 1))
    g.set(xlim=(0, 5100))
    g.set_xlabels("Chunk Size (nt)")
    g.set_ylabels("Gene-Start Error Rate")
    g.add_legend()

    for ax, (_, subdata) in zip(g.axes[0], df2_tidy.groupby('Condition')):
        ax2 = ax.twinx()
        subdata = subdata.sort_values("Chunk Size")
        for hue in hue_order:
            subdata_hue = subdata[subdata["Tool"] == hue]
            ax2.plot(subdata_hue["Chunk Size"], subdata_hue["Found"], label=hue, linestyle="dashed")
            # ax2.set_ylim(40, 100)
        # subdata.plot(x='data_sondage', y='impossible', ax=ax2, legend=False, color='r')

    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    ###################### 2-level facetgrid ######################

    # show 5p error by condition (combine all tools)
    df2 = df_tidy.groupby(["Chunk Size", "Tool"], as_index=False).sum()
    df2["Number of Comp Found"] += df2["Number of IC3p Found"]
    df2["Number of Comp Match"] += df2["Number of IC3p Match"]

    df2_tidy = pd.melt(
        df2, id_vars=["Chunk Size", "Tool"],
        value_vars=[f"Number of {cond} Found" for cond in ["IC5p", "Comp"]],
        var_name="Condition", value_name="Score"
    )
    df2_tidy["Condition"] = df2_tidy["Condition"].apply(lambda x: x.split()[2])
    df2_tidy["Metric"] = "Found"

    df_tmp = pd.melt(
        df2, id_vars=["Chunk Size", "Tool"],
        value_vars=[f"Number of {cond} Match" for cond in ["IC5p", "Comp"]],
        var_name="Condition", value_name="Match"
    )
    df_tmp["Condition"] = df_tmp["Condition"].apply(lambda x: x.split()[2])
    df_tmp = reduce(lambda df1, df2: pd.merge(df1, df2, on=["Chunk Size", "Condition", "Tool"],
                                              how="outer"), [df2_tidy, df_tmp])

    df_tmp[f"Score"] = (df_tmp[f"Score"] - df_tmp[f"Match"]) / df_tmp[f"Score"]
    df_tmp["Metric"] = "Error Rate"

    df2_tidy = pd.concat([df2_tidy, df_tmp])

    df2_tidy["Condition"].replace({
        "IC5p": "Incomplete at Gene Start",
        # "IC3p": "Incomplete at Gene End",
        "Comp": "Complete at Gene Start"
    }, inplace=True)
    df2_tidy = df2_tidy[df2_tidy["Chunk Size"] <= 5000]

    hue_order = sorted(df2_tidy["Tool"].unique())
    g = seaborn.FacetGrid(
        df2_tidy, col="Condition", hue="Tool", sharey="row", palette=CM.get_map("tools"),
        row="Metric", hue_order=hue_order
    )

    g.map(plt.plot, "Chunk Size", f"Score")
    g.set_titles("{col_name}", style="italic")
    # g.set(ylim=(0, 1))
    # g.set(xlim=(0, 5100))
    g.set_xlabels("Chunk Size (nt)")
    # g.set_ylabels("Gene-Start Error Rate")

    for i, axes_row in enumerate(g.axes):
        for j, axes_col in enumerate(axes_row):

            if j == 0:
                if i == 0:
                    axes_col.set_ylabel("Number of Genes Found")
                else:
                    axes_col.set_ylabel("Gene-Start Error Rate")

    g.add_legend()
    plt.tight_layout(rect=[0, 0, 0.8, 1])
    plt.savefig(next_name(env["pd-work"]))
    plt.show()


def _helper_join_reference_and_tidy_data(env, df_per_gene, tools, list_ref):
    # type: (Environment, pd.DataFrame, List[str], List[str]) -> [str, pd.DataFrame]
    reference = _helper_df_joint_reference(df_per_gene, list_ref)
    df_per_gene = update_dataframe_with_stats(df_per_gene, tools, reference).copy()

    #### Genome Level
    # compute stats per genome
    df_stats_gcfid = list()
    for _, df_group in df_per_gene.groupby("Chunk Size", as_index=False):
        df_stats_gcfid.append(get_stats_at_gcfid_level_with_reference(df_group, tools, reference))
    df_per_genome = pd.concat(df_stats_gcfid, ignore_index=True, sort=False)

    df_tidy = tidy_genome_level(env, df_per_genome)
    df_tidy = df_tidy[df_tidy["Tool"].apply(lambda x: x.lower()).isin(tools + [reference])]

    return reference, df_tidy


def viz_stats_3p_missed_vs_length(env, df_per_gene, reference):
    # type: (Environment, pd.DataFrame, str) -> None

    pass



def viz_stats_3p(env, df_per_gene, tools, list_ref):
    # type: (Environment, pd.DataFrame, List[str], List[str]) -> None
    """Visualize statistics at 3prime level"""
    reference, df_tidy = _helper_join_reference_and_tidy_data(env, df_per_gene, tools, list_ref)

    ########## Genome Level ##########
    # Number of Predictions, number of found
    viz_stats_3p_number_of_predictions_number_of_found(env, df_tidy, reference)

    # Number of Predictions, Precision
    viz_stats_3p_number_of_predictions_precision(env, df_tidy, reference)

    # Sensitivity Specificity
    viz_stats_3p_sensitivity_specificity(env, df_tidy, reference)


    ########## Gene Level ##########

    # Missed vs reference length
    viz_stats_3p_missed_vs_length(env, df_per_gene, reference)



def viz_stats_5p(env, df_per_gene, tools, list_ref):
    # type: (Environment, pd.DataFrame, List[str], List[str]) -> None
    """Visualize statistics at 5prime level"""
    reference, df_tidy = _helper_join_reference_and_tidy_data(env, df_per_gene, tools, list_ref)

    # Number of 5p Errors, number of found
    viz_stats_5p_number_of_errors_number_of_found(env, df_tidy, reference)
    viz_stats_5p_error_rate(env, df_tidy, reference)
    viz_stats_5p_error_rate_partial(env, df_tidy, reference)


def viz_stats_per_gene(env, df_per_gene, tools, list_ref_5p, list_ref_3p):
    # type: (Environment, pd.DataFrame, List[str], List[str], List[str]) -> None

    viz_stats_3p(env, df_per_gene, tools, list_ref_3p)
    viz_stats_5p(env, df_per_gene, tools, list_ref_5p)


def tools_match_for_dataframe_row(r, tools):
    # type: (pd.Series, Iterable[str]) -> bool

    # check all tools make a prediction for current gene
    list_5ps = list()

    for t in tools:
        if r[f"5p-{t}"] is None:
            return False
        list_5ps.append(r[f"5p-{t}"])

    return all_elements_equal(list_5ps)


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

    viz_stats_per_gene(env, df, tools, args.ref_5p, args.ref_3p)


if __name__ == "__main__":
    main(my_env, parsed_args)
