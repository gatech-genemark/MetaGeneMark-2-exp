# Author: Karl Gemayel
# Created: 7/8/20, 3:49 PM

import logging
import seaborn
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import fix_names, next_name

parser = argparse.ArgumentParser("DRIVER DESCRIPTION.")

parser.add_argument('--pf-stats', required=True)

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


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = pd.read_csv(args.pf_stats)


    df["Genome"] = df.apply(fix_names, axis=1)
    df = df.sort_values("Chunk Size").copy()
    df["Number of Errors"] = (df["Error"] * df["Number of genes found"] / 100.0).astype(int)
    # sns.lineplot(df_stats, "Chunk Size", "Error", hue="Genome")

    genes_found = {
        g: 0 for g in df["Genome"].unique()
    }

    genes_found = {
        df.at[idx, "Genome"]: max(genes_found[df.at[idx, "Genome"]], df.at[idx, "Number of genes found"])
        for idx in df.index
    }

    df["title"] = df.apply(lambda r: f"{r['Genome']} (from ({genes_found[r['Genome']]}))", axis=1)

    g = seaborn.FacetGrid(df, col="Genome", hue="Tool", col_wrap=4)
    g.map(plt.scatter, "Chunk Size", "Error", alpha=.7)
    g.add_legend()
    plt.show()

    g = seaborn.FacetGrid(df, col="Genome", hue="Tool", col_wrap=4)
    g.map(plt.plot, "Chunk Size", "Error", alpha=.7)
    g.add_legend()
    g.set_xlabels("Chunk Size (nt)")
    plt.xlim(0, 100000)

    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    g.set_titles("{col_name}")
    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    g = seaborn.FacetGrid(df, col="title", hue="Tool", col_wrap=4)
    g.map(plt.plot, "Chunk Size", "Number of Errors", alpha=.7)
    g.add_legend()
    g.set_xlabels("Chunk Size (nt)")
    plt.xlim(0, 100000)

    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    g.set_titles("{col_name}")
    plt.savefig(next_name(env["pd-work"]))
    plt.show()

    g = seaborn.FacetGrid(df, col="Genome", hue="Tool", col_wrap=4)
    g.map(plt.plot, "Chunk Size", "Number of genes found", alpha=.7)
    g.add_legend()
    g.set_xlabels("Chunk Size (nt)")
    plt.xlim(0, 100000)

    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    g.set_titles("{col_name}" )
    plt.savefig(next_name(env["pd-work"]))
    plt.show()


if __name__ == "__main__":
    main(my_env, parsed_args)
