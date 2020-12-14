# Author: Karl Gemayel
# Created: 9/4/20, 10:28 AM

import logging
import argparse
import os

import pandas as pd
from typing import *

import seaborn
from Bio import SeqIO
import matplotlib.pyplot as plt

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.genome_list import GenomeInfoList
from mg_viz.general import square_subplots, set_size

from mg_container.gms2_mod import GMS2Mod
from mg_general import Environment, add_env_args_to_parser
from mg_general.general import get_value, os_join, next_name
from mg_general.shelf import compute_gc
from mg_io.general import save_obj, load_obj
from mg_options.parallelization import ParallelizationOptions
from mg_parallelization.pbs import PBS
from mg_pbs_data.mergers import merge_identity
from mg_pbs_data.splitters import split_gil
from mg_argparse.parallelization import add_parallelization_options



# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #


parser = argparse.ArgumentParser("Plot the distribution of GMS2 groups across GC.")

parser.add_argument('--pf-gil', required=True)
parser.add_argument('--pf-checkpoint')
add_parallelization_options(parser)

parser.add_argument('--dn-gms2', default="gms2", required=False)


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


def helper_read_group_data(env, gil, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, Any]) -> pd.DataFrame
    dn_gms2 = get_value(kwargs, "dn_gms2", "gms2")

    list_entries = list()

    for gi in gil:
        try:
            dn_gms2 = "gms2" if int(gi.genetic_code) == 11 else "gms2_known"

            pf_mod = os_join(env["pd-runs"], gi.name, dn_gms2, "GMS2.mod")
            mod = GMS2Mod.init_from_file(pf_mod)

            gc = compute_gc(SeqIO.to_dict(SeqIO.parse(
                    os_join(env["pd-data"], gi.name, "sequence.fasta"), "fasta"))
            )

            group = mod.items["GENOME_TYPE"].split("-")[1].upper()

            list_entries.append({
                "Genome": gi.name,
                "Name": gi.attributes.get("name"),
                "Ancestor": gi.attributes.get("ancestor"),
                "GC": gc,
                "Group": group,
                "Genetic Code": gi.genetic_code
            })
        except FileNotFoundError:
            continue

    return pd.DataFrame(list_entries)

def read_group_data(env, gil, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, Any]) -> pd.DataFrame

    prl_options = get_value(kwargs, "prl_options", None)

    if not prl_options or not prl_options["use-pbs"]:
        df = helper_read_group_data(env, gil, **kwargs)
    else:
        pbs = PBS(env, prl_options, splitter=split_gil, merger=merge_identity)
        list_results = pbs.run(
            gil,
            helper_read_group_data,
            {"env": env, **kwargs}
        )

        df = pd.concat(list_results, ignore_index=True, sort=False)

    return df


def viz_gms2_groups_over_gc(env, df):
    # type: (Environment, pd.DataFrame) -> None

    df = df.sort_values(["Group", "GC"], ignore_index=True).copy()
    df["Genetic Code"] = df["Genetic Code"].astype(int)
    df[df["Group"] == "D2"] = "D"
    indeces = (df["Ancestor"] == "Archaea") & (df["Genetic Code"] == 11)
    #df.loc[indeces, "Group"] = df[indeces]["Group"].apply(lambda x: f"{x}*")

    figsize = set_size("thesis", subplots=(2, 3), titles=True)
    #figsize = (figsize[0], figsize[1]*2)

    fig, axes = plt.subplots(1, 2, sharex="all", sharey="all", figsize=figsize)

    for i, gcode in enumerate([4, 11]):
        ax = axes[i]
        curr_df = df[df["Genetic Code"] == gcode]

        ax.set_title(f"Genetic Code {gcode}")

        # manual filtering
        if gcode == 4:
            curr_df = curr_df[curr_df["Group"].isin({"A", "C"})]

        for g in sorted(curr_df["Group"].unique()):
            df_group = curr_df[curr_df["Group"] == g]
            seaborn.kdeplot(df_group["GC"], label=f"{g} ({len(df_group)})", ax=ax,
                    color={"C": "green", "A": "blue", "B": "orange", "D": "purple", "A*": "gray", "X": "red"}[g])

        if i == 0:
            ax.set_ylabel("Density")
    plt.xlim(10, 90)
    plt.tight_layout()
    # handles, labels = ax.get_legend_handles_labels()

    # leg = fig.legend(handles, labels, bbox_to_anchor=(0.5, 0.1), loc='upper center', ncol=2,
    #                  bbox_transform=fig.transFigure, frameon=False)
    #                  # fontsize=fontsize)

    # for lh in leg.legendHandles:
    #     lh0
    plt.savefig(next_name(env["pd-work"]))
    plt.show()


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    if not args.pf_checkpoint or not os.path.isfile(args.pf_checkpoint):
        gil = GenomeInfoList.init_from_file(args.pf_gil)
        prl_options = ParallelizationOptions.init_from_dict(env, args.pf_parallelization_options, vars(args))
        df = read_group_data(env, gil, dn_gms2=args.dn_gms2, prl_options=prl_options)

        if args.pf_checkpoint:
            save_obj(df, args.pf_checkpoint)
    else:
        df = load_obj(args.pf_checkpoint)

    viz_gms2_groups_over_gc(env, df)

if __name__ == "__main__":
    main(my_env, parsed_args)
