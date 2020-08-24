# Author: Karl Gemayel
# Created: 8/22/20, 3:08 PM

import logging
import argparse
import os

import seaborn
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO
import logomaker as lm
from tqdm import tqdm

from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.genome_list import GenomeInfoList, GenomeInfo
from mg_container.gms2_mod import GMS2Mod
from mg_general.general import os_join, get_value, next_name
from mg_general.shelf import compute_gc
from mg_general import Environment, add_env_args_to_parser
from mg_io.general import save_obj, load_obj

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_models.gms2_noncoding import GMS2Noncoding
from mg_models.motif_model import MotifModel
from mg_viz.general import square_subplots, set_size

parser = argparse.ArgumentParser("Visualize GMS2 models over GC.")

parser.add_argument('--pf-gil', required=True)
parser.add_argument('--pf-checkpoint')

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


def viz_genome_type_per_gc(env, list_gi, list_mod, list_gc):
    # type: (Environment, List[GenomeInfo], List[GMS2Mod], List[float]) -> None
    df = pd.DataFrame({
        "GC": list_gc,
        "Group": [mod.items["GENOME_TYPE"] for mod in list_mod],
        "Mod": list_mod,
        "Name": [gi.attributes.get("name") for gi in list_gi]
    })
    df = df.sort_values(["Group", "GC"], ignore_index=True)

    df2 = df[df["Group"].isin({"group-a", "group-c"})]
    df2["Group"] = df2["Group"].apply(lambda x: x.split('-')[1].upper())
    seaborn.catplot("GC", "Group", data=df2 )
    plt.show()

    # kde plot
    for g in sorted(df["Group"].unique()):
        seaborn.kdeplot(df[df["Group"] == g]["GC"], label=g)
    plt.show()

    df_c = df[df["Group"] == "group-c"]

    # remove manually
    df_c = df_c.iloc[[0,7,10, 22]]
    num_c = len(df_c)
    num_rows, num_columns = square_subplots(num_c)

    fig, axes = plt.subplots(num_rows, num_columns, sharey="all")
    fig_rbs, axes_rbs = plt.subplots(num_rows, num_columns, sharey="all")

    list_df = list()
    for i in range(num_rows):
        for j in range(num_columns):

            index = i*num_columns + j
            # if index in {1, 4, 5, 9, 13, 21, 25, 31, 12, 20}:
            #     continue


            ax = axes[i][j]

            ax.set_xticklabels([])

            if index >= num_c:
                break
            mod = df_c.at[df_c.index[i*num_columns + j], "Mod"]

            motif_mat = MotifModel(mod.items[f"PROMOTER_MAT"], mod.items[f"PROMOTER_POS_DISTR"])
            nonc_mat = GMS2Noncoding(mod.items[f"NON_MAT"])

            list_df.append(pd.DataFrame({
                "Distance": range(len(motif_mat._spacer)),
                "Probability": motif_mat._spacer,
                "Hue": [df_c.at[df_c.index[index], "GC"]] * len(motif_mat._spacer)
            }))

            df_motif_mat = motif_mat.pwm_to_df()
            np_nonc_mat = nonc_mat.pwm_to_array(0)

            df_rel = lm.transform_matrix(
                df_motif_mat, from_type="probability", to_type="information", background=np_nonc_mat
            )

            lm.Logo(df_rel, ax=ax)
            ax.set_ylim(0, 2.5)

            ax = axes_rbs[i][j]
            motif_mat = MotifModel(mod.items[f"RBS_MAT"], mod.items[f"RBS_POS_DISTR"])
            nonc_mat = GMS2Noncoding(mod.items[f"NON_MAT"])
            df_motif_mat = motif_mat.pwm_to_df()
            np_nonc_mat = nonc_mat.pwm_to_array(0)

            df_rel = lm.transform_matrix(
                df_motif_mat, from_type="probability", to_type="information", background=np_nonc_mat
            )

            lm.Logo(df_rel, ax=ax)
            ax.set_ylim(0, 2.5)

    fig.show()
    fig_rbs.show()

    df = pd.concat(list_df, ignore_index=True)
    seaborn.lineplot("Distance", "Probability", data=df, hue="Hue")
    plt.show()


    figsize = set_size("thesis", subplots=(3, num_c), titles=True)
    fig, axes = plt.subplots(3, num_c, sharex="row", sharey="row", figsize=figsize)
    label_fs = "small"
    for i, idx in enumerate(df_c.index):

        mod = df_c.at[df_c.index[i], "Mod"]

        # get models
        rbs_model = MotifModel(mod.items[f"RBS_MAT"], mod.items[f"RBS_POS_DISTR"])
        prom_model = MotifModel(mod.items[f"PROMOTER_MAT"], mod.items[f"PROMOTER_POS_DISTR"])
        nonc_mat = GMS2Noncoding(mod.items[f"NON_MAT"])
        np_nonc_mat = nonc_mat.pwm_to_array(0)

        # rbs
        ax = axes[0][i]
        name = df_c.at[idx, "Name"].split()
        name = f"{name[0][0]}. {name[1]}"
        ax.set_title(name, fontsize="small")
        df_motif_mat = rbs_model.pwm_to_df()
        df_rel = lm.transform_matrix(
            df_motif_mat, from_type="probability", to_type="information", background=np_nonc_mat
        )
        lm.Logo(df_rel, ax=ax)
        ax.set_ylim(0, 2.5)

        ax.set_xticklabels([])
        ax.set_xticks([])
        if i == 0:
            ax.set_ylabel("Relative Entropy", fontsize=label_fs)

        # promoter
        ax = axes[1][i]
        df_motif_mat = prom_model.pwm_to_df()
        df_rel = lm.transform_matrix(
            df_motif_mat, from_type="probability", to_type="information", background=np_nonc_mat
        )
        lm.Logo(df_rel, ax=ax)
        ax.set_ylim(0, 2.5)

        ax.set_xticklabels([])
        ax.set_xticks([])
        # ax.set_yticklabels([])
        if i == 0:
            ax.set_ylabel("Relative Entropy", fontsize=label_fs)

        # spacers
        ax = axes[2][i]
        df_sp = pd.DataFrame({
                "Distance": list(range(len(rbs_model._spacer))) + list(range(len(prom_model._spacer))),
                "Probability": rbs_model._spacer + prom_model._spacer,
                "Type": ["RBS"] * len(rbs_model._spacer) + ["PROMOTER"] * len(prom_model._spacer)
        })

        for t in df_sp["Type"].unique():
            df_curr = df_sp[df_sp["Type"] == t]
            ax.plot(df_curr["Distance"], df_curr["Probability"], label=t)
        # seaborn.lineplot("Distance", "Probability", data=df_sp, hue="Type", ax=ax)
        if i == 0:
            ax.set_ylabel("Frequency", fontsize=label_fs)
        ax.set_xlabel("Spacer Length", fontsize=label_fs)

    ax = axes[2,  num_c-1]
    fig.subplots_adjust(bottom=0.2)
    handles, labels = ax.get_legend_handles_labels()

    leg = fig.legend(handles, labels, bbox_to_anchor=(0.5, 0.1), loc='upper center', ncol=2,
                     bbox_transform=fig.transFigure, frameon=False)
                     # fontsize=fontsize)

    for lh in leg.legendHandles:
        lh.set_alpha(1)
        # lh.set_sizes([18] * 2)

    for i in range(num_c):
        fig.align_ylabels(axes[:, i])

    fig.tight_layout(rect=[0, 0.1, 1, 1])
    fig.savefig(next_name(env["pd-work"]), bbox_extra_artists=(leg,))  # bbox_inches='tight'

    fig.show()




def viz_gms2_models_over_gc(env, list_gi, list_mod, list_gc, **kwargs):
    # type: (Environment, List[GenomeInfo], List[GMS2Mod], List[float], Dict[str, Any]) -> None
    viz_genome_type_per_gc(env, list_gi, list_mod, list_gc)


def read_genome_data(env, gil, **kwargs):
    dn_gms2 = get_value(kwargs, "dn_gms2", "gms2")

    list_gi = list()
    list_mod = list()
    list_gc = list()

    for gi in tqdm(gil, total=len(gil)):
        try:
            pf_mod = os_join(env["pd-runs"], gi.name, dn_gms2, "GMS2.mod")
            mod = GMS2Mod.init_from_file(pf_mod)

            list_mod.append(mod)
            list_gi.append(gi)
            list_gc.append(
                compute_gc(SeqIO.to_dict(SeqIO.parse(
                    os_join(env["pd-data"], gi.name, "sequence.fasta"), "fasta")))
            )
        except FileNotFoundError:
            continue

    return list_gi, list_mod, list_gc


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    if not args.pf_checkpoint or not os.path.isfile(args.pf_checkpoint):
        gil = GenomeInfoList.init_from_file(args.pf_gil)
        list_gi, list_mod, list_gc = read_genome_data(env, gil, dn_gms2=args.dn_gms2)

        if args.pf_checkpoint:
            save_obj([list_gi, list_mod, list_gc], args.pf_checkpoint)
    else:
        list_gi, list_mod, list_gc = load_obj(args.pf_checkpoint)

    viz_gms2_models_over_gc(env, list_gi, list_mod, list_gc)


if __name__ == "__main__":
    main(my_env, parsed_args)
