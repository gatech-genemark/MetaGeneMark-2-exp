# Author: Karl Gemayel
# Created: 6/25/20, 4:11 PM
import copy
import re
import logging
import argparse
from tempfile import mkstemp

import numpy as np
import pandas as pd
from typing import *

from Bio import SeqIO
from tqdm import tqdm
import matplotlib.pyplot as plt

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.genome_list import GenomeInfoList, GenomeInfo
from mg_container.gms2_mod import GMS2Mod
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import os_join, next_name, fix_names, get_value
from mg_general.labels import Labels
from mg_general.labels_comparison_detailed import LabelsComparisonDetailed
from mg_io.general import mkdir_p, write_to_file, remove_p, save_obj, load_obj
from mg_io.labels import read_labels_from_file, read_lst
from mg_io.shelf import convert_multi_fasta_into_single_fasta
from mg_models.gms2_noncoding import GMS2Noncoding
from mg_models.motif_model import MotifModel
from mg_models.shelf import run_gms2_prediction_with_model, train_gms2_model, relative_entropy
from mg_parallelization.generic_threading import run_one_per_thread, run_n_per_thread
from mg_viz import sns
from mg_viz.general import FigureOptions

parser = argparse.ArgumentParser("Test RBS model if built by highest scoring motifs.")

parser.add_argument('--pf-gil', required=True, help="Genome info list")
parser.add_argument('--pf-trained')

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

def short_name(genome):
    # type: (str) -> str
    return f"{genome[0]}. {genome.split('_')[1]}"

def train_predict_score(env, gi, mod_original, pf_sequence, labels_training, labels_reference):
    # type: (Environment, GenomeInfo, GMS2Mod, str, Labels, Labels) -> Dict[str, Any]

    _, pf_labels = mkstemp()
    _, pf_mod = mkstemp()
    _, pf_prediction = mkstemp()

    # write labels to file
    write_to_file(labels_training.to_string_lst(shift_coordinates_by=0), pf_labels)

    # train and put in new RBS in original model
    mod = train_gms2_model(env, pf_sequence, pf_labels, pf_mod)
    new_mod = copy.deepcopy(mod_original)
    new_mod.items["RBS_MAT"] = mod.items["RBS_MAT"]
    new_mod.items["RBS_POS_DISTR"] = mod.items["RBS_POS_DISTR"]

    new_mod.to_file(pf_mod)

    # run prediction and compute error
    run_gms2_prediction_with_model(env, pf_sequence, pf_mod, pf_prediction)

    # compare with reference
    lcd = LabelsComparisonDetailed(labels_reference, read_labels_from_file(pf_prediction, shift=0))

    # cleanup
    remove_p(pf_labels, pf_mod, pf_prediction)

    return {
        "Mod": mod,
        "Error": 100 - 100 * len(lcd.match_3p_5p('a')) / len(lcd.match_3p('a'))
    }


def test_rbs_built_from_most_conserved_motifs_for_gi(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> Dict[str, Any]
    pd_work = env["pd-work"]
    mkdir_p(pd_work)

    pf_prediction = os_join(pd_work, "gms2.ext")
    pf_sequence = os_join(env['pd-data'], gi.name, "sequence.fasta")
    pf_verified = os_join(env['pd-data'], gi.name, "verified.gff")
    pf_mod = os_join(env['pd-runs'], gi.name, "gms2", "GMS2.mod")

    # convert all labels and sequences to single fasta (for GMS2 training)
    sequences, labels_verified = convert_multi_fasta_into_single_fasta(
        SeqIO.to_dict(SeqIO.parse(pf_sequence, "fasta")),
        read_labels_from_file(pf_verified, shift=0),
        "anydf"
    )

    # write sequence and labels to files
    f_lab, pf_verified_lst = mkstemp(suffix=".gff")
    f_seq, pf_sequence_single = mkstemp(suffix=".fasta")
    SeqIO.write(sequences.values(), open(pf_sequence_single, "w"), "fasta")
    write_to_file(labels_verified.to_string_lst(shift_coordinates_by=0), pf_verified_lst)


    # run gms2 with extended format to get RBS scores
    run_gms2_prediction_with_model(env, pf_sequence_single, pf_mod, pf_prediction, format="ext")

    labels = read_labels_from_file(pf_prediction, key_value_delimiter=" ", shift=0)
    labels_verified = read_lst(pf_verified_lst, shift=0)
    scores = [float(l.get_attribute_value("rbs_score")) for l in labels]

    sns.distplot(pd.DataFrame({"Scores": scores}), "Scores", figure_options=FigureOptions(
        ylabel="PDF", title=f"Distribution of motif scores in {short_name(gi.name)}", xlabel="Motif Score",
        save_fig=next_name(env["pd-work"])
    ))
    # import sys
    # sys.exit()
    sorted_labels = sorted(labels, key=lambda l: float(l.get_attribute_value("rbs_score")))

    result = dict()
    current_index = 0  # index from which labels are selected (above)

    gms2_mod = GMS2Mod.init_from_file(pf_mod)

    list_entries = list()
    for thresh in tqdm(np.arange(-2.0, 6, 0.2)):

        while current_index < len(labels) and float(
                sorted_labels[current_index].get_attribute_value("rbs_score")) < thresh:
            current_index += 1

        if current_index >= len(labels):
            break

        label_subset_upper = Labels(sorted_labels[current_index:])
        score_upper = train_predict_score(env, gi, gms2_mod, pf_sequence_single, label_subset_upper, labels_verified)

        label_subset_lower = Labels(sorted_labels[:current_index])
        score_lower = train_predict_score(env, gi, gms2_mod, pf_sequence_single, label_subset_lower, labels_verified)

        list_entries.append({
            "Threshold": thresh,
            "Genome": gi.name,
            "Upper Percentage of labels": 100 * len(label_subset_upper) / len(labels),
            **{f"Upper {x}": score_upper[x] for x in score_upper.keys()},
            "Lower Percentage of labels": 100 * len(label_subset_lower) / len(labels),
            **{f"Lower {x}": score_lower[x] for x in score_lower.keys()}
        })

    df = pd.DataFrame(list_entries)

    remove_p(pf_sequence_single, pf_verified_lst)

    return {
        "df": df,
        "Scores": scores,
    }


def test_rbs_built_from_most_conserved_motifs(env, gil, **kwargs):
    # type: (Environment, GenomeInfoList, Dict[str, Any]) -> None

    # list_results = run_n_per_thread(
    #     [x for x in gil],
    #     test_rbs_built_from_most_conserved_motifs_for_gi,
    #     "gi",
    #     {"env": env, **kwargs},
    #     n=1
    # )
    pf_trained = get_value(kwargs, "pf_trained", None)


    if pf_trained is None:
        list_df = list()
        for gi in gil:
            # if "denitrificans" not in gi.name:
            #     continue
            result = test_rbs_built_from_most_conserved_motifs_for_gi(
                env.duplicate({"pd-work": os_join(env["pd-work"], gi.name)}),
                gi
            )
        # for result, gi in zip(list_results, gil):
            list_df.append(result["df"])

            sns.distplot(pd.DataFrame({"Scores": result["Scores"]}), "Scores", figure_options=FigureOptions(
                ylabel="Frequency",
                save_fig=next_name(env["pd-work"])
            ))

        df = pd.concat(list_df, ignore_index=True, sort=True)
        df["Genome"] = df.apply(fix_names, axis=1)

        save_obj(df, "summary.pkl")
    else:
        df = load_obj(pf_trained)


    # sns.lmplot(df, "Threshold", "Percentage of labels", hue="Genome", sns_kwargs={"lowess": True},
    #            figure_options=FigureOptions(save_fig=next_name(env["pd-work"])))
    # sns.lmplot(df, "Threshold", "Error", hue="Genome", sns_kwargs={"lowess": True},
    #            figure_options=FigureOptions(save_fig=next_name(env["pd-work"])))
    # sns.lmplot(df, "Threshold", "Error", hue="Genome", sns_kwargs={"lowess": True, "scatter_kws": {"s": 0}},
    #            figure_options=FigureOptions(
    #                ylim=[0, 15], save_fig=next_name(env["pd-work"]))
    #            )
    # sns.lmplot(df, "Percentage of labels", "Error", hue="Genome", sns_kwargs={"lowess": True},
    #            figure_options=FigureOptions(save_fig=next_name(env["pd-work"])))


    for genome, df_group in df.groupby("Genome", as_index=False):

        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        for state, ax in zip(["Lower", "Upper"], axes.ravel()):
            df_group = df_group.copy()


            list_mm = [MotifModel(r[f"{state} Mod"].items["RBS_MAT"]) for _, r in df_group.iterrows()]
            list_non = [GMS2Noncoding(r[f"{state} Mod"].items["NON_MAT"]) for _, r in df_group.iterrows()]
            df_group["Relative Entropy"] = [relative_entropy(m, n) for m, n in zip(list_mm, list_non)]
            import seaborn
            seaborn.regplot("Threshold", f"{state} Error", df_group, ax=ax, lowess=True, scatter_kws={"s": 2},
                            label="Error")
            ax.set_ylim(0, 20)
            ax.set_xlabel(f"Threshold: {'Minimum' if state == 'Upper' else 'Maximum'} motif score")
            ax.set_ylabel("Gene-start error rate")
            ax2 = ax.twinx()
            seaborn.regplot("Threshold", "Relative Entropy", df_group, ax=ax2,
                            lowess=True, scatter_kws={"s": 2},color="red",
                            label="Relative Entropy")
            ax2.set_ylim(0, 12)
            # plt.legend(handles=[ax.lines, ax2.lines], labels=["Error", "Relative Entropy"])
            # added these three lines
            lines, labels = ax.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()
            ax2.legend(lines + lines2, labels + labels2, loc=0)

            ax.set_title(state)
        fig.suptitle(genome)

        plt.savefig(next_name(env["pd-work"]))
        # plt.title(genome)
        # labs = [l.get_label() for l in lns]
        # ax.legend(lns, labs, loc=0)
        plt.show()

    # list_df = [r["df"] for r in list_results]
    # df = pd.concat(list_df, ignore_index=True, sort=False)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    gil = GenomeInfoList.init_from_file(args.pf_gil)
    test_rbs_built_from_most_conserved_motifs(env, gil, pf_trained=args.pf_trained)


if __name__ == "__main__":
    main(my_env, parsed_args)
