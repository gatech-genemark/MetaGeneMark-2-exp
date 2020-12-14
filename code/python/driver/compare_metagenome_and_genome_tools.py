# Author: Karl Gemayel
# Created: 8/14/20, 10:26 AM

import logging
import argparse
import os

import pandas as pd
from typing import *

import matplotlib.pyplot as plt
import seaborn
from Bio import SeqIO

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.genome_list import GenomeInfoList, GenomeInfo
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import os_join, next_name
from mg_general.labels_comparison_detailed import LabelsComparisonDetailed
from mg_general.shelf import compute_gc
from mg_io.general import save_obj, load_obj
from mg_io.labels import read_labels_from_file
from mg_options.parallelization import ParallelizationOptions
from mg_parallelization.generic_threading import run_n_per_thread
from mg_parallelization.pbs import PBS
from mg_pbs_data.mergers import merge_identity
from mg_pbs_data.splitters import split_gil

parser = argparse.ArgumentParser("Compare Prodigal vs MetaProdigal, and GMS2 vs MGM2.")

parser.add_argument('--pf-gil', required=True, help="List of genomes")
parser.add_argument('--pf-parallelization-options')
parser.add_argument('--pf-checkpoint')

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


def compare_metagenome_and_standard_tool(env, gi, tool_meta, tool_standard, tag):
    # type: (Environment, GenomeInfo, str, str, str) -> pd.Series

    pf_meta = os_join(env["pd-runs"], gi.name, tool_meta, "prediction.gff")
    pf_standard = os_join(env["pd-runs"], gi.name, tool_standard, "prediction.gff")

    label_meta = read_labels_from_file(pf_meta)
    label_standard = read_labels_from_file(pf_standard)

    lcd = LabelsComparisonDetailed(label_standard, label_meta)

    error = 0
    if float(len(lcd.match_3p('a'))) > 0:
        error = len(lcd.match_3p_not_5p('a')) / float(len(lcd.match_3p('a')))

    return pd.Series({
        "Genome": gi.name,
        "Tag": tag,
        "Clade": gi.attributes.get("ancestor"),
        "Match 3p 5p": len(lcd.match_3p_5p('a')),
        "Match 3p": len(lcd.match_3p('a')),
        "Error Rate 5p": error,
        "Sensitivity": len(lcd.match_3p('a')) / len(label_standard) if len(label_standard) > 0 else 0,
        "Specificity": len(lcd.match_3p('a')) / len(label_meta) if len(label_meta) > 0 else 0
    })


def compare_for_gi(env, gi, **kwargs):
    # type: (Environment, GenomeInfo, Dict[str, Any]) -> pd.DataFrame

    pf_sequence = os_join(env["pd-data"], gi.name, "sequence.fasta")
    sequences = SeqIO.to_dict(SeqIO.parse(pf_sequence, "fasta"))
    genome_gc = compute_gc(sequences)

    se_prod = compare_metagenome_and_standard_tool(env, gi, "mprodigal", "prodigal", "Prodigal")
    se_gms2 = compare_metagenome_and_standard_tool(env, gi, "mgm2", "gms2", "GeneMarkS-2")

    se_prod["GC"] = genome_gc
    se_gms2["GC"] = genome_gc

    return pd.DataFrame([se_prod, se_gms2])

def compare_for_gil(env, gil, **kwargs):

    list_df = list()
    for gi in gil:
        list_df.append(compare_for_gi(env, gi, **kwargs))

    return pd.concat(list_df, ignore_index=True, sort=False)


def plot_5p_difference_versus_gc(env, df):
    # type: (Environment, pd.DataFrame) -> None
    df = df[df["Sensitivity"] > 0.4].copy()

    # g = seaborn.FacetGrid(data=df, col="Clade", hue="Tag", sharex=True)
    # g.map(seaborn.regplot, )

    scatter_kws={"s": 2, "alpha": 0.3}

    seaborn.lmplot("GC", "Error Rate 5p", data=df,
                   hue="Tag", col="Clade",
                   lowess=True,
                   scatter_kws=scatter_kws)

    plt.show()

    df_tidy = pd.melt(df, id_vars=["Tag", "Clade", "GC"],
                      value_vars=["Sensitivity", "Specificity"],
                      var_name="Metric", value_name="Score")
    g = seaborn.lmplot("GC", "Score", data=df_tidy,
                   hue="Tag", col="Metric", lowess=True,
                   scatter_kws=scatter_kws)
    g.set(ylim=[0.8, 1.001])

    plt.show()

    # print sensitivity less than 0.4
    print(df[df["Sensitivity"] < 0.4]["Genome"].unique())





def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    pf_checkpoint = args.pf_checkpoint

    prl_options = ParallelizationOptions.init_from_dict(
        env, args.pf_parallelization_options, vars(args)
    )

    if not pf_checkpoint or not os.path.isfile(pf_checkpoint):
        gil = GenomeInfoList.init_from_file(args.pf_gil)

        if not prl_options["use-pbs"]:
            df_list = run_n_per_thread(
                [g for g in gil],
                compare_for_gi,
                data_arg_name="gi",
                func_kwargs={"env": env}
            )
        else:
            pbs = PBS(env, prl_options,
                      splitter=split_gil,
                      merger=merge_identity)
            df_list = pbs.run(gil, compare_for_gil, {"env": env})

        df = pd.concat(df_list, ignore_index=True, sort=False)

        if pf_checkpoint:
            save_obj(df, pf_checkpoint)
    else:
        df = load_obj(pf_checkpoint)

    df.to_csv(next_name(env["pd-work"], ext="csv"), index=False)

    plot_5p_difference_versus_gc(env, df)



if __name__ == "__main__":
    main(my_env, parsed_args)
