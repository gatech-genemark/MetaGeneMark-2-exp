# Author: Karl Gemayel
# Created: 7/7/20, 10:12 AM

import logging
import argparse
import pandas as pd
from typing import *
import matplotlib.pyplot as plt

import seaborn
from tqdm import tqdm

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import os_join, fix_names
from mg_general.labels_comparison_detailed import LabelsComparisonDetailed
from mg_io.labels import read_labels_from_file
from mg_viz import sns

parser = argparse.ArgumentParser("Collect statistics for experiment running tools on "
                                 "genome chunks.")

parser.add_argument('--pf-summary', required=True)
parser.add_argument('--pf-output', required=True)

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


def stats_tools_on_chunks(env, df):
    # type: (Environment, pd.DataFrame) -> pd.DataFrame

    list_entries = list()

    for idx in tqdm(df.index, total=len(df)):
        pf_prediction = df.at[idx, "Predictions"]
        pf_verified = os_join(env["pd-data"], df.at[idx, "Genome"], "verified.gff")

        labels = read_labels_from_file(pf_prediction, shift=-1)
        labels_ref = read_labels_from_file(pf_verified)

        lcd = LabelsComparisonDetailed(labels_ref, labels)
        list_entries.append({
            "Error": 100 - 100 * len(lcd.match_3p_5p('a')) / len(lcd.match_3p('a')),
            "Number of genes found": len(lcd.match_3p('a')),
            **df.loc[idx, :].to_dict(),
        })

    return pd.DataFrame(list_entries)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = pd.read_csv(args.pf_summary)
    df_stats = stats_tools_on_chunks(env, df)
    df_stats.to_csv(args.pf_output, index=False)

    df_stats["Genome"] = df_stats.apply(fix_names, axis=1)

    sns.lineplot(df_stats, "Chunk Size", "Error", hue="Genome")


    g = seaborn.FacetGrid(df_stats, col="Genome", hue="Tool", col_wrap=4)
    g.map(plt.scatter, "Chunk Size", "Error", alpha=.7)
    g.add_legend()

    plt.show()

if __name__ == "__main__":
    main(my_env, parsed_args)
