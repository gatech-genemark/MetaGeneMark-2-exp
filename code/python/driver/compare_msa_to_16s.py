# Author: Karl Gemayel
# Created: 10/2/20, 10:28 AM

import logging
import argparse
import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.genome_list import GenomeInfoList
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import get_value, next_name, os_join
from mg_io.general import load_obj
from mg_models.building import build_mgm_motif_model_for_gc_v2
from mg_models.shelf import fix_genome_type, get_consensus_sequence, create_numpy_for_column_with_extended_motif
from mg_viz import sns
from mg_viz.general import FigureOptions

parser = argparse.ArgumentParser("Compare results from alignment to MSA.")


parser.add_argument('--pf-data', required=True, help="File containing GMS2 models")
parser.add_argument('--seq-16s', required=True, help="DNA-level 16S tail")

parser.add_argument('--gc-range', nargs=2, type=float, help="Range of allowed GC")
parser.add_argument('--group', choices=list("ABCDX"), type=str.upper, help="Allowed GMS2 group.")

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


def load_gms2_models_from_pickle(pf_mods):
    # type: (str) -> pd.DataFrame
    df = load_obj(pf_mods)  # type: pd.DataFrame

    df["Type"] = "Bacteria"

    df.reset_index(inplace=True)

    df["GENOME_TYPE"] = df.apply(lambda r: r["Mod"].items["GENOME_TYPE"], axis=1)
    # df["GC"] = df.apply(lambda r: r["Mod"].items["GC"], axis=1)
    fix_genome_type(df)

    return df


def best_match_to_16s(seq_16s, consensus):
    # type: (str, str) -> int

    bm_pos = 0
    bm_score = 0
    for curr_pos in range(len(seq_16s)-len(consensus)):

        curr_score = sum([1 for i in range(len(consensus)) if consensus[i] == seq_16s[curr_pos+i]])

        if curr_score > bm_score:
            bm_score = curr_score
            bm_pos = curr_pos

    return bm_pos

def best_shifts_to_16s(seq_16s, df, col):
    # type: (str, pd.DataFrame, str) -> None

    list_pos = list()
    for idx in df.index:
        consensus = df.at[idx, "CONSENSUS_RBS_MAT"]
        list_pos.append(best_match_to_16s(seq_16s, consensus))

    df[col] = list_pos
    df[col] = df[col] - min(df[col])


def compare_msa_to_16s(env, df, seq_16s, **kwargs):
    # type: (Environment, pd.DataFrame, str, Dict[str, any]) -> None
    split_by_gc = get_value(kwargs, "split_by_gc", False)
    gc_range = get_value(kwargs, "gc_range", [0, 100])
    group = get_value(kwargs, "group")

    name = "RBS_MAT"

    df[f"CONSENSUS_{name}"] = df.apply(lambda r: get_consensus_sequence(r["Mod"].items[name]), axis=1)

    # df = df[(df["GC"] > gc_range[0]) & (df["GC"] < gc_range[1])]
    print([x for x in df["CONSENSUS_RBS_MAT"].value_counts().items()])
    df = df[df["CONSENSUS_RBS_MAT"] == "AAAAAA"]
    # if group:
    #     df = df[df["GENOME_TYPE"] == group]

    df = df[df.groupby(f"CONSENSUS_{name}")[f"CONSENSUS_{name}"].transform(len) > 5]

    # RBS models from MSA
    array, update_shifts = create_numpy_for_column_with_extended_motif(env, df, "RBS_MAT")

    # RBS
    print(update_shifts)

    df["MSA"] = update_shifts
    best_shifts_to_16s(seq_16s, df, "16S")

    print(abs(df["MSA"] - df["16S"]).sum())


    # find peak distribution by shift
    shifts = sorted(df["MSA"].unique())

    def get_peak_pos(r):
        pos_distr = r["Mod"].items["RBS_POS_DISTR"] # type: Dict[str, float]

        return max(pos_distr, key=lambda key: pos_distr[key])

    for s in shifts:
        df_s = df[df["MSA"] == s].copy()

        df_s["Peak"] = df_s.apply(get_peak_pos, axis=1)
        consensus = df_s.loc[df_s.index[0], "CONSENSUS_RBS_MAT"]
        consensus = df_s["CONSENSUS_RBS_MAT"].mode().values[0]

        pos_to_count = {x[0]: x[1] for x in df_s["Peak"].value_counts().items()}

        for i in range(15):
            if i not in pos_to_count:
                pos_to_count[i] = 0

        sns.barplot(
            pd.DataFrame(
                {"Number of Models": [pos_to_count[k] for k in sorted(pos_to_count.keys())],
                "Peak of Spacer Distribution": sorted(pos_to_count.keys())}
            ),
            "Peak of Spacer Distribution", "Number of Models",
            figure_options=FigureOptions(
                title=f"Shift: {s}, Example Consensus: {consensus}",
                # save_fig=next_name(env["pd-work"])
                save_fig=os_join(f"{env['pd-work']}/{consensus}.pdf")
            )
        )


        for consensus in sorted(df_s["CONSENSUS_RBS_MAT"].unique()):
            df_c = df_s[df_s["CONSENSUS_RBS_MAT"] == consensus]
            pos_to_count = {x[0]: x[1] for x in df_c["Peak"].value_counts().items()}

            for i in range(15):
                if i not in pos_to_count:
                    pos_to_count[i] = 0

            sns.barplot(
                pd.DataFrame(
                    {"Number of Models": [pos_to_count[k] for k in sorted(pos_to_count.keys())],
                     "Peak of Spacer Distribution": sorted(pos_to_count.keys())}
                ),
                "Peak of Spacer Distribution", "Number of Models",
                figure_options=FigureOptions(
                    title=f"Shift: {s}, Example Consensus: {consensus}",
                    save_fig=next_name(env["pd-work"])
                )
            )

        # import matplotlib.pyplot as plt
        # plt.savefig(next_name(env["pd-work"]))

        # sns.distplot(df_s, "Peak", sns_kwargs={"kde": False}, figure_options=FigureOptions(
        #     title=f"Shift: {s}, Consensus: {consensus}"
        # ))







def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    df = load_gms2_models_from_pickle(args.pf_data)

    compare_msa_to_16s(env, df, args.seq_16s,
                       gc_range=args.gc_range,
                       group=args.group)


if __name__ == "__main__":
    main(my_env, parsed_args)
