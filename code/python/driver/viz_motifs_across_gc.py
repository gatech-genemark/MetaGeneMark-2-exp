# Author: Karl Gemayel
# Created: 7/8/20, 8:56 AM

import umap
import umap.plot
import logging
import numpy as np
import argparse
from typing import *
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import next_name
from mg_models.shelf import read_archaea_bacteria_inputs
from mg_viz.general import save_figure, FigureOptions
from mg_viz.shelf import create_mappable_for_colorbar

parser = argparse.ArgumentParser("DRIVER DESCRIPTION.")

parser.add_argument('--pf-bac', required=True, help="Collected GMS2 model files for bacteria")
parser.add_argument('--pf-arc', required=True, help="Collected GMS2 model files for archaea")

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


def create_numpy_for_column(df, col):
    # type: (pd.DataFrame, str) -> np.ndarray

    example = df.at[df.index[0], "Mod"].items[col]

    n = len(df)                         # number of examples
    w = len(next(iter(example.values())))        # width (numbere of positions)
    b = len(example)                    # number of bases (letters)

    mat = np.empty((n,w,b), dtype=float)

    # fill the array
    for n_pos, idx in enumerate(df.index):
        dict_arr = df.at[idx, "Mod"].items[col]


        # for each base
        for b_pos, letter in enumerate(sorted(dict_arr.keys())):
            for w_pos, value in enumerate(dict_arr[letter]):
                if w_pos == 6:
                    print("HI")
                mat[n_pos, w_pos, b_pos] = value

    return mat


def visualize_matrix_column(env, df, col):
    # type: (Environment, pd.DataFrame, str) -> None

    fp = FontProperties()
    fp.set_family("monospace")

    # create N x 6 x 4 matrix for RBS
    mat = create_numpy_for_column(df, col)
    mat = mat.reshape((mat.shape[0], mat.shape[1] * mat.shape[2]))

    # get interesting features to view data by
    gc = df["GC"]
    group = df["GENOME_TYPE"]

    for r in range(1):

        reducer = umap.UMAP(random_state=r)
        reducer = reducer.fit(mat)
        embedding = reducer.embedding_
        print(embedding.shape)

        # fig, ax = plt.subplots()
        #
        # plt.scatter(embedding[:, 0], embedding[:, 1], c=gc, marker="+")
        # plt.colorbar()
        # plt.show()
        # themes = ["fire", "viridis", "inferno", "blue", "red", "green", "darkblue", "darkred", "darkgreen"]
        # fig, axes = plt.subplots(3, 3)
        # for ax, theme in zip(axes.ravel(), themes):
        #     fig, ax = plt.subplots()
        #     umap.plot.points(reducer, values=gc, theme=theme, )
        #     plt.show()
        ax = umap.plot.points(reducer, values=gc, cmap="viridis")
        mappable = create_mappable_for_colorbar(gc, "viridis")
        plt.colorbar(mappable)
        plt.title(col.replace("_", " "))
        plt.tight_layout()
        save_figure(FigureOptions(
            save_fig=next_name(env["pd-work"])
        ))
        plt.show()

        umap.plot.points(reducer, labels=group.values, color_key_cmap="Paired")
        plt.title(col.replace("_", " "))
        plt.tight_layout()
        save_figure(FigureOptions(
            save_fig=next_name(env["pd-work"])
        ))
        plt.show()

        # umap.plot.points(reducer, labels=group.values, color_key_cmap="Dark2")
        # plt.title(col)
        # save_figure(FigureOptions(
        #     save_fig=next_name(env["pd-work"])
        # ))
        # plt.show()


        umap.plot.points(reducer, labels=df["Type"])
        plt.title(col.replace("_", " "))
        plt.tight_layout()
        save_figure(FigureOptions(
            save_fig=next_name(env["pd-work"])
        ))
        plt.show()

    # fig, ax = plt.subplots()
    # plt.scatter(embedding[:, 0], embedding[:, 1], c=group.values, cmap=cm.brg)
    # plt.show()


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df = read_archaea_bacteria_inputs(args.pf_arc, args.pf_bac)
    logger.debug(f"Removing genetic code 4: {(df['Genetic Code'] == 4).sum()}")
    df = df[df["Genetic Code"] != 4]
    df = df.convert_dtypes().copy()

    visualize_matrix_column(env, df, "RBS_MAT")
    visualize_matrix_column(env, df[(df["Type"] == "Bacteria") & (df["GENOME_TYPE"] == "C")], "PROMOTER_MAT")

    visualize_matrix_column(env, df[(df["Type"] == "Archaea") & (df["GENOME_TYPE"] == "D")], "PROMOTER_MAT")
    # visualize_matrix_column(env, df[(df["GENOME_TYPE"].isin({"C", "D"}))], "PROMOTER_MAT")


if __name__ == "__main__":
    main(my_env, parsed_args)
