# Author: Karl Gemayel
# Created: 8/26/20, 8:23 PM

import logging
import argparse
import pandas as pd
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


parser = argparse.ArgumentParser("DRIVER DESCRIPTION.")

parser.add_argument('--pf-data', required=True)

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

    df = pd.read_csv(args.pf_data)

    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    import seaborn as sns

    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.plot_trisurf(df['p4'], df['p11'], df['Match Rate'], linewidth=0.2)
    # plt.show()
    # df_sample = df[df["Tool"] == "mgm2"].sample(1000)
    # df = pd.concat([df_sample, df[df["Tool"] == "mprodigal"]])



    import numpy as np
    df["Params"] = np.nan
    all_vals = df.groupby(["p4", "p11"], as_index=False).sum()      # type: pd.DataFrame
    all_vals = all_vals.sort_values(["p4", "p11"])

    for counter, i in enumerate(all_vals.index):
        p4 = all_vals.at[i, "p4"]
        p11 = all_vals.at[i, "p11"]
        df.loc[(df["Tool"] == "mgm2") & (df["p4"] == p4) & (df["p11"] == p11), "Params"] = counter

    df.loc[df["Tool"] == "mprodigal", "Params"] = -1

    df_avg = df.groupby(["Params", "Chunk Size", "Tool"], as_index=False).mean()
    df_avg = df_avg.sort_values("Chunk Size")

    fig, ax = plt.subplots()
    sns.lineplot("Chunk Size", "Match Rate", data=df_avg[df_avg["Tool"] == "mgm2"], hue="Params", ax=ax, legend=False)
    df_prod = df_avg[df_avg["Tool"] == "mprodigal"]
    ax.plot(df_prod["Chunk Size"], df_prod["Match Rate"], label="MProdigal")

    plt.show()

    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # ax.scatter(df['Chunk Size'], df['p4,p11'], df['Match Rate'], s=0.2)
    # ax.set_xlabel("Chunk Size")
    # ax.set_ylabel("Parameters")
    # ax.set_zlabel("Match Rate")
    # plt.show()
    return

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(df['Chunk Size'], df['p4,p11'], df['Match Rate'], linewidth=0.2)
    ax.set_xlabel("Chunk Size")
    ax.set_ylabel("Parameters")
    ax.set_zlabel("Match Rate")
    plt.show()

    df2 = df[df["Tool"] == "mgm2"].groupby(["p4", "p11"], as_index=False).mean()

    idx = df2["Match Rate"].argmax()
    p4 = df2.at[idx, "p4"]
    p11 = df2.at[idx, "p11"]
    df_best = df[(df["p4"] == p4) & (df["p11"] == p11)]
    df_alex = df[(df["p4"] == 10) & (df["p11"] == 20)]
    fig, ax = plt.subplots()
    sns.lineplot("Chunk Size", "Match Rate", data=df_best, label="Optimized")
    sns.lineplot("Chunk Size", "Match Rate", data=df[df["Tool"] == "mprodigal"], label="MProdigal")
    sns.lineplot("Chunk Size", "Match Rate", data=df_alex, label="Original")
    ax.set_ylim(0, 1)
    plt.show()




if __name__ == "__main__":
    main(my_env, parsed_args)
