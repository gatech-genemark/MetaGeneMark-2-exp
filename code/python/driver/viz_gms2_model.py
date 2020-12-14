# Author: Karl Gemayel
# Created: 8/22/20, 1:51 PM

import logging
import argparse
import seaborn
from typing import *
import logomaker as lm
import matplotlib.pyplot as plt

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.gms2_mod import GMS2Mod
from mg_general.general import next_name
from mg_models.motif_model import MotifModel
from mg_models.gms2_noncoding import GMS2Noncoding
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_viz.general import set_size

parser = argparse.ArgumentParser("DRIVER DESCRIPTION.")

parser.add_argument('--pf-mod', required=True)
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


def viz_motif_and_spacer(env, mod, motif, **kwargs):
    # type: (Environment, GMS2Mod, str, Dict[str, Any]) -> None

    motif_mat = MotifModel(mod.items[f"{motif}_MAT"], mod.items[f"{motif}_POS_DISTR"])
    nonc_mat = GMS2Noncoding(mod.items[f"NON_MAT"])

    df_motif_mat = motif_mat.pwm_to_df()
    np_nonc_mat = nonc_mat.pwm_to_array(0)

    df_rel = lm.transform_matrix(
        df_motif_mat, from_type="probability", to_type="information", background=np_nonc_mat
    )

    figsize = set_size("thesis", subplots=(1, 2), titles=True)
    fig, (ax_motif, ax_spacer) = plt.subplots(1, 2, figsize=figsize)

    # Plot motif
    ax = ax_motif
    ax.set_title(motif)
    lm.Logo(df_rel, ax=ax)
    ax.set_ylim(0, 2)

    # Plot spacer
    ax = ax_spacer
    x = [a for a in range(len(motif_mat._spacer))]
    y = motif_mat._spacer
    seaborn.lineplot(x, y, ax=ax)

    fig.tight_layout()
    fig.savefig(next_name(env["pd-work"]))
    fig.show()


def viz_rbs_and_spacer(env, mod, **kwargs):
    # type: (Environment, GMS2Mod, Dict[str, Any]) -> None
    if "RBS" not in mod.items:
        return
    viz_motif_and_spacer(env, mod, "RBS", **kwargs)


def viz_promoter_and_spacer(env, mod, **kwargs):
    # type: (Environment, GMS2Mod, Dict[str, Any]) -> None
    if "PROMOTER" not in mod.items:
        return
    viz_motif_and_spacer(env, mod, "PROMOTER", **kwargs)


def viz_gms2_model(env, mod, **kwargs):
    # type: (Environment, GMS2Mod, Dict[str, Any]) -> None

    viz_rbs_and_spacer(env, mod, **kwargs)
    viz_promoter_and_spacer(env, mod, **kwargs)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    mod = GMS2Mod.init_from_file(args.pf_mod)
    viz_gms2_model(env, mod)


if __name__ == "__main__":
    main(my_env, parsed_args)
