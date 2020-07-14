# Author: Karl Gemayel
# Created: 7/3/20, 10:51 AM

import logging
import argparse
from typing import *
import matplotlib.pyplot as plt
import logomaker as lm

# noinspection All
import seaborn
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm

import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.mgm_model import MGMModel, MGMModelGC
from mg_general import Environment, add_env_args_to_parser
from mg_general.general import os_join, get_value
from mg_io.general import mkdir_p
from mg_models.gms2_noncoding import GMS2Noncoding
from mg_models.motif_model import MotifModel

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_viz import sns
from mg_viz.general import square_subplots

parser = argparse.ArgumentParser("Visualize an MGM model file.")

parser.add_argument('--pf-mgm-mod', required=True, help="MGM Model file")

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


def get_gms2_groups_for_genome_type(genome_type):
    # type: (str) -> List[str]
    return {
        "Archaea": ["A", "D"],
        "Bacteria": ["A", "B", "C", "X"],
    }[genome_type]


def visualize_motif(env, mgm, genome_type, tag, **kwargs):
    # type: (Environment, MGMModel, str, str, Dict[str, Any]) -> None

    pdf = get_value(kwargs, "pdf", None)  # type: PdfPages

    num_gc = len(mgm.items_by_species_and_gc[genome_type[0]])
    num_rows, num_cols = square_subplots(num_gc)

    all_gc_tags = sorted(mgm.items_by_species_and_gc[genome_type[0]])

    gc_mod_pair = dict()  # type: Dict[float, MGMModelGC]

    for i, gc_tag in enumerate(all_gc_tags):
        mgm_mod_gc = mgm.items_by_species_and_gc[genome_type[0]][gc_tag]

        if f"{tag}_MAT" not in mgm_mod_gc.items.keys():
            continue

        gc_mod_pair[gc_tag] = mgm_mod_gc

    if len(gc_mod_pair) == 0:
        return

    fig, axes = plt.subplots(num_rows, num_cols, sharex="all", sharey="all")
    fig_sp, axes_sp = plt.subplots(num_rows, num_cols, sharex="all", sharey="all")

    for i, gc_tag in enumerate(all_gc_tags):
        if gc_tag not in gc_mod_pair.keys():
            continue

        mgm_mod_gc = gc_mod_pair[gc_tag]
        motif_mat = MotifModel(mgm_mod_gc.items[f"{tag}_MAT"], mgm_mod_gc.items[f"{tag}_POS_DISTR"])
        nonc_mat = GMS2Noncoding(mgm_mod_gc.items[f"NON_MAT"])

        df_motif_mat = motif_mat.pwm_to_df()
        np_nonc_mat = nonc_mat.pwm_to_array(0)

        df_rel = lm.transform_matrix(
            df_motif_mat, from_type="probability", to_type="information", background=np_nonc_mat
        )

        ax = axes.ravel()[i]
        ax.set_title(gc_tag)
        lm.Logo(df_rel, ax=ax)

        ax.set_ylim(0, 2)

        # spacer
        ax = axes_sp.ravel()[i]
        x = [a for a in range(len(motif_mat._spacer))]
        y = motif_mat._spacer
        seaborn.lineplot(x, y, ax=ax)

    if pdf:
        pdf.savefig(fig)
        pdf.savefig(fig_sp)
    else:
        fig.show()
        fig_sp.show()


def visualize_start_context(env, mgm, genome_type, tag, **kwargs):
    # type: (Environment, MGMModel, str, str, Dict[str, Any]) -> None

    pdf = get_value(kwargs, "pdf", None)  # type: PdfPages

    all_gc_tags = sorted(mgm.items_by_species_and_gc[genome_type[0]])

    # get all possible words in start contexts
    all_words = sorted(set(
        w for gc_tag in all_gc_tags
        for w in mgm.items_by_species_and_gc[genome_type[0]][gc_tag].items[f"{tag}_MAT"]
        if f"{tag}_MAT" in mgm.items_by_species_and_gc[genome_type[0]][gc_tag].items

    ))

    all_positions = sorted(set(
        p for gc_tag in all_gc_tags
        for w in mgm.items_by_species_and_gc[genome_type[0]][gc_tag].items[f"{tag}_MAT"]
        if f"{tag}_MAT" in mgm.items_by_species_and_gc[genome_type[0]][gc_tag].items
        for p in range(len(mgm.items_by_species_and_gc[genome_type[0]][gc_tag].items[f"{tag}_MAT"][w]))

    ))

    if len(all_words) == 0:
        return

    num_rows, num_cols = square_subplots(len(all_words))

    for p in tqdm(all_positions, "Start Context", total=len(all_positions), leave=True, position=1):
        fig, axes = plt.subplots(num_rows, num_cols, sharex="all", sharey="all")

        for i, w in tqdm(enumerate(all_words), f"Words in position {p}", total=len(all_words), leave=True, position=2):
            ax = axes.ravel()[i]

            x = list()
            y = list()

            for gc_tag in all_gc_tags:
                mod = mgm.items_by_species_and_gc[genome_type[0]][gc_tag]
                if f"{tag}_MAT" in mod.items:
                    x.append(gc_tag)
                    y.append(mod.items[f"{tag}_MAT"][w][p])

            seaborn.lineplot(x, y, ax=ax)

        fig.suptitle(p)
        if pdf:
            pdf.savefig(fig)
        else:
            fig.show()


def visualize_start_codons(env, mgm, genome_type, gms2_group, **kwargs):
    # type: (Environment, MGMModel, str, str, Dict[str, Any]) -> None
    pdf = get_value(kwargs, "pdf", None)  # type: PdfPages

    starts = ["ATG", "GTG", "TTG"]
    all_gc_tags = sorted(mgm.items_by_species_and_gc[genome_type[0]])

    fig, axes = plt.subplots()

    for s in starts:
        x = list()
        y = list()

        for gc_tag in all_gc_tags:
            mod = mgm.items_by_species_and_gc[genome_type[0]][gc_tag]

            name = f"{s}_{gms2_group}"

            prob = mod.items.get(name)

            if prob is not None:
                x.append(float(gc_tag))
                y.append(float(prob))

        seaborn.lineplot(x, y, label=s)

    if pdf:
        pdf.savefig(fig)
    else:
        fig.show()


def viz_mgm_model_for_type_and_group(env, mgm, genome_type, gms2_group, **kwargs):
    # type: (Environment, MGMModel, str, str, Dict[str, Any]) -> None

    # open pdf
    with PdfPages(f"graphs_{genome_type}_{gms2_group}.pdf") as pdf:
        # visualize RBS
        visualize_motif(env, mgm, genome_type, f"RBS_{gms2_group}_0", pdf=pdf, **kwargs)

        # visualize Promoter
        visualize_motif(env, mgm, genome_type, f"PROMOTER_{gms2_group}_0", pdf=pdf, **kwargs)
        #
        # visualize RBS Start Context
        visualize_start_context(env, mgm, genome_type, f"SC_RBS_{gms2_group}", pdf=pdf, **kwargs)

        # visualize Start codons
        visualize_start_codons(env, mgm, genome_type, gms2_group, pdf=pdf, **kwargs)


def viz_mgm_model_for_type(env, mgm, genome_type, **kwargs):
    # type: (Environment, MGMModel, str, Dict[str, Any]) -> None

    valid_groups = get_gms2_groups_for_genome_type(genome_type)

    # visualize each group separately
    for g in tqdm(valid_groups, f"{genome_type}", total=len(valid_groups)):
        pd_work = os_join(env["pd-work"], g)
        mkdir_p(pd_work)
        curr_env = env.duplicate({"pd-work": pd_work})
        viz_mgm_model_for_type_and_group(curr_env, mgm, genome_type, g, **kwargs)


def viz_mgm_model(env, mgm, **kwargs):
    # type: (Environment, MGMModel, Dict[str, Any]) -> None

    viz_mgm_model_for_type(env, mgm, "Bacteria", **kwargs)
    viz_mgm_model_for_type(env, mgm, "Archaea", **kwargs)


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    mgm = MGMModel.init_from_file(args.pf_mgm_mod)
    viz_mgm_model(env, mgm)


if __name__ == "__main__":
    main(my_env, parsed_args)
