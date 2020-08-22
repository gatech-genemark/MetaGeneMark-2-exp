# Author: Karl Gemayel
# Created: 8/22/20, 3:08 PM

import logging
import argparse
import os

import seaborn
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO

from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.genome_list import GenomeInfoList, GenomeInfo
from mg_container.gms2_mod import GMS2Mod
from mg_general.general import os_join, get_value
from mg_general.shelf import compute_gc
from mg_general import Environment, add_env_args_to_parser
from mg_io.general import save_obj, load_obj

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

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
        "Group": [mod.items["GENOME_TYPE"] for mod in list_mod]
    })
    seaborn.catplot("Group", "GC", data=df)
    plt.show()


def viz_gms2_models_over_gc(env, list_gi, list_mod, list_gc, **kwargs):
    # type: (Environment, List[GenomeInfo], List[GMS2Mod], List[float], Dict[str, Any]) -> None
    viz_genome_type_per_gc(env, list_gi, list_mod, list_gc)


def read_genome_data(env, gil, **kwargs):
    dn_gms2 = get_value(kwargs, "dn_gms2", "gms2")

    list_gi = list()
    list_mod = list()
    list_gc = list()

    for gi in gil:
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
