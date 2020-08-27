# Author: Karl Gemayel
# Created: 8/21/20, 11:45 AM

import logging
import argparse
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.gms2_mod import GMS2Mod
from mg_container.mgm_model import MGMModel
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #


parser = argparse.ArgumentParser("Update target model file with components from source model file.")

parser.add_argument('--pf-source', required=True)
parser.add_argument('--pf-target', required=True)
parser.add_argument('--pf-output', required=True)
parser.add_argument('--groups', required=True, nargs="+", choices=["A", "B", "C", "D", "X"], type=str.upper)
parser.add_argument('--components', required=True, nargs="+",
                    choices=["RBS", "PROMOTER", "SC_RBS", "SC_PROMOTER", "STARTS"],
                    type=str.upper)

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

def get_tags_for_motif(motif, group):
    # type: (str, str) -> List[str]
    motif = motif.upper()
    group = group.upper()
    prefix = f"{motif}_{group}"

    entries = list()
    for index in [0,1,2]:

        p_index = f"{prefix}_{index}"

        entries += [
            f"{p_index}",
            f"{p_index}_MARGIN",
            f"{p_index}_MAT",
            f"{p_index}_MAX_DUR",
            f"{p_index}_ORDER",
            f"{p_index}_POS_DISTR",
            f"{p_index}_SHIFT",
            f"{p_index}_WIDTH",
        ]
    return entries

def get_tags_for_start_context(sc_tag, group):
    # type: (str, str) -> List[str]
    sc_tag = sc_tag.upper()
    group = group.upper()

    return [
        f"{sc_tag}_{group}",
        f"{sc_tag}_{group}_MARGIN",
        f"{sc_tag}_{group}_MAT",
        f"{sc_tag}_{group}_ORDER",
        f"{sc_tag}_{group}_WIDTH",
    ]

def component_to_tags(component, group):
    # type: (str, str) -> List[str]
    if component.upper() == "RBS" or component.upper() == "PROMOTER":
        return get_tags_for_motif(component, group)
    if component.upper() == "STARTS":
        return [f"{s}_{group}" for s in ["ATG", "GTG", "TTG"]]
    if component.upper() in {"SC_RBS", "SC_PROMOTER"}:
        return get_tags_for_start_context(component.upper(), group)
    return list()


def get_tags_from_components(list_components, list_groups):
    # type: (List[str], List[str]) -> List[str]

    list_tags = list()
    for c in list_components:
        for g in list_groups:
            list_tags += component_to_tags(c, g)

    return list_tags


def update_target_model_from_source(mod_source, mod_target, list_tags):
    # type: (MGMModel, MGMModel, List[str]) -> MGMModel
    import copy
    mod_output = copy.deepcopy(mod_target)

    for species in mod_source.items_by_species_and_gc.keys():
        for gc in mod_source.items_by_species_and_gc[species].keys():
            target = mod_output.items_by_species_and_gc[species][gc]
            source = mod_source.items_by_species_and_gc[species][gc]

            for tag in list_tags:
                if tag in source.items.keys():
                    target.items[tag] = source.items[tag]

    return mod_output


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    mod_source = MGMModel.init_from_file(args.pf_source)
    mod_target = MGMModel.init_from_file(args.pf_target)

    list_groups = args.groups
    list_tags = get_tags_from_components(args.components, list_groups)

    mod_output = update_target_model_from_source(mod_source, mod_target, list_tags)
    mod_output.to_file(args.pf_output)




if __name__ == "__main__":
    main(my_env, parsed_args)
