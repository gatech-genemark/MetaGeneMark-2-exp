# Author: Karl Gemayel
# Created: 6/27/20, 5:54 PM

import logging
import argparse
import os
from tempfile import mkstemp
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.genome_list import GenomeInfoList, GenomeInfo
from mg_container.mgm_model import MGMModel
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import get_value, os_join
from mg_io.general import remove_p, mkdir_p
from mg_models.shelf import run_mgm
from mg_options.parallelization import ParallelizationOptions
from mg_parallelization.generic_threading import run_n_per_thread
from mg_parallelization.pbs import PBS
from mg_pbs_data.mergers import merge_identity
from mg_pbs_data.splitters import split_gil

parser = argparse.ArgumentParser("Run MGM on list of genomes.")

parser.add_argument('--pf-gil', required=True, help="File containing list of genomes")
parser.add_argument('--pf-mgm-mod', required=True, help="MGM model file")
parser.add_argument('--dn-run', default="mgm")

parser.add_argument('--tag-value-pair', required=False, nargs="+", help="Set tags in MGM (e.g. RBS 0 PROMOTER 1)")

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


def parse_tags_from_list(list_tag_value_pairs):
    # type: (List[str]) -> List[Tuple]
    if list_tag_value_pairs is None:
        return list()

    if len(list_tag_value_pairs) % 2 != 0:
        raise ValueError("Tag/value pairs list must have a length multiple of 2.")

    list_parsed = list()
    for i in range(0, len(list_tag_value_pairs), 2):
        list_parsed.append((list_tag_value_pairs[i], list_tag_value_pairs[i + 1]))
    return list_parsed


def parse_and_set_tags(list_tag_value_pairs, mgm):
    # type: (List[str], MGMModel) -> None
    list_tag_value_pairs = parse_tags_from_list(list_tag_value_pairs)

    for tv in list_tag_value_pairs:
        for species in mgm.items_by_species_and_gc.keys():
            for gc in mgm.items_by_species_and_gc[species].keys():
                mgm.items_by_species_and_gc[species][gc].items[tv[0]] = tv[1]


def run_mgm_on_gi(env, gi, pf_mgm_mod, **kwargs):
    # type: (Environment, GenomeInfo, str, Dict[str, Any]) -> None
    pf_sequence = os_join(env["pd-data"], gi.name, "sequence.fasta")
    dn_run = get_value(kwargs, "dn_run", "mgm", valid_type=str)

    pd_gi_run = os_join(env["pd-runs"], gi.name, dn_run)
    mkdir_p(pd_gi_run)
    pf_prediction = os_join(pd_gi_run, "prediction.gff")

    run_mgm(env, pf_sequence, pf_mgm_mod, pf_prediction)


def helper_run_mgm_on_genome_list(env, gil, pf_mgm_mod, **kwargs):
    # type: (Environment, GenomeInfoList, str, Dict[str, Any]) -> None
    for gi in gil:
        run_mgm_on_gi(env, gi, pf_mgm_mod, **kwargs)


def run_mgm_on_genome_list(env, gil, pf_mgm_mod, **kwargs):
    # type: (Environment, GenomeInfoList, str, Dict[str, Any]) -> None
    """
    Run MGM on list of genomes
    :param env: Environment variable
    :param gil: List of genomes
    :param pf_mgm_mod: Path to MGM Model file
    :param kwargs: Options for running MGM
        prl-options: ParallelizationOptions object that determines how to parallelize runs
        dn-run: Name of run directory (e.g. $PATH_TO_RUN/GENOME/{dn-run}
    :return: None
    """

    prl_options = get_value(kwargs, "prl_options", None)  # type: ParallelizationOptions

    # No parallelization
    if prl_options is None:
        helper_run_mgm_on_genome_list(env, gil, pf_mgm_mod, **kwargs)
    else:
        # PBS Parallelization
        if prl_options["use-pbs"]:
            # setup PBS jobs
            pbs = PBS(env, prl_options, splitter=split_gil, merger=merge_identity)
            pbs.run(
                data={"gil": gil},
                func=helper_run_mgm_on_genome_list,
                func_kwargs={"env": env, "pf_mgm_mod": pf_mgm_mod, **kwargs}
            )
        # Multithreading parallelization
        else:
            # parallel using threads
            run_n_per_thread(
                list(gil), run_mgm_on_gi, "gi",
                {"env": env, "pf_mgm_mod": pf_mgm_mod, **kwargs},
                simultaneous_runs=prl_options.safe_get("num-processors")
            )


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    gil = GenomeInfoList.init_from_file(os.path.abspath(args.pf_gil))
    mgm_mod = MGMModel.init_from_file(args.pf_mgm_mod)
    parse_and_set_tags(args.tag_value_pair, mgm_mod)

    # write to temporary file
    _, pf_mgm_mod = mkstemp(prefix=env["pd-work"], suffix="_mgm.mod")
    mgm_mod.to_file(pf_mgm_mod)
    run_mgm_on_genome_list(env, gil, pf_mgm_mod, dn_run=args.dn_run)

    # delete temporary mgm file
    remove_p(pf_mgm_mod)


if __name__ == "__main__":
    main(my_env, parsed_args)
