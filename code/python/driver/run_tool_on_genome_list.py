# Author: Karl Gemayel
# Created: 6/27/20, 5:54 PM

import os
import logging
import argparse
from subprocess import CalledProcessError
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_argparse.parallelization import add_parallelization_options
from mg_io.general import mkdir_p
from mg_container.genome_list import GenomeInfo, GenomeInfoList
from mg_general import add_env_args_to_parser, Environment
from mg_general.general import get_value, os_join, run_shell_cmd
from mg_models.shelf import run_gms2, run_mgm2, run_mgm, run_prodigal, run_meta_prodigal
from mg_options.parallelization import ParallelizationOptions
from mg_parallelization.generic_threading import run_n_per_thread
from mg_parallelization.pbs import PBS
from mg_pbs_data.mergers import merge_identity
from mg_pbs_data.splitters import split_gil

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Run external prediction tools on genome list.")

parser.add_argument('--pf-gil', required=True, help="List of genomes")
parser.add_argument('--type', required=True, choices=["archaea", "bacteria", "auto"], help="Is the list archaea or bacteria")
parser.add_argument('--dn-run', required=False, help="Name of directory that will contain the run")
parser.add_argument('--skip-if-exists', action="store_true", default=False, help="If set, tool isn't run if predictions.gff file exists")
parser.add_argument('--tool', choices=["gms2", "prodigal", "mgm", "mgm2", "mprodigal"], required=True, help="Tool used for prediction")

parser.add_argument('--pf-mgm-mod', help="Path to MGM model file", type=os.path.abspath)
parser.add_argument('--pf-mgm2-mod', help="Path to MGM model file", type=os.path.abspath)

add_parallelization_options(parser)

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


def run_tool_on_gi(env, gi, tool, **kwargs):
    # type: (Environment, GenomeInfo, str, Dict[str, Any]) -> None
    dn_run = get_value(kwargs, "dn_run", tool, default_if_none=True)
    skip_if_exists = get_value(kwargs, "skip_if_exists", False)

    pd_work = os_join(env["pd-runs"], gi.name, dn_run)
    curr_env = env.duplicate({"pd-work": pd_work})
    mkdir_p(pd_work)

    pf_sequence = os_join(env["pd-data"], gi.name, "sequence.fasta")
    pf_prediction = os_join(pd_work, "prediction.gff")

    if skip_if_exists and os.path.isfile(pf_prediction):
        return

    try:
        if tool == "gms2":
            run_gms2(curr_env, pf_sequence, pf_prediction, **kwargs)
        elif tool == "prodigal":
            run_prodigal(curr_env, pf_sequence, pf_prediction, **kwargs)
        elif tool == "mprodigal":
            run_meta_prodigal(curr_env, pf_sequence, pf_prediction, **kwargs)
        elif tool == "mgm2":
            run_mgm2(curr_env, pf_sequence, kwargs.get("pf_mgm2_mod"), pf_prediction)
        elif tool == "mgm":
            run_mgm(curr_env, pf_sequence, kwargs.get("pf_mgm_mod"), pf_prediction)
        else:
            raise NotImplementedError()
    except CalledProcessError:
        logger.warning(f"Could not run {tool} on {gi.name}")

# def run_tool_on_gil(env, gil, tool, **kwargs):
#     # type: (Environment, GenomeInfoList, str, Dict[str, Any]) -> None
#
#     logger.info("Running tool {} on {} genomes".format(tool, len(gil)))
#     func = {
#         "gms2": run_tool_on_gi,
#         "prodigal": run_tool_on_gi,
#     }[tool]
#
#     for gi in gil:
#         func(env, gi, tool, **kwargs)


def helper_run_tool_on_genome_list(env, gil, tool, **kwargs):
    # type: (Environment, GenomeInfoList, str, Dict[str, Any]) -> None
    for gi in gil:
        run_tool_on_gi(env, gi, tool, **kwargs)


def run_tool_on_genome_list(env, gil, tool, **kwargs):
    # type: (Environment, GenomeInfoList, str, Dict[str, Any]) -> None
    """
    Run MGM on list of genomes
    :param env: Environment variable
    :param gil: List of genomes
    :param kwargs: Options for running MGM
        prl-options: ParallelizationOptions object that determines how to parallelize runs
        dn-run: Name of run directory (e.g. $PATH_TO_RUN/GENOME/{dn-run}
    :return: None
    """

    prl_options = get_value(kwargs, "prl_options", None)  # type: ParallelizationOptions

    # No parallelization
    if prl_options is None:
        helper_run_tool_on_genome_list(env, gil, tool, **kwargs)
    else:
        # PBS Parallelization
        if prl_options["use-pbs"]:
            # setup PBS jobs
            pbs = PBS(env, prl_options, splitter=split_gil, merger=merge_identity)
            pbs.run(
                data=gil,
                func=helper_run_tool_on_genome_list,
                func_kwargs={"env": env, "tool": tool, **kwargs}
            )
        # Multithreading parallelization
        else:
            # parallel using threads
            run_n_per_thread(
                list(gil), run_tool_on_gi, "gi",
                {"env": env, "tool": tool, **kwargs},
                simultaneous_runs=prl_options.safe_get("num-processors")
            )


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    gil = GenomeInfoList.init_from_file(args.pf_gil)
    prl_options = ParallelizationOptions.init_from_dict(env, args.pf_parallelization_options, vars(args))

    # write to temporary file

    run_tool_on_genome_list(env, gil, args.tool, prl_options=prl_options,
                            dn_run=args.dn_run, genome_type=args.type,
                            skip_if_exists=args.skip_if_exists,
                            pf_mgm_mod=args.pf_mgm_mod, pf_mgm2_mod=args.pf_mgm2_mod)


if __name__ == "__main__":
    main(my_env, parsed_args)
