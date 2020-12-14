# Author: Karl Gemayel
# Created: 7/10/20, 7:40 AM

import logging
import argparse
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

parser.add_argument('--pf-gil', required=True, help="List of genomes")
parser.add_argument('--dn-tools', nargs="+", required=False, help="Name of directories that will contain the run. "
                                                                "If set, this should be a parallel list of names "
                                                                "to the '--tools' argument. If not set, the "
                                                                "values of '--tools' are used as directory names.")
parser.add_argument('--tools', nargs="+",
                    choices=["gms2", "prodigal", "mgm", "mgm2", "mprodigal"], required=True,
                    help="Tool used for prediction")

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
    tools = args.tools
    dn_tools = args.dn_tools if args.dn_tools is not None else tools

    # check that both have the same length
    if len(tools) != len(dn_tools):
        raise ValueError(f"The 'tools' and 'dn-tools' arguments"
                         f" must have equal lengths: {len(tools)} != {len(dn_tools)}")


if __name__ == "__main__":
    main(my_env, parsed_args)
