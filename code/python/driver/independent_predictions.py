# Author: Karl Gemayel
# Created: 6/24/20, 9:08 AM

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

# FILL IN ARGUMENTS

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


def ratio_false_true(p_true):
    # type: (float) -> float
    if p_true == 0:
        return float('inf')
    return 1.0 / p_true - 1


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    pass


if __name__ == "__main__":
    main(my_env, parsed_args)
