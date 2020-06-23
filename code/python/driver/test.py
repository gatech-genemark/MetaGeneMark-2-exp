# Author: karl
# Created: 2020-06-21, 8:29 a.m.

import logging
import argparse
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.mgm_model import MGMModel
from mg_general import Environment, add_env_args_to_parser
from mg_general.shelf import test_log_level

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_options.mg_models import MGModelsOptions

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


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    test_log_level()
    logger.debug("New")

    model = MGMModel.init_from_file("mgm_11.mod")

    model.to_file("output.mod")


if __name__ == "__main__":
    main(my_env, parsed_args)
