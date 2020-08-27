# Author: Karl Gemayel
# Created: 8/26/20, 8:53 PM

import logging
import argparse
import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_container.genome_list import GenomeInfoList, GenomeInfo
from mg_general import Environment, add_env_args_to_parser

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #


parser = argparse.ArgumentParser("DRIVER DESCRIPTION.")

parser.add_argument('--pf-gil', required=True)

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


def build_features_from_gi(env, gi):
    # type: (Environment, GenomeInfo) -> pd.Series
    dict_entries = dict()

    pf_mgm2 = pd.DataFrame()

    return pd.Series(dict_entries)

def build_features_from_gil(env, gil):
    # type: (Environment, GenomeInfoList) -> pd.DataFrame

    list_entries = list()
    for gi in gil:
        entry = build_features_from_gi(env, gi)
        list_entries.append(entry)

    return pd.DataFrame(list_entries)



def build_mgm2_dataset(env, gil):
    # type: (Environment, GenomeInfoList) -> None


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    gil = GenomeInfoList.init_from_file(args.pf_gil)

    df = build_mgm2_dataset(env, gil)


if __name__ == "__main__":
    main(my_env, parsed_args)
