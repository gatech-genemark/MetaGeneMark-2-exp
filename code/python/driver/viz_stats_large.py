# Author: Karl Gemayel
# Created: 8/5/20, 7:59 AM

import logging
import argparse
import pandas as pd
from typing import *

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import mg_log  # runs init in mg_log and configures logger

# Custom imports
from mg_argparse.stats import add_stats_options
from mg_general import Environment, add_env_args_to_parser
from mg_general.general import os_join, get_value
from mg_io.general import mkdir_p
from mg_stats.shelf import check_tools_and_reference_lists, read_small_stats_per_gene

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_viz.stats_large import viz_stats_large_3p, viz_stats_large_5p

parser = argparse.ArgumentParser("Visualize stats collected for large set of genomes"
                                 " (e.g. > 20 genomes).")

add_stats_options(parser)

parser.add_argument('--pf-checkpoint-3p')
parser.add_argument('--pf-checkpoint-5p')

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


def viz_stats_large(env, df_per_gene, tools, ref_5p, ref_3p, **kwargs):
    # type: (Environment, pd.DataFrame, List[str], List[str], List[str], Dict[str, Any]) -> None
    pd_3p = os_join(env["pd-work"], "figures_3p")
    pd_5p = os_join(env["pd-work"], "figures_5p")
    mkdir_p(pd_3p)
    mkdir_p(pd_5p)

    viz_stats_large_3p(env.duplicate({"pd-work": pd_3p}), df_per_gene, tools, ref_3p,
                       pf_checkpoint=kwargs.get("pf_checkpoint_3p"))
    viz_stats_large_5p(env.duplicate({"pd-work": pd_5p}), df_per_gene, tools, ref_5p,
                       pf_checkpoint=kwargs.get("pf_checkpoint_5p"))


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None
    df_per_gene = read_small_stats_per_gene(args.pf_data, args.parse_names)
    # df_per_gene = df_per_gene.head(100000).copy()
    logger.debug(f"Num genomes = {len(list(df_per_gene['Genome'].unique()))}")
    tools = check_tools_and_reference_lists(df_per_gene, args.tools, args.ref_5p, args.ref_3p)

    viz_stats_large(env, df_per_gene, tools, args.ref_5p, args.ref_3p,
                    pf_checkpoint_3p=args.pf_checkpoint_3p,
                    pf_checkpoint_5p=args.pf_checkpoint_5p)


if __name__ == "__main__":
    main(my_env, parsed_args)
