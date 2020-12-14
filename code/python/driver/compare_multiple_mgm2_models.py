# Author: Karl Gemayel
# Created: 6/22/20, 4:43 PM

import os
import logging
import argparse
import pandas as pd
from typing import *
import matplotlib.pyplot as plt

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
from mg_general.general import os_join, fix_names, next_name
from mg_io.general import mkdir_p
from mg_models.shelf import run_gms2_with_component_toggles_and_get_accuracy, component_in_model_file, \
    run_mgm_and_get_accuracy, run_mgm2_and_get_accuracy
from mg_viz import sns
from mg_viz.general import FigureOptions
from mg_viz.shelf import get_order_by_rank

parser = argparse.ArgumentParser("Test MGM2 models on verified set.")

parser.add_argument('--pf-gil', required=True, help="Genome information list")

group_mgm2 = parser.add_mutually_exclusive_group(required=True)
group_mgm2.add_argument('--pf-mgm2', nargs="+")
group_mgm2.add_argument('--pd-mgm2', type=str, help="Directory containing all mgm2 models (ending in .mod)")

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


def test_component_for_gi(env, gi, list_pf_mgm):
    # type: (Environment, GenomeInfo, List[str]) -> pd.DataFrame

    list_entries = list()

    pd_gi = os_join(env["pd-work"], gi.name)
    mkdir_p(pd_gi)

    env_dup = env.duplicate({"pd-work": pd_gi})

    for pf_mgm in list_pf_mgm:
        basename = os.path.basename(pf_mgm)
        logger.info(f"Model file: {pf_mgm}")
        results = run_mgm2_and_get_accuracy(env_dup, gi, pf_mgm)
        list_entries.append({"Tool": basename, **results})

    return pd.DataFrame(list_entries)


def test_mgm2_models_on_verified(env, gil, list_pf_mgm2):
    # type: (Environment, GenomeInfoList, List[str]) -> None

    list_df = list()

    for gi in gil:
        logger.info(f"Genome: {gi.name}")
        list_df.append(test_component_for_gi(env, gi, list_pf_mgm2))

        list_df[-1]["Genome"] = gi.name

    df = pd.concat(list_df, ignore_index=True, sort=False)

    df["Genome"] = df.apply(fix_names, axis=1)
    print(df.to_csv())

    hue_order = get_order_by_rank(df, "Genome", "Error", "Tool")

    fig, ax = plt.subplots(figsize=(12, 4))
    sns.barplot(df, "Genome", "Error", hue="Tool",
                ax=ax,
                sns_kwargs={"hue_order": hue_order},
                figure_options=FigureOptions(
                    save_fig=next_name(env["pd-work"])
                ),
                )

    df_piv = df.pivot(index="Genome", columns="Tool", values=["Error", "Number of Errors"])
    print(df_piv.to_csv())


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_gil)

    if args.pf_mgm2:
        pf_mgm2 = args.pf_mgm2
    else:
        pd_mgm2 = os.path.abspath(args.pd_mgm2)
        pf_mgm2 = list()
        for file in os.listdir(pd_mgm2):
            if file.endswith(".mod"):
                pf_mgm2.append(os.path.join(pd_mgm2, file))

    test_mgm2_models_on_verified(env, gil, pf_mgm2)


if __name__ == "__main__":
    main(my_env, parsed_args)
