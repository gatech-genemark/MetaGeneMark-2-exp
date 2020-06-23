import logging
import argparse
import matplotlib.pyplot as plt
import pandas as pd

# noinspection All
import pathmagic

# noinspection PyUnresolvedReferences
import sbsp_log  # runs init in sbsp_log and configures logger

# Custom imports
from mg_container.genome_list import GenomeInfoList, GenomeInfo
from mg_general import Environment

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #
from mg_general.general import os_join, next_name, fix_names
from mg_io.general import mkdir_p
from mg_models.shelf import run_gms2_with_component_toggles_and_get_accuracy, \
    component_in_model_file
from mg_viz import sns
from mg_viz.general import FigureOptions
from mg_viz.colormap import ColorMap as CM

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-genome-list', required=True, help="Genome information list")

parser.add_argument('--pd-work', required=False, default=None, help="Path to working directory")
parser.add_argument('--pd-data', required=False, default=None, help="Path to data directory")
parser.add_argument('--pd-results', required=False, default=None, help="Path to results directory")
parser.add_argument("-l", "--log", dest="loglevel", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    help="Set the logging level", default='WARNING')

parsed_args = parser.parse_args()

# ------------------------------ #
#           Main Code            #
# ------------------------------ #

# Load environment variables
my_env = Environment(pd_data=parsed_args.pd_data,
                     pd_work=parsed_args.pd_work,
                     pd_results=parsed_args.pd_results)

# Setup logger
logging.basicConfig(level=parsed_args.loglevel)
logger = logging.getLogger("logger")  # type: logging.Logger


def analyze_gms2_components_on_verified_set_for_gi(env, gi):
    # type: (Environment, GenomeInfo) -> pd.DataFrame

    list_entries = list()

    start_components = {
        "Start Codons", "Start Context", "RBS", "Promoter",
    }

    pd_gi = os_join(env["pd-work"], gi.name)
    mkdir_p(pd_gi)

    # for each component to keep on
    for component_on in sorted(start_components) + ["MGM2*", "MGM", "GMS2"]:
        components_off = start_components.difference({component_on})

        if component_on == "MGM2*" or component_on == "GMS2":
            components_off = set()
        elif component_on == "MGM":
            pass
        elif not component_in_model_file(env, gi, component_on) and component_on not in {"MGM2*", "MGM", "GMS2"}:
            continue

        native_coding_off = False if component_on == "GMS2" else True

        pd_gi_component = os_join(pd_gi, component_on).replace(" ", "")
        mkdir_p(pd_gi_component)

        env_dup = env.duplicate({"pd-work": pd_gi_component})

        if component_on == "Start Context":
            component_on = {component_on}  # "rbs", "promoter"}
            components_off.remove("RBS")
            components_off.remove("Promoter")
        else:
            component_on = {component_on}


        results = run_gms2_with_component_toggles_and_get_accuracy(env_dup, gi, components_off,
                                                                   native_coding_off=native_coding_off)

        list_entries.append({
            "Genome": gi.name,
            "Component": next(iter(component_on)).replace("_", "-"),
            # **{con: True for con in component_on},                             # current component is on
            # **{coff: False for coff in components_off},     # all others are off
            **results
        })



    return pd.DataFrame(list_entries)


def analyze_gms2_components_on_verified_set(env, gil):
    # type: (Environment, GenomeInfoList) -> None

    # run different components
    list_df = list()
    for gi in gil:
        list_df.append(
            analyze_gms2_components_on_verified_set_for_gi(env, gi)
        )

    df = pd.concat(list_df, ignore_index=True, sort=False)
    df["Genome"] = df.apply(fix_names, axis=1)
    print(df.to_csv())


    fig, ax = plt.subplots(figsize=(12,4))
    sns.barplot(df, "Genome", "Error", hue="Component",
                ax=ax,
                figure_options=FigureOptions(
                    save_fig=next_name(env["pd-work"])
                ),
                sns_kwargs={
                    "hue_order": reversed(["GMS2", "MGM2*", "Start Context", "RBS", "Start Codons", "Promoter", "MGM"]),
                    "palette": CM.get_map("gms2_components")

                })




def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    gil = GenomeInfoList.init_from_file(args.pf_genome_list)

    analyze_gms2_components_on_verified_set(env, gil)


if __name__ == "__main__":
    main(my_env, parsed_args)
