# Author: Karl Gemayel
# Created: 8/5/20, 8:06 AM

import logging
import numpy as np
import pandas as pd
from typing import *

from mg_general import Environment
from mg_parallelization.generic_threading import run_one_per_thread
from mg_stats.shelf import _helper_df_joint_reference, update_dataframe_with_stats, tidy_genome_level

log = logging.getLogger(__name__)


def get_stats_at_gcfid_level_with_reference(df, tools, reference):
    # type: (pd.DataFrame, List[str], str) -> pd.DataFrame

    list_entries = list()

    for gcfid, df_group in df.groupby("Genome", as_index=False):

        result = dict()
        for t in tools:

            tag = ",".join([t, reference])
            tag_eq = "=".join([t, reference])

            if df_group[f"3p:Match({tag_eq})"].sum() == 0:
                result[f"Match({tag})"] = np.nan
                result[f"Number of Error({tag})"] = np.nan
                result[f"Number of Found({tag})"] = np.nan
                result[f"Number of Predictions({tag})"] = np.nan
            else:
                result[f"Match({tag})"] = 100 * df_group[f"5p:Match({tag_eq})"].sum() / float(
                    df_group[f"3p:Match({tag_eq})"].sum())
                result[f"Number of Error({tag})"] = df_group[f"3p:Match({tag_eq})"].sum() - df_group[
                    f"5p:Match({tag_eq})"].sum()
                result[f"Number of Match({tag})"] = df_group[f"5p:Match({tag_eq})"].sum()

                result[f"Number of Found({tag})"] = df_group[f"3p:Match({tag_eq})"].sum()
                result[f"Number of Predictions({tag})"] = df[f"3p-{t}"].count()

                result[f"Sensitivity({t},{reference})"] = result[f"Number of Found({t},{reference})"] / df_group[
                    f"5p-{reference}"].count()
                result[f"Specificity({t},{reference})"] = result[f"Number of Found({t},{reference})"] / df_group[
                    f"5p-{t}"].count()

        result["Genome"] = gcfid
        result["Genome GC"] = df_group.at[df_group.index[0], "Genome GC"]
        result["Number in Reference"] = df_group[f"5p-{reference}"].count()

        list_entries.append(result)

    return pd.DataFrame(list_entries)


def _helper_join_reference_and_tidy_data(env, df_per_gene, tools, list_ref):
    # type: (Environment, pd.DataFrame, List[str], List[str]) -> [str, pd.DataFrame]

    reference = _helper_df_joint_reference(df_per_gene, list_ref)
    df_per_gene = update_dataframe_with_stats(df_per_gene, tools, reference).copy()

    #### Genome Level: compute stats per genome
    df_stats_gcfid = list()
    for _, df_group in df_per_gene.groupby("Genome", as_index=False):
        df_stats_gcfid.append(get_stats_at_gcfid_level_with_reference(df_group, tools, reference))
    df_per_genome = pd.concat(df_stats_gcfid, ignore_index=True, sort=False)

    ### Tidy Data and filter out those not present in tools or reference
    df_tidy = tidy_genome_level(env, df_per_genome)
    df_tidy = df_tidy[df_tidy["Tool"].apply(lambda x: x.lower()).isin(tools + [reference])]

    log.debug(f"Parsing: {', '.join(df_tidy['Genome'].unique())}")

    return reference, df_tidy

def prl_join_reference_and_tidy_data(env, df_per_gene, tools, list_ref):
    def yield_data(gb_gen):
        counter = 0
        for a in gb_gen:
            log.debug(f"Parsing: {', '.join(a[1]['Genome'].unique())}, {counter}")
            counter += 1
            yield a[1]
    list_ref_df = run_one_per_thread(yield_data(df_per_gene.groupby("Genome", as_index=False)),
                       _helper_join_reference_and_tidy_data,
                       "df_per_gene", {"env": env, "tools": tools, "list_ref": list_ref},
                       simultaneous_runs=7)

    if len(list_ref_df) > 0:
        list_df_genome = [x[1] for x in list_ref_df]
        reference = list_ref_df[0][0]
        return reference, pd.concat(list_df_genome, ignore_index=True, sort=False)
    else:
        return [None, None]

