import itertools
import logging
import pandas as pd
from typing import *

from mg_general.general import get_value
from mg_io.general import remove_p

log = logging.getLogger(__name__)

T = TypeVar('T')


def merge_identity(list_output_data, **kwargs):
    # type: (List[T]) -> List[T]
    return list_output_data


def merge_dataframes(list_df):
    # type: (List[pd.DataFrame]) -> pd.DataFrame
    return pd.concat(list_df, ignore_index=True, sort=False)

def merge_dataframes_to_file(dfs, **kwargs):
    # type: (Iterable[pd.DataFrame], Dict[str, Any]) -> None
    pf_output = get_value(kwargs, "pf_output", required=True)
    remove_p(pf_output)
    counter = 0
    header = None
    for df in dfs:
        df.sort_index(axis=1, inplace=True)
        if header == None:
            header = list(df.columns.values)
        else:

            if header != list(df.columns.values):
                logger.warning("Could not append dataframe to file. Header inconsistent")
                continue

        print(len(sorted(df.columns.values)))
        df.to_csv(pf_output, index=False, mode="a", header=counter==0)
        counter += 1


def merge_lists(list_lists):
    # type: (List[List[T]]) -> List[T]
    merged = list()
    for l in list_lists:
        merged += l

    return merged
