import itertools
import logging
import pandas as pd
from typing import *

from mg_general.general import get_value
from mg_io.general import remove_p

log = logging.getLogger(__name__)

T = TypeVar('T')


def merge_identity(list_output_data):
    # type: (List[T]) -> List[T]
    return list_output_data


def merge_dataframes(list_df):
    # type: (List[pd.DataFrame]) -> pd.DataFrame
    return pd.concat(list_df, ignore_index=True, sort=False)

def merge_dataframes_to_file(dfs, **kwargs):
    # type: (Iterable[pd.DataFrame], Dict[str, Any]) -> None
    pf_output = get_value(kwargs, "pf_output", required=True)
    remove_p(pf_output)

    for df in dfs:
        df.to_csv(pf_output, index=False, mode="a")


def merge_lists(list_lists):
    # type: (List[List[T]]) -> List[T]
    merged = list()
    for l in list_lists:
        merged += l

    return merged
