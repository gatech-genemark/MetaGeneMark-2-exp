import itertools
import logging
import pandas as pd
from typing import *

log = logging.getLogger(__name__)

T = TypeVar('T')


def merge_identity(list_output_data):
    # type: (List[T]) -> List[T]
    return list_output_data


def merge_dataframes(list_df):
    # type: (List[pd.DataFrame]) -> pd.DataFrame
    return pd.concat(list_df, ignore_index=True, sort=False)


def merge_lists(list_lists):
    # type: (List[List[T]]) -> List[T]

    return list(itertools.chain(*list_lists))
