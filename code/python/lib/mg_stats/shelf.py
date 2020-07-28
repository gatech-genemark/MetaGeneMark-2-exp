# Author: Karl Gemayel
# Created: 7/23/20, 8:43 PM

import logging
from typing import *

log = logging.getLogger(__name__)


def all_columns_equal(df, columns=None):
    # type: (pd.DataFrame, List[str]) -> pd.Series
    """Return True/False series shown rows where all columns have the same value"""

    if columns is None:
        columns = df.columns.values

    # create condition list
    conditions = list()
    for i in range(1, len(columns)):
        conditions.append(f"(df[{columns[i-1]}] == df[{columns[i]}])")

    return eval(" & ".join(conditions))


def create_joint_reference_from_list(df, list_reference):
    # type: (pd.DataFrame, List[str]) -> str

    reference = "=".join(list_reference)

    reference_rows = all_columns_equal(df, [f"'5p-{r}'" for r in list_reference])
    reference_values = df.loc[reference_rows, f"5p-{list_reference[0]}"]
    df.loc[reference_rows, f"5p-{reference}"] = reference_values

    reference_rows = all_columns_equal(df, [f"'3p-{r}'" for r in list_reference])
    reference_values = df.loc[reference_rows, f"3p-{list_reference[0]}"]
    df.loc[reference_rows, f"3p-{reference}"] = reference_values

    return reference