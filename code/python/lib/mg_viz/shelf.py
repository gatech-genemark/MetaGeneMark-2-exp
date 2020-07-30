# Author: Karl Gemayel
# Created: 7/8/20, 9:22 AM

import logging
import os
from tempfile import mkstemp
from matplotlib.ticker import FuncFormatter

import pandas as pd
from typing import *

log = logging.getLogger(__name__)


def create_mappable_for_colorbar(values, cmap):
    import matplotlib.colors
    import matplotlib.cm
    vmin = min(values)
    vmax = max(values)

    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array(values)
    return mappable


def get_order_by_rank(df, col_id, col_value, col_hue):
    # type: (pd.DataFrame, str, str, str) -> List[str]
    """Returns the order of hues based on which have the highest value
    """

    # Note: Orders are computed in increasing order.

    unique_hues = df[col_hue].unique()

    total_score_per_hue = {h: 0 for h in unique_hues}

    for val_id, df_group in df.groupby(col_id, as_index=False):  # type: Any, pd.DataFrame

        hue_order_for_id = list(df_group.sort_values(col_value)[col_hue])

        # add any missing hues to the start (since these are assumed to have value -inf)
        empty_hues = list(set(unique_hues).difference(hue_order_for_id))
        hue_order_for_id = empty_hues + hue_order_for_id

        # assign scores. The higher in rank, the more score it gets
        for i in range(len(hue_order_for_id)):
            total_score_per_hue[hue_order_for_id[i]] += i

    # return based on order
    return [a[0] for a in sorted(total_score_per_hue.items(), key=lambda x: x[1])]


def mkstemp_closed(**kwargs):
    # type: (Dict[str, Any]) -> str
    fd, pf = mkstemp(**kwargs)
    os.close(fd)
    return pf


def number_formatter(number, pos=None):
    """Convert a number into a human readable format."""
    magnitude = 0
    while abs(number) >= 1000:
        magnitude += 1
        number /= 1000.0
    return '%.0f%s' % (number, ['', 'K', 'M', 'B', 'T', 'Q'][magnitude])
