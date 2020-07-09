# Author: Karl Gemayel
# Created: 7/8/20, 9:22 AM

import logging
from typing import *

log = logging.getLogger(__name__)

def create_mappable_for_colorbar(values, cmap):
    import matplotlib.colors
    import matplotlib.cm
    vmin=min(values)
    vmax = max(values)

    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    mappable = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array(values)
    return mappable


