# Author: Karl Gemayel
# Created: 8/5/20, 8:54 AM

import argparse
import logging
from typing import *

logger = logging.getLogger(__name__)


def add_stats_options(parser):
    # type: (argparse.ArgumentParser) -> None
    parser.add_argument('--pf-data', required=True)
    parser.add_argument('--ref-3p', required=True, nargs='+', help="List of tools to be used as 3prime reference "
                                                                   "(ground-truth). Note: If more than one provided, "
                                                                   "their intersection (in terms of 5' end is taken as"
                                                                   " reference.")

    parser.add_argument('--ref-5p', required=True, nargs='+', help="List of tools to be used as 5prime reference "
                                                                   "(ground-truth). Note: If more than one provided, "
                                                                   "their intersection (in terms of 5' end is taken as"
                                                                   " reference.")

    parser.add_argument('--tools', nargs="+", help="If set, only compare these tools. Otherwise all tools are chosen")
    parser.add_argument('--parse-names', action='store_true', help="If set, try to shorten genome names. Useful only "
                                                                   "genome ID's in the data are actually names")
