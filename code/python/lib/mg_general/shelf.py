# Author: karl
# Created: 2020-06-21, 8:56 a.m.

import logging
from typing import *

log = logging.getLogger(__name__)


def test_log_level():
    log.debug(f"Test")
    log.info(f"Test")
    log.warning(f"Test")
    log.critical(f"Test")


def list_find_first(a_list, a_filter):
    # type: (List[Any], Callable) -> Any
    for x in a_list:
        if a_filter(x):
            return x
    return None
