# Author: karl
# Created: 2020-06-21, 9:34 a.m.

import os
import errno
import logging
from typing import *

import cloudpickle as dill

log = logging.getLogger(__name__)


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def remove_p(*path):
    # type: (List[str]) -> None
    for p in path:
        if os.path.isfile(p):
            os.remove(p)


def write_to_file(a_string, pf_out, mode="w"):
    # type: (str, str, str) -> None
    with open(pf_out, mode) as f:
        f.write(a_string)


def save_obj(obj, name):
    # type: (object, str) -> None
    with open(name + '.pkl', 'wb') as f:
        dill.dump(obj, f)


def load_obj(name):
    with open(name + '.pkl', 'rb') as f:
        loaded_data = dill.load(f)
        return loaded_data