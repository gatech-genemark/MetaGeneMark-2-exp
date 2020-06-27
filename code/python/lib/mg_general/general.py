# Author: karl
# Created: 2020-06-21, 9:03 a.m.
import copy
import logging
import os
import subprocess
import pandas as pd
from typing import *

log = logging.getLogger(__name__)


def _get_value_helper(kv_pairs, key, default, value_type=None):
    # type: (Dict[str, Any], str, Any, Type) -> Any
    value = kv_pairs[key] if key in kv_pairs else default

    if value is not None and value_type is not None:
        try:
            value = value_type(value)
        except ValueError:
            raise ValueError(f"Key {key} has value {value} of type {type(value)}. Must have type {value_type}.")

    return value


def get_value(kv_pairs, key, default=None, **kwargs):
    # type: (Dict[str, Any], str, Any, Dict[str, Any]) -> Any

    value_type = _get_value_helper(kwargs, "type", None)  # TODO: figure out value_type for this
    choices = _get_value_helper(kwargs, "choices", None, value_type=list)
    required = _get_value_helper(kwargs, "required", False, value_type=bool)
    perform_copy = _get_value_helper(kwargs, "copy", False, value_type=bool)
    default_if_none = _get_value_helper(kwargs, "default_if_none", True, value_type=bool)
    default_value_callable = _get_value_helper(kwargs, "default_value_callable", None, value_type=Callable)

    value = _get_value_helper(kv_pairs, key, default)
    if value is None and default_if_none:
        if default is not None:
            value = default
        elif default_value_callable is not None:
            value = default_value_callable()

    # perform Checks
    if value is not None and value_type is not None and not isinstance(value, value_type):
        try:
            value = value_type(value)
        except ValueError:
            raise ValueError(f"Key {key} has value {value} of type {type(value)}. Must have type {value_type}.")

    if required and not default_if_none and value is None:
        raise ValueError(f"Key {key} is required.")

    if choices is not None and value not in choices:
        raise ValueError(f"Key {key} has value {value}, not in: " + ", ".join(choices))



    if perform_copy:
        value = copy.deepcopy(value)

    return value


def except_if_not_valid(value, allowed_choices=None, invalid_choices=None):
    # type: (Any, Iterable[Any], Iterable[Any]) -> None
    if allowed_choices is not None:
        if value not in allowed_choices:
            raise ValueError(f"Value {value} not in allowed set: " + ", ".join(allowed_choices))

    if invalid_choices is not None:
        if value in invalid_choices:
            raise ValueError(f"Value {value} not in invalid set: " + ", ".join(invalid_choices))


def run_shell_cmd(cmd, do_not_log=False):
    # type: (str, bool) -> str
    if not do_not_log:
        log.debug(cmd)

    return subprocess.check_output(cmd, shell=True).decode("utf-8")


def os_join(*args):
    # type: (List[Any]) -> str
    return os.path.join(*args)


def all_elements_equal(elements):
    # type: (Iterable[Any]) -> bool

    old_val = None
    already_set_first = False
    for v in elements:
        if already_set_first:
            if v != old_val:
                return False
        else:
            old_val = v
            already_set_first = True

    return True


def next_name(pd_work, **kwargs):
    # type: (str, Dict[str, Any]) -> str

    ext = get_value(kwargs, "ext", "pdf")
    if "counter" not in next_name.__dict__: next_name.counter = -1
    next_name.counter += 1
    return os_join(pd_work, "{}.{}".format(next_name.counter, ext))


def create_gene_key(genome=None, accession=None, left=None, right=None, strand=None, delimiter=";"):
    # type: (object, object, object, object, object, str) -> str

    return "{}{}{}{}{}{}{}{}{}".format(
        genome, delimiter,
        accession, delimiter,
        left, delimiter,
        right, delimiter,
        strand
    )

def fix_names(r):
    # type: (pd.Series) -> str
    return "{}. {}".format(
        r["Genome"][0], r["Genome"].split("_")[1]
    )