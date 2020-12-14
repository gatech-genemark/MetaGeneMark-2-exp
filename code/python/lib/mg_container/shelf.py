# Author: Karl Gemayel
# Created: 6/22/20, 2:27 PM

import logging
from typing import *

from mg_general.general import get_value

log = logging.getLogger(__name__)


def read_value_for_tag(words, position, **kwargs):
    # type: (List[str], int, Dict[str, Any]) -> Tuple[List[str], int]

    stop_if_starts_with = get_value(kwargs, "stop_if_starts_with", {"$"}, valid_type=Set[str])

    num_words = len(words)
    result = list()
    while position < num_words:
        curr_word = words[position]

        should_stop = False
        for s in stop_if_starts_with:
            if curr_word.startswith(s):
                should_stop = True
                break

        if should_stop:
            break

        result.append(curr_word)

        position += 1

    return result, position


def convert_to_position_distribution(words):
    # type: (List[str]) -> Dict[int, float]

    num_words = len(words)
    result = dict()

    if num_words % 2 != 0:
        raise ValueError("Position distribution should have equal number of positions and probabilities")

    for index in range(0, num_words, 2):
        pos = words[index]
        prob = words[index + 1]

        try:
            float_prob = float(prob)
            int_pos = int(pos)
            result[int_pos] = float_prob
        except ValueError:
            raise ValueError(f"Unknown value/probability pair: {pos, prob}")

    return result


def convert_to_matrix(words):
    # type: (List[str]) -> Dict[str, List[float]]

    num_words = len(words)
    result = dict()
    key = None
    for position in range(num_words):
        curr_word = words[position]
        if curr_word.startswith("$"):
            break

        try:
            float_word = float(curr_word)
            # number
            if key is None:
                raise ValueError(f"Readingn value {curr_word} without key")

            result[key].append(float_word)
        except ValueError:

            key = curr_word

            if key in result:
                raise ValueError(f"Reading same key multiple times {key}")

            result[key] = list()

    return result


def gms2_model_matrix_to_string(value):
    # type: (Dict[str, List[float]]) -> str

    out = ""
    for key in sorted(value.keys()):
        out += key
        for v in value[key]:
            out += f" {v:.5f}"

        out += "\n"

    return out


def gms2_model_position_distribution_to_string(value):
    # type: (List[float]) -> str

    out = ""

    for i, v in enumerate(value):
        out += f"{i} {v}\n"

    return out


