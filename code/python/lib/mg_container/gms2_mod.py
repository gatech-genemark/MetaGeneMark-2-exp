import logging
from typing import *

from mg_container.shelf import read_value_for_tag, convert_to_matrix, convert_to_position_distribution

log = logging.getLogger(__name__)


class GMS2Mod:

    def __init__(self, items):
        # type: (Dict[str, Any]) -> None
        self.items = items

    @classmethod
    def init_from_file(cls, pf_mod):
        # type: (str) -> GMS2Mod
        try:
            f = open(pf_mod, "r")
        except FileNotFoundError:
            raise ValueError(f"Can't open file {pf_mod}")

        words = [word for line in f for word in line.split()]
        result = dict()
        position = 0
        while position < len(words) - 1:
            curr_word = words[position]
            if curr_word.startswith("$"):
                tag = curr_word.split("$", maxsplit=1)[1]
                value, position = read_value_for_tag(words, position + 1)

                if len(value) == 1:
                    result[tag] = value[0]
                else:
                    if tag.endswith("_MAT"):
                        result[tag] = convert_to_matrix(value)
                    elif tag.endswith("_POS_DISTR"):
                        result[tag] = convert_to_position_distribution(value)
                    else:
                        log.warning(f"Unknown format for tag: {tag}")
            else:
                pass
                position += 1
                # raise ValueError("Error in reading file")

        return cls(result)


