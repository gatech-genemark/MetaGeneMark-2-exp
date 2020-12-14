import logging
from typing import *

from mg_container.shelf import read_value_for_tag, convert_to_matrix, convert_to_position_distribution, \
    gms2_model_matrix_to_string, gms2_model_position_distribution_to_string
from mg_io.general import write_to_file

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

    def to_string(self):
        # type: () -> str

        ordered_tags = ["NAME", "GCODE", "NON_DURATION_DECAY", "COD_DURATION_DECAY", "COD_P_N", "NON_P_N", "ATG", "GTG",
                        "TTG", "TAA", "TAG", "TGA", "GENE_MIN_LENGTH", "NON_ORDER", "COD_ORDER", "NON_MAT", "COD_MAT"]

        remaining_tags = sorted(set(self.items.keys()).difference(ordered_tags))

        out = "__NATIVE\n"

        for tag in ordered_tags + remaining_tags:

            out += f"${tag}"
            if tag in self.items.keys():
                if tag.endswith("_MAT"):
                    out += "\n"
                    val_string = gms2_model_matrix_to_string(self.items[tag])
                elif tag.endswith("POS_DISTR"):
                    out += "\n"
                    val_string = gms2_model_position_distribution_to_string(self.items[tag])
                else:
                    out += " "
                    val_string = str(self.items[tag]) + "\n"

                out += f"{val_string}"

        return out

    def to_file(self, pf_output):
        # type: (str) -> None
        write_to_file(self.to_string(), pf_output)


