# Author: Karl Gemayel
# Created: 6/22/20, 1:26 PM

import logging
from typing import *

from mg_container.shelf import read_value_for_tag, convert_to_position_distribution, convert_to_matrix, \
    gms2_model_matrix_to_string, gms2_model_position_distribution_to_string
from mg_io.general import write_to_file

log = logging.getLogger(__name__)


class MGMModelGC:
    def __init__(self, items, gc):
        # type: (Dict[str, Any], float) -> None

        self.items = items
        self.gc = gc

    def to_string(self):
        # type: () -> str

        ordered_tags = ["NAME", "GCODE", "NON_DURATION_DECAY", "COD_DURATION_DECAY", "COD_P_N", "NON_P_N", "ATG", "GTG",
                        "TTG", "TAA", "TAG", "TGA", "GENE_MIN_LENGTH", "NON_ORDER", "COD_ORDER", "NON_MAT", "COD_MAT"]

        remaining_tags = sorted(set(self.items.keys()).difference(ordered_tags))

        out = ""

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


class MGMModel:
    def __init__(self, items_by_species_and_gc):
        # type: (Dict[str, Dict[float, MGMModelGC]]) -> None

        self.items_by_species_and_gc = items_by_species_and_gc  # type: Dict[str, Dict[float, MGMModelGC]]

    def to_string(self):
        # type: () -> str

        out = ""
        for species in self.items_by_species_and_gc.keys():

            for gc, model_gc in self.items_by_species_and_gc[species].items():
                out += f"__{species}_GC_{int(gc)}\n"

                out += model_gc.to_string()

        return out

    def to_file(self, pf_output):
        # type: (str) -> None
        write_to_file(self.to_string(), pf_output)

    @classmethod
    def init_from_file(cls, pf_mod):
        # type: (str) -> MGMModel

        # GC bins are split by __1_GC_2 tags, where
        # 1: B for bacteria, A for archaea
        # 2: Integer showing GC bin

        try:
            f = open(pf_mod, "r")
        except FileNotFoundError:
            raise ValueError(f"Can't open file {pf_mod}")

        words = [word for line in f for word in line.split()]
        result = dict()
        position = 0

        while position < len(words) - 1:

            curr_word = words[position]

            if curr_word.startswith("__"):
                species, _, gc = curr_word[2:].split("_")
                position += 1
                mgm_model_gc, position = MGMModel._read_mgm_model_gc(words, position, gc)

                if species not in result.keys():
                    result[species] = dict()

                result[species][gc] = mgm_model_gc
            else:
                position += 1

        return cls(result)

    @staticmethod
    def _read_value(words, position):
        # type: (List[str], int) -> Tuple[List[str], int]

        """Read all words until next tag that starts with a '$' sign

        Return a list of words start from 'position', up until but excluding the next tag.
        Also returned is the position of the next tag
        """

        num_words = len(words)
        result = list()

        while position < num_words:
            curr_word = words[position]
            if curr_word.startswith("$"):
                break

            result.append(curr_word)
            position += 1

        return result, position

    @staticmethod
    def _read_mgm_model_gc(words, position, gc):
        # type: (List[str], int, float) -> [MGMModelGC, int]

        """Read all tags up to next GC tag, and put in MGMMOdelGC"""
        result = dict()

        while position < len(words) - 1:
            curr_word = words[position]
            if curr_word.startswith("$"):
                tag = curr_word.split("$", maxsplit=1)[1]
                value, position = read_value_for_tag(words, position + 1, stop_if_starts_with={"$", "__"})

                if len(value) == 1:
                    result[tag] = value[0]
                else:
                    if tag.endswith("_MAT"):
                        result[tag] = convert_to_matrix(value)
                    elif tag.endswith("_POS_DISTR"):
                        result[tag] = convert_to_position_distribution(value)
                    else:
                        log.warning(f"Unknown format for tag: {tag}")
            elif curr_word.startswith("__"):
                break
            else:
                position += 1

        return [MGMModelGC(result, gc), position]
