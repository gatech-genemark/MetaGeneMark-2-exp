import logging
from typing import *

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
                value, position = GMS2Mod._read_value(words, position + 1)

                if len(value) == 1:
                    result[tag] = value[0]
                else:
                    if tag.endswith("_MAT"):
                        result[tag] = GMS2Mod._convert_to_matrix(value)
                    elif tag.endswith("_POS_DISTR"):
                        result[tag] = GMS2Mod._convert_to_position_distribution(value)
                    else:
                        log.warning(f"Unknown format for tag: {tag}")
            else:
                pass
                position += 1
                # raise ValueError("Error in reading file")

        return cls(result)

    @staticmethod
    def _convert_to_position_distribution(words):
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

    @staticmethod
    def _convert_to_matrix(words):
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

    @staticmethod
    def _read_value(words, position):
        # type: (List[str], int) -> Tuple[List[str], int]

        num_words = len(words)
        result = list()
        while position < num_words:
            curr_word = words[position]
            if curr_word.startswith("$"):
                break

            result.append(curr_word)

            position += 1

        return result, position
