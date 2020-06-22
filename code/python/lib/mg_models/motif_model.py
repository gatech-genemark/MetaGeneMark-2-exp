import logging
from typing import *
import pandas as pd

from sbsp_general.general import get_value

log = logging.getLogger(__name__)

class MotifModel:
    """Motif model, holding the composition matrix and position (spacer) distribution"""

    def __init__(self, motif, spacer=None):
        # type: (Dict[str, List[float]], Union[None, Dict[int, float]]) -> None

        self._motif = motif     # type: Dict[str, List[float]]
        self._motif_width = max([len(motif[x]) for x in self._motif.keys()])

        self._spacer = MotifModel._init_spacer(spacer)      # type: Union[None, List[float]]

    def score(self, fragment, **kwargs):
        # type: (str, Dict[str, Any]) -> float

        begin = get_value(kwargs, "begin", None)
        use_log = get_value(kwargs, "use_log", False)
        component = get_value(kwargs, "component", "both", choices=["both", "motif", "spacer"])

        if begin is None and len(fragment) != self._motif_width:
            raise ValueError("If 'begin' not specified, fragment length should equal motif width")
        elif begin is not None and begin + self._motif_width > len(fragment):
            raise ValueError("Not enough space in fragment")

        if begin is None:
            begin = 0

        score = 0 if use_log else 1
        if component != "spacer":
            for i in range(self._motif_width):
                if fragment[begin + i] == "N":
                    if use_log:
                        score += 0.25
                    else:
                        score *= 0.25
                else:
                    if use_log:
                        score += self._motif[fragment[begin + i]][i]
                    else:
                        score *= self._motif[fragment[begin + i]][i]

        if component != "motif":
            if self._spacer is not None:
                distance_from_start = len(fragment) - (begin + self._motif_width)
                if use_log:
                    score += self._spacer[distance_from_start]
                else:
                    score *= self._spacer[distance_from_start]

        return score

    def find_best_position_and_score(self, fragment, **kwargs):
        # type: (str, Dict[str, Any]) -> Tuple[int, float]

        return max(
            [(pos, self.score(fragment, begin=pos, **kwargs))
             for pos in range(len(fragment) - self._motif_width)],
            key=lambda x: x[1]
        )

    def motif_width(self):
        # type: () -> int
        return self._motif_width

    @staticmethod
    def _init_spacer(spacer=None):
        # type: (Union[None, Dict[int, float], List[float]]) -> Union[List[float], None]

        if spacer is None:
            return None

        if isinstance(spacer, dict):
            max_position = max([int(x) for x in spacer.keys()])
            result = [0] * (max_position + 1)
            for i in spacer.keys():
                result[int(i)] = spacer[i]

            return result

        if isinstance(spacer, list):
            return spacer.copy()

        raise ValueError(f"Unknown spacer type: {type(spacer)}")


    def pwm_to_df(self):
        # type: () -> pd.DataFrame

        keys = sorted(self._motif.keys())

        list_entries = list()
        for p in range(self.motif_width()):
            list_entries.append(
                [self._motif[k][p] for k in keys]
            )

        return pd.DataFrame(list_entries, columns=keys)

    def to_string(self):
        # type: () -> str
        out = ""
        for letter in sorted(self._motif.keys()):
            out += letter + " ".join([str(x) for x in self._motif[letter]]) + "\n"
        return out


