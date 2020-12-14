import logging
import math
import pandas as pd
from typing import *

from mg_general.general import get_value
from mg_models.motif_model import MotifModel

log = logging.getLogger(__name__)


class MGMMotifModel(MotifModel):
    """Motif model, holding the composition matrix and position (spacer) distribution"""

    def __init__(self, shift_prior, motif, motif_width, spacer=None, **kwargs):
        # type: (Dict[int, float], Dict[str, List[float]], int, Dict[int, Dict[int, float]], Dict[str, Any]) -> None

        super().__init__(motif)
        self._motif_width = motif_width
        self._shift_prior = shift_prior

        self._spacer = MGMMotifModel._init_spacer(spacer)  # type: Union[None, Dict[int, List[float]]]

        self._kwargs = kwargs.copy()

    def score(self, fragment, **kwargs):
        # type: (str, Dict[str, Any]) -> float

        begin = get_value(kwargs, "begin", None)
        use_log = get_value(kwargs, "use_log", False)
        component = get_value(kwargs, "component", "both", choices=["both", "motif", "spacer"])
        prior = get_value(kwargs, "prior", True)

        if begin is None and len(fragment) != self._motif_width:
            raise ValueError("If 'begin' not specified, fragment length should equal motif width")
        elif begin is not None and begin + self._motif_width > len(fragment):
            raise ValueError("Not enough space in fragment")

        if begin is None:
            begin = 0

        score_per_shift = list()
        for s in self._shift_prior:
            s = int(s)
            # shift prior
            score = 0 if use_log else 1
            if prior:
                score = math.log(self._shift_prior[s]) if use_log else self._shift_prior[s]

            # motif
            if component != "spacer":
                for i in range(self._motif_width):
                    if fragment[begin + i] == "N":
                        if use_log:
                            score += 0.25
                        else:
                            score *= 0.25
                    else:
                        if use_log:
                            score += math.log(self._motif[fragment[begin + i]][s + i])
                        else:
                            score *= self._motif[fragment[begin + i]][s + i]

            # spacer
            if component != "motif":
                if self._spacer is not None and s in self._spacer.keys():
                    distance_from_start = len(fragment) - (begin + self._motif_width)
                    if use_log:
                        score += math.log(self._spacer[s][distance_from_start])
                    else:
                        score *= self._spacer[s][distance_from_start]

            score_per_shift.append(score)

        return max(score_per_shift)

    def find_best_position_and_score(self, fragment, **kwargs):
        # type: (str, Dict[str, Any]) -> Tuple[int, float, float]

        v_list = [
            (pos, self.score(fragment, begin=pos, **kwargs), self.score(fragment, begin=pos, prior=False, **kwargs)) for
            pos in range(len(fragment) - self._motif_width)]
        return max(
            v_list,
            key=lambda x: x[1]
        )

    @staticmethod
    def _init_spacer(spacer=None):
        # type: (Union[None, Dict[int, Dict[int, float]]], Dict[int, List[float]]) -> Union[Dict[int, List[float]], None]

        if spacer is None:
            return None

        new_spacer = dict()
        for shift in spacer.keys():
            shift = int(shift)

            if isinstance(spacer[shift], dict):
                if len(spacer[shift].keys()) == 0:
                    continue
                max_position = max([int(x) for x in spacer[shift].keys()])
                result = [0] * (max_position + 1)
                for i in spacer[shift].keys():
                    result[int(i)] = spacer[shift][i]

                new_spacer[shift] = result

            elif isinstance(spacer[shift], list):
                new_spacer[shift] = spacer[shift].copy()

            else:
                raise ValueError(f"Unknown spacer type: {type(spacer[shift])}")

        return new_spacer

    def pwm_to_df(self):
        # type: () -> pd.DataFrame

        keys = sorted(self._motif.keys())

        list_entries = list()
        for p in range(len(next(iter(self._motif.values())))):
            list_entries.append(
                [self._motif[k][p] for k in keys]
            )

        return pd.DataFrame(list_entries, columns=keys)
