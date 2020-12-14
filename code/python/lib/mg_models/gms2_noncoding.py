import logging
import numpy as np
from typing import *

log = logging.getLogger(__name__)


class GMS2Noncoding:

    def __init__(self, pwm):
        # type: (Dict[str, List[float]]) -> None

        pwm_no_list = {
            x: pwm[x][0] for x in pwm
        }
        self._pwm_by_order = GMS2Noncoding._get_all_orders(pwm_no_list)
        self._max_order = len(next(iter(pwm.keys())))

    @staticmethod
    def _get_all_orders(pwm):
        # type: (Dict[str, float]) -> Dict[int, Dict[str, float]]

        max_order = len(next(iter(pwm.keys()))) - 1

        result = dict()

        result[max_order] = pwm.copy()

        curr_order = max_order - 1
        while curr_order >= 0:

            result[curr_order] = GMS2Noncoding._lower_order_by_one(result[curr_order+1])

            curr_order -= 1

        return result

    @staticmethod
    def _lower_order_by_one(pwm):
        # type: (Dict[str, float]) -> Dict[str, float]

        result = dict()

        for key, value in pwm.items():

            key_marginal = key[1:]
            if key_marginal not in result.keys():
                result[key_marginal] = 0

            result[key_marginal] += value

        return result

    def prob(self, key):
        # type: (str) -> float

        if len(key) > self._max_order + 1:
            raise ValueError(f"Maximum order is {self._max_order}. Can't process: {key}")

        return self._pwm_by_order[len(key)-1][key]

    def pwm_to_array(self, order):
        # type: (int) -> np.ndarray
        pwm = self._pwm_by_order[order]
        sorted_keys = sorted(pwm.keys())

        arr = np.zeros(len(sorted_keys), dtype=float)
        for i, k in enumerate(sorted_keys):
            arr[i] = pwm[k]

        return arr
