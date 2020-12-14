import logging
from typing import *
from intervaltree import Interval, IntervalTree

from mg_models.mgm_motif_model import MGMMotifModel
from mg_models.mgm_motif_model_v2 import MGMMotifModelV2

log = logging.getLogger(__name__)


class MGMMotifModelAllGC:
    """Class that holds all MGM motif models, and indexes them by their GC range"""

    def __init__(self, models):
        # type: (List[List[float, float, MGMMotifModelV2]]) -> None

        self._models = MGMMotifModelAllGC._create_interval_structure(models)

    @staticmethod
    def _create_interval_structure(models):
        # type: (List[List[float, float, MGMMotifModelV2]]) -> IntervalTree

        sorted_models = sorted(models, key=lambda x: x[0])
        sorted_models[0][0] = -float('inf')
        sorted_models[len(sorted_models)-1][1] = float('inf')
        list_intervals = [
            Interval(x[0], x[1], x[2]) for x in sorted_models
        ]
        return IntervalTree(list_intervals)

    def get_model_by_gc(self, gc):
        # type: (float) -> MGMMotifModelV2
        # if len(self._models[gc]) == 0:
        #     print(gc)
        #     print(self._models)
        # self._models[gc].pop().data
        return self._models[gc].pop().data
