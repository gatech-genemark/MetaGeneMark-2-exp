# Author: Karl Gemayel
# Created: 7/14/20, 9:41 AM

import logging
from typing import *
from mg_general import Environment
from mg_general.general import os_join
from mg_options.options import Options

log = logging.getLogger(__name__)


class LearnFromOptions(Options):
    """Maps which group of genomes to learn from, for each feature
    E.g. It can specify that Start codons for archaea group A should be learned only from
    group A, or from group A and D.

    For now, allowed groups are hard coded: Bacteria (A, B, C, X) and Archaea (A, D)
    """

    def __init__(self, env, pf_options_custom, **kwargs):
        # type: (Environment, str, Dict[str, Any]) -> None
        super(LearnFromOptions, self).__init__(env, pf_options_custom, **kwargs)

    def path_to_default_options_file(self, env):
        # type: (Environment) -> str
        return os_join(env["pd-config"], "learn_from_defaults.conf")

    def required(self):  # type: () -> Union[Set[str], None]
        return None
# class LearnFrom:
#
#
#     def __init__(self, info=None):
#         # type: (Dict[str, Dict[str, Dict[str, Set[str]]]]) -> None
#         self._values = LearnFrom._default_values()
#
#         if info is not None:
#             # update values based on input
#             for component, c_vals in self._values.items():
#                 if component in info:
#                     for gtype, g_vals in c_vals.items():
#                         if gtype in info[component]:
#                             for group, gr_vals in g_vals.items():
#                                 if group in info[component][gtype]:
#                                     self._values[component][gtype][group] = info[component][gtype][group]

    # def __getitem__(self, item):
    #     # type: (str) -> Dict[str, Dict[str, Set[str]]]
    #     return self._values[item]
    #
    # @staticmethod
    # def _default_values():
    #     # type: () -> Dict[str, Dict[str, Dict[str, Set[str]]]]
    #     return {
    #         component: {
    #             "Archaea": {"A": {"A"}, "D": {"D"}},
    #             "Bacteria": {"A": {"A"}, "B": {"B"}, "C": {"C"}, "X": {"X"}}
    #         } for component in {"RBS", "PROMOTER", "Start Context", "Start Codons", "Stop Codons"}
    #     }
    #
    # @classmethod
    # def init_from_file(cls, pf_config):
    #     # type: (str) -> LearnFrom
    #     try:
    #         f = open(pf_config, "r")
    #         return LearnFrom(yaml.load(f, Loader=yaml.FullLoader))
    #     except IOError:
    #         logger.warning(f"Configuration File Not Found: {pf_config}. "
    #                        f"Using defaults.")
    #         return LearnFrom()