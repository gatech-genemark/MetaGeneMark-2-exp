# Author: Karl Gemayel
# Created: 2020-06-21, 2:48 p.m.

import logging
from typing import *

from mg_general import Environment
from mg_general.general import os_join
from mg_options.options import Options

log = logging.getLogger(__name__)


class MGModelsOptions(Options):
    """Options for parallelization of program executions, including PBS systems"""

    def __init__(self, env, pf_options_custom, **kwargs):
        # type: (Environment, str, Dict[str, Any]) -> None
        super(MGModelsOptions, self).__init__(env, pf_options_custom, **kwargs)

    def path_to_default_options_file(self, env):
        # type: (Environment) -> str
        return os_join(env["pd-config"], "mg_models_defaults.conf")

    def required(self):  # type: () -> Union[Set[str], None]
        return None
