import os
import logging
from typing import *
from mg_general import Environment
from mg_general.general import os_join
from mg_options.options import Options

log = logging.getLogger(__name__)


class ParallelizationOptions(Options):
    """Options for parallelization of program executions, including PBS systems"""

    def __init__(self, env, pf_options_custom, **kwargs):
        # type: (Environment, str, Dict[str, Any]) -> None
        super(ParallelizationOptions, self).__init__(env, pf_options_custom, **kwargs)

        self._check_pbs_paths_for_directories()

    def path_to_default_options_file(self, env):
        # type: (Environment) -> str
        return os_join(env["pd-config"], "parallelization_defaults.conf")

    def required(self):  # type: () -> Union[Set[str], None]
        return None

    def _check_pbs_paths_for_directories(self):
        # type: () -> None

        # if pbs head node directory not specified, use current working directory
        if self._options["pbs-pd-head"] is None:
            self._options["pbs-pd-head"] = self.env["pd-work"]

        # make sure it's an absolute path
        self._options["pbs-pd-head"] = os.path.abspath(self._options["pbs-pd-head"])

        # if path to compute not specified, use PBS head directory
        if self._options["pbs-pd-root-compute"] is None:
            self._options["pbs-pd-root-compute"] = self._options["pbs-pd-head"]

    def update_env(self, env):
        # type: (Environment) -> None
        self.env = env.duplicate()
        self._options["pbs-pd-head"] = os.path.abspath(self.env["pd-work"])
