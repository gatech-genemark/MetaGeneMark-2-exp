# Author: karl
# Created: 2020-06-21, 8:31 a.m.

import os
import copy
import logging
import argparse
from typing import *

logger = logging.getLogger(__name__)


def add_env_args_to_parser(parser):
    # type: (argparse.ArgumentParser) -> None

    parser.add_argument('--pd-work', required=False, default=None, help="Path to working directory")
    parser.add_argument('--pd-data', required=False, default=None, help="Path to data directory")
    parser.add_argument("-l", "--log", dest="loglevel", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                        help="Set the logging level", default='WARNING')


class Environment:
    """A class representing the project environment variables, including
    paths to useful directories (data, code, etc...)
    """

    def __init__(self, pd_data=None, pd_work=None, **kwargs):
        # type: (str, str, Dict[str, Any]) -> None

        self._env = Environment.load_environment_variables(pd_data, pd_work, **kwargs)

    def __getitem__(self, item):
        # type: (str) -> Any
        return self._env[item]

    def __setitem__(self, key, value):
        # type: (str, Any) -> None
        self._env[key] = value

    def duplicate(self, new_values=None):
        # type: (Dict[str, Any]) -> Environment
        """Creates a copy of the environment, with update variables
        """
        new_env = copy.deepcopy(self)
        if new_values is not None:
            for item in new_values.keys():
                new_env[item] = new_values[item]

        return new_env

    @classmethod
    def init_from_argparse(cls, parser):
        # type: (argparse.Namespace) -> Environment

        return cls(**Environment.load_environment_variables(
            parser.pd_data, parser.pd_work))

    @staticmethod
    def load_environment_variables(pd_data=None, pd_work=None, **kwargs):
        # type: (str, str, Dict[str, Any]) -> Dict[str, str]

        # path to current file
        pd_current_file = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

        # Structure:
        pd_base = os.path.abspath(os.path.join(pd_current_file, "../../../../"))        # get base level
        pd_bin = os.path.join(pd_base, "bin")
        pd_tmp = os.path.join(pd_base, "tmp")
        pd_runs = os.path.join(pd_base, "runs")
        pd_code = os.path.join(pd_base, "code")
        pd_config = os.path.join(pd_base, "config")
        pd_bin_external = os.path.join(pd_base, "bin_external")

        pd_data = os.path.abspath(pd_data) if pd_data is not None else os.path.join(pd_base, "data")
        pd_work = os.path.abspath(pd_work) if pd_work is not None else os.path.join(pd_base, ".")

        if not os.path.exists(pd_work):
            os.makedirs(pd_work)

        env = {
            "pd-base": pd_base, "pd-bin": pd_bin, "pd-tmp": pd_tmp, "pd-runs": pd_runs, "pd-code": pd_code,
            "pd-config": pd_config, "pd-bin-external": pd_bin_external, "pd-data": pd_data, "pd-work": pd_work
        }

        import copy
        global ENV
        ENV = copy.deepcopy(env)

        return env


ENV = Environment()
