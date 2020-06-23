# Karl Gemayel
# Georgia Institute of Technology
#
# Created:

import logging
import random
import argparse
from typing import *

# noinspection PyUnresolvedReferences
import pathmagic                        # add path to custom library

# Custom library imports
from mg_general import Environment
import mg_general
from mg_general.general import run_shell_cmd
from mg_io.general import mkdir_p
from mg_parallelization.pbs import PBSJobPackage

# ------------------------------ #
#           Parse CMD            #
# ------------------------------ #

parser = argparse.ArgumentParser("Description of driver.")

parser.add_argument('--pf-job-input', required=True)
parser.add_argument('--pf-job-output', required=True)


parser.add_argument('--pd-work', required=False, default=None, help="Path to working directory")
parser.add_argument('--pd-data', required=False, default=None, help="Path to data directory")
parser.add_argument('--pd-results', required=False, default=None, help="Path to results directory")
parser.add_argument("-l", "--log", dest="loglevel", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                    help="Set the logging level", default='WARNING')


parsed_args = parser.parse_args()

# ------------------------------ #
#           Main Code            #
# ------------------------------ #

# Load environment variables
my_env = Environment(pd_data=parsed_args.pd_data,
                     pd_work=parsed_args.pd_work,
                     pd_results=parsed_args.pd_results)

# Setup logger
logging.basicConfig(level=parsed_args.loglevel)
logger = logging.getLogger("logger")                    # type: logging.Logger


def main(env, args):
    # type: (Environment, argparse.Namespace) -> None

    pbs_package = PBSJobPackage.load(args.pf_job_input)
    func = pbs_package["func"]
    func_args = pbs_package["func_kwargs"]

    if "sbsp_options" in func_args:

        rs = func_args["sbsp_options"].safe_get("random-seed")
        if rs is None:
            random.seed(100)
        else:
            random.seed(int(rs))
            logger.critical("Random-seed: {}".format(rs))

    else:
        random.seed(100)

    if "env" in func_args:
        if args.pd_work is not None:
            func_args["env"] = func_args["env"].duplicate({"pd-work": args.pd_work})
            logger.critical("{}".format(func_args["env"]["pd-work"]))

    # Update pd-work to create a tmp directory
    mkdir_p(func_args["env"]["pd-work"])
    func_args["env"]["pd-work"] = run_shell_cmd(
        "mktemp --tmpdir={} -d".format(func_args["env"]["pd-work"])
    ).strip()

    # logger.critical("{}\n{}".format(func, func_args))
    output = {
        "data": func(**func_args)
    }

    PBSJobPackage.save(output, args.pf_job_output)


if __name__ == "__main__":
    main(my_env, parsed_args)
