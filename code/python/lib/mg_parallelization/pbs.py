# Author: karl
# Created: 2020-06-21, 9:34 a.m.

import os
import copy

import time
import logging
import numpy as np
from typing import *

from mg_general import Environment
from mg_io.general import mkdir_p, write_to_file, remove_p

from mg_general.general import get_value, run_shell_cmd
from mg_options.parallelization import ParallelizationOptions
from mg_parallelization.pbs_job_package import PBSJobPackage

log = logging.getLogger(__name__)


class FunctionArguments:

    def __init__(self, **kwargs):
        self._kwargs = kwargs

    def get_arguments(self, data):
        # type: (Dict[str, Any]) -> Dict[str, Any]
        new_kwargs = copy.deepcopy(self._kwargs)
        new_kwargs.update(data)

        return new_kwargs


class PBS:
    """Runs any function on input data using PBS scheduler"""

    def __init__(self, env, prl_options, splitter, merger, **kwargs):
        # type: (Environment, ParallelizationOptions, Callable, Callable, Dict[str, Any]) -> None

        """Create a PBS instance that can run a function on
        :param env: Environment
        :param prl_options: Parallelization options that contains PBS configuration
        :param splitter: How to split input data
        :param merger: How to merge back all PBS job outputs
        :param kwargs: Other arguments
        """

        self._rng = np.random.RandomState(int(time.time()))

        self._dry_run = get_value(kwargs, "dry_run", False)
        self._env = env

        if prl_options is None:
            raise ValueError("prl_options cannot be None")

        self._prl_options = copy.deepcopy(prl_options)
        self._splitter = splitter
        self._merger = merger

    def _setup_pbs_run(self):
        mkdir_p(self._prl_options["pbs-pd-head"])

    def run(self, data, func, func_kwargs, **kwargs):
        # type: (Dict[str, Any], Callable, Dict[str, Any], Dict[str, Any]) -> Any
        """
        Run function on data using PBS scheduler
        :param data: Dictionary containing data name and value
        :param func: function to call on data with arguments in func_kwargs
        :param func_kwargs: Additional arguments to be passed
        :param kwargs:
        :return: output of merger function
        """

        job_name = get_value(kwargs, "job_name", "JOBNAME")
        pf_input_package_template_formatted = os.path.join(
            os.path.abspath(self._prl_options["pbs-pd-head"]), "input_package_{}"
        )

        num_jobs = self._prl_options["pbs-jobs"]

        # 1) Parse all PBS arguments
        self._setup_pbs_run()

        # 2) Create input packages files, one for every PBS run
        list_pf_input_job = self.create_input_package_files(
            data, func, func_kwargs, num_jobs,
            pf_package_template_formatted=pf_input_package_template_formatted,
            **kwargs
        )

        num_jobs = len(list_pf_input_job)

        # 3) Run all
        list_pf_output_job_packages = self.execute_function_on_input_packages(
            pf_input_package_template_formatted,
            job_name=job_name, num_jobs=num_jobs
        )

        # 4) Merge end-results
        data_output = None
        if not self._dry_run:
            data_output = self.merge_output_package_files(list_pf_output_job_packages)

        # 5) Clean
        if self._prl_options.safe_get("pbs-clean"):
            remove_p(*list_pf_input_job)
            remove_p(*list_pf_output_job_packages)

        return data_output

    def create_input_package_files(self, data, func, func_kwargs, num_splits, **kwargs):
        # type: (Dict, Callable, Dict[str, Any], int, Dict[str, Any]) -> List[str]
        """
        Run a function on the data using PBS
        :param data: the entire data
        :type data: DataHandler.D
        :param func: the function to execute on the (split) data
        :type func: Callable
        :param func_kwargs: the remaining arguments (i.e. not data) to be passed to the function
        :type func_kwargs: Dict[str, Any]
        :param num_splits: number of job splits
        :type num_splits: int
        :param kwargs:
        :return: List of paths to input package files
        :rtype: List[str]
        """

        pd_work_pbs = self._prl_options["pbs-pd-head"]

        pf_package_template_formatted = get_value(
            kwargs, "pf_package_template_formatted", os.path.join(pd_work_pbs, "input_package_{}")
        )

        # Split data
        list_split_data = self._splitter(data, num_splits)

        # Write package to disk
        list_pf_data = self._package_and_save_list_data(list_split_data, func, func_kwargs,
                                                        pf_package_template_formatted)

        # return list of filenames
        return list_pf_data

    def execute_function_on_input_packages(self, pf_input_package_template_formatted, job_name, num_jobs):
        # type: (str, str, int) -> List[str]
        """
        Create PBS file for run and execute it, returning the paths to all the job output packages
        :param pf_input_package_template_formatted:
        :param job_name:
        :param num_jobs:
        :returns: list of paths to output file packages
        :rtype: str
        """

        pd_head = self._prl_options["pbs-pd-head"]

        pf_pbs = os.path.join(self._prl_options["pbs-pd-head"], "run.pbs")
        pf_input_package_template = pf_input_package_template_formatted.format("${PBS_ARRAYID}")

        # create pbs file
        pf_output_package_template = "{}_output".format(pf_input_package_template)
        self._create_pbs_file(job_name, num_jobs, pf_pbs, pf_input_package_template, pf_output_package_template)

        # run
        if not self._dry_run:
            array_job_name = PBS._qsub(pf_pbs)

            # wait for jobs to end
            self._wait_for_job_array(array_job_name, pd_head)

        # collect all output files
        list_pf_outputs = []
        for x in range(1, num_jobs + 1):
            if os.path.isfile(PBS.create_concrete_from_template(pf_output_package_template + ".pkl", x)):
                list_pf_outputs.append(PBS.create_concrete_from_template(pf_output_package_template, x))

        # write summary file
        pf_pbs_summary = os.path.join(self._prl_options["pbs-pd-head"], self._prl_options["pbs-fn-summary"])
        write_to_file("\n".join(list_pf_outputs), pf_pbs_summary)

        return list_pf_outputs

    @staticmethod
    def _qsub(pf_pbs):
        # type: (str) -> str
        return run_shell_cmd("qsub  -V " + pf_pbs, do_not_log=True).strip()

    def _read_data_from_output_packages(self, list_pf_output_packages):

        list_data = list()

        for pf_output_package in list_pf_output_packages:
            list_data.append(PBSJobPackage.load(pf_output_package)["data"])

        return list_data

    def merge_output_package_files(self, list_pf_output_packages):

        list_output_data = self._read_data_from_output_packages(list_pf_output_packages)

        # 4-a) Merge data while loading packages one by one
        data_output = self._merger(list_output_data)

        return data_output

    def _package_and_save_data(self, data, func, func_kwargs, pf_package):
        # type: (Dict[str, Any], Callable, Dict[str, Any], str) -> None

        complete_func_kwargs = FunctionArguments(**func_kwargs).get_arguments(data)

        PBSJobPackage.save(
            {
                "func": func,
                "func_kwargs": complete_func_kwargs
            },
            pf_package
        )

    def _package_and_save_list_data(self, list_data, func, func_kwargs, pf_package_template_formatted):
        # type: (List[Dict[str, Any]], Callable, Dict[str, Any], str) -> List[str]

        list_pf = list()
        file_number = 1

        for data in list_data:
            pf_save = pf_package_template_formatted.format(file_number)

            self._package_and_save_data(data, func, func_kwargs, pf_save)
            list_pf.append(pf_save)

            file_number += 1

        return list_pf

    def _create_pbs_file(self, jobname, num_jobs, pf_pbs, pf_input_package_template, pf_output_package_template):
        """
        Create PBS file for runnning all input jobs
        :param jobname: Name of job
        :param num_jobs:
        :param pf_pbs:
        :param pf_input_package_template:
        :return:
        """

        # create unique compute directory
        pd_compute = None  # run_shell_cmd("mktemp --tmpdir={}".format(self._prl_options["pbs-pd-root-compute"]))

        pbs_text = PBS._generate_pbs_header_array(num_jobs, jobname, self._prl_options, pd_compute=pd_compute)

        pbs_text += "\n{}\n".format(
            PBS._generate_call_command(self._env,
                                       pf_input_package_template,
                                       pf_output_package_template,
                                       self._prl_options,
                                       pd_compute=pd_compute
                                       )
        )

        # write to file
        write_to_file(pbs_text, pf_pbs)



    @staticmethod
    def _generate_call_command(env, pf_job_input, pf_job_output, prl_options, pd_compute):

        pd_compute = os.path.abspath(os.path.join(prl_options["pbs-pd-root-compute"], prl_options["pbs-dn-compute"]))
        pd_job_template = os.path.join(pd_compute, "job_${PBS_ARRAYID}")

        cmd = "{} --pf-job-input {} --pf-job-output {} --pd-work {} -l {}".format(
            "python {}".format(os.path.join(env["pd-code"], "python/driver", "run-pbs-job.py")),
            pf_job_input,
            pf_job_output,
            pd_job_template,
            log.level
        )

        return cmd

    @staticmethod
    def create_concrete_from_template(pf_template, file_number):
        """Create a concrete file name based on template and file number
        e.g. Calling the function with filename_${PBS_ARRAYID}.txt, 5 returns
        filename_5.txt

        :param pf_template: template of file
        :type pf_template: str
        :param file_number: the file's number
        :returns: a concrete filename
        """

        return pf_template.replace("${PBS_ARRAYID}", str(file_number))

    def _wait_for_job_array(self, array_jobname, pd_work):
        # type: (str, str) -> None

        import string

        def _create_dummy_pbs_file(pf_dummy, jobname_dummy, pd_work):
            # type: (str, str, str) -> None
            pbs_text = PBS.generate_pbs_header(jobname_dummy, pd_work, 1, 1, "00:00:01")
            write_to_file(pbs_text, pf_dummy)

        def _cmd_run_dummy_and_wait(pf_dummy, jobname_dummy, jobname_array):
            cmd = "qsub -W depend=afteranyarray:{} {} \n".format(
                jobname_array,
                pf_dummy
            )

            cmd += r'while [ $(qstat -a | grep " R\|Q\|H " | grep ' + jobname_dummy + \
                   r'  | wc -l) != 0 ]; do sleep 60 ; done'

            return cmd

        # generate a random filename for the dummy job
        fn_dummy = ''.join(self._rng.choice(list(string.ascii_lowercase)) for _ in range(10))
        pf_dummy = os.path.join(pd_work, fn_dummy)

        # create an dummy pbs job that waits for the array to finish
        _create_dummy_pbs_file(pf_dummy, fn_dummy, pd_work)

        # generate pbs command to wait for job-array to finish and then run this dummy job
        cmd = _cmd_run_dummy_and_wait(pf_dummy, fn_dummy, array_jobname)

        # run command that waits
        run_shell_cmd(cmd, do_not_log=True)

    @staticmethod
    def generate_pbs_header(job_name, working_dir=".", num_nodes=1, ppn=1, walltime="00:30:00"):
        pbs_text = ""

        pbs_text += "#PBS -N " + str(job_name) + "\n"
        pbs_text += "#PBS -o " + str(working_dir) + "\n"
        pbs_text += "#PBS -j oe" + "\n"
        pbs_text += "#PBS -l nodes=" + str(num_nodes) + ":ppn=" + str(ppn) + "\n"
        pbs_text += "#PBS -l walltime=" + str(walltime) + "\n"

        pbs_text += "#PBS -W umask=002" + "\n"

        pbs_text += "set PBS_O_WORKDIR = " + str(working_dir) + "\n"
        pbs_text += "cd $PBS_O_WORKDIR \n"

        pbs_text += "echo The working directory is `echo $PBS_O_WORKDIR`" + "\n"
        pbs_text += "echo This job runs on the following nodes:" + "\n"
        pbs_text += "echo `cat $PBS_NODEFILE`" + "\n"

        return pbs_text

    @staticmethod
    def _generate_pbs_header_array(num_jobs, job_name, prl_options, pd_compute):
        """

        :param num_jobs:
        :param job_name:
        :param prl_options:
        :type prl_options: ParallelizationOptions
        :return:
        """

        num_nodes = prl_options["pbs-nodes"]
        ppn = prl_options["pbs-ppn"]
        walltime = prl_options["pbs-walltime"]

        pd_compute = os.path.abspath(os.path.join(prl_options["pbs-pd-root-compute"], prl_options["pbs-dn-compute"]))

        pd_job_template = os.path.join(pd_compute, "job_${PBS_ARRAYID}")

        pd_pbs_logs = os.path.join(prl_options["pbs-pd-head"], "pbs_logs")
        mkdir_p(pd_pbs_logs)

        node_property = prl_options.safe_get("pbs-node-property")
        if node_property is not None:
            node_property = ":" + node_property
        else:
            node_property = ""

        pbs_text = ""

        pbs_text += "#PBS -N " + str(job_name) + "\n"
        pbs_text += "#PBS -o " + "{}/{}".format(pd_pbs_logs, "error_${PBS_ARRAYID}") + "\n"
        pbs_text += "#PBS -j oe" + "\n"
        pbs_text += "#PBS -l nodes=" + str(num_nodes) + ":ppn=" + str(ppn) + "{}\n".format(node_property)
        pbs_text += "#PBS -l walltime=" + str(walltime) + "\n"

        if prl_options:
            array_param = "1-{}".format(num_jobs)
            if prl_options["pbs-concurrent-nodes"]:
                total_concurrent_jobs = prl_options["pbs-concurrent-nodes"] * int(8 / ppn)
                array_param = "{}%{}".format(array_param, total_concurrent_jobs)

            pbs_text += "#PBS -t {}".format(array_param) + "\n"

        pbs_text += "#PBS -W umask=002" + "\n"

        pbs_text += "export PATH=\"/home/karl/anaconda/envs/biogem_sbsp/bin:$PATH\"\n"

        pbs_text += "mkdir -p {}".format(pd_job_template) + "\n"

        pbs_text += "PBS_O_WORKDIR=" + pd_job_template + "\n"
        pbs_text += "cd $PBS_O_WORKDIR \n"
        pbs_text += "sleep 60\n"

        pbs_text += "echo The working directory is `echo $PBS_O_WORKDIR`" + "\n"
        pbs_text += "echo This job runs on the following nodes:" + "\n"
        pbs_text += "echo `cat $PBS_NODEFILE`" + "\n"

        return pbs_text