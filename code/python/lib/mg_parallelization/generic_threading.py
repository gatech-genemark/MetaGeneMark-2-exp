import logging
import math
import threading
import time
from typing import *

from mg_general.general import get_value

log = logging.getLogger(__name__)
logger=log


def wait_for_all(active_threads):
    # type: (List[threading.Thread]) -> List[threading.Thread]

    done_threads = list()
    while True:
        if len(active_threads) == 0:
            break
        for i, t in enumerate(active_threads):
            if not t.is_alive():
                del active_threads[i]
                done_threads.append(t)

        time.sleep(30)

    return done_threads


def wait_for_any(active_threads):
    # type: (List[threading.Thread]) -> Union[threading.Thread, None]

    while True:
        if len(active_threads) == 0:
            return None
        for i, t in enumerate(active_threads):
            if not t.is_alive():
                del active_threads[i]
                return t

        time.sleep(1)


class GenericThread(threading.Thread):
    def __init__(self, func, func_kwargs, **kwargs):
        # type: (Callable, Dict[str, Any], Dict[str, Any]) -> None
        threading.Thread.__init__(self)
        self._func = func
        self._func_kwargs = func_kwargs
        self._output = get_value(kwargs, "output", None, valid_type=dict())
        self._thread_id = get_value(kwargs, "thread_id", self.ident)


    def run(self):
        out = self._func(**self._func_kwargs)
        if self._output is not None:
            logger.debug(f"Self ident: {self._thread_id}")
            self._output[self._thread_id] = out


class GenericThreadN(threading.Thread):
    def __init__(self, func, list_func_kwargs, **kwargs):
        # type: (Callable, List[Dict[str, Any]], Dict[str, Any]) -> None
        threading.Thread.__init__(self)
        self._func = func
        self._list_func_kwargs = list_func_kwargs
        self._output = get_value(kwargs, "output", None, valid_type=dict())
        self._thread_id = get_value(kwargs, "thread_id", self.ident)

    def run(self):
        list_outputs = list()
        for func_kwargs in self._list_func_kwargs:
            output = self._func(**func_kwargs)

            if self._output is not None:
                list_outputs.append(output)

        if self._output is not None:
            self._output[self._thread_id] = list_outputs


def run_one_per_thread(data, func, data_arg_name, func_kwargs, **kwargs):
    # type: (Iterable[Any], Callable, str, Dict[str, Any], Dict[str, Any]) -> Any

    simultaneous_runs = get_value(kwargs, "simultaneous_runs", 8)

    active_threads = list()
    thread_id = 0
    output = dict()  # type: Dict[Any, List[Any]]


    for dp in data:

        # Create a thread for genome and run
        thread = GenericThread(func, {data_arg_name: dp, **func_kwargs},
                               output=output, thread_id=thread_id)
        thread.start()
        thread_id += 1
        active_threads.append(thread)
        logger.debug(f"Number of active threads: {len(active_threads)}")

        # wait until number of active threads is low
        if len(active_threads) >= simultaneous_runs:
            wait_for_any(active_threads)

        time.sleep(1)

    wait_for_all(active_threads)

    return [
        l for l in output.values()
    ]


def run_n_per_thread(data, func, data_arg_name, func_kwargs, **kwargs):
    # type: (Collection[Any], Callable, str, Dict[str, Any], Dict[str, Any]) -> Union[List[Any], None]
    simultaneous_runs = get_value(kwargs, "simultaneous_runs", 8)
    n = get_value(kwargs, "n", math.ceil(len(data) / simultaneous_runs), valid_type=int)
    arg_name_threadid = get_value(kwargs, "arg_name_threadid", None, type=str)

    output = dict()  # type: Dict[Any, List[Any]]

    data = list(data)

    if n * simultaneous_runs > len(data):
        n = math.ceil(len(data) / simultaneous_runs)

    active_threads = list()
    thread_id = 0

    thread_kwargs = dict()

    i = 0
    while i < len(data):
        if arg_name_threadid is not None:
            thread_kwargs[arg_name_threadid] = thread_id
            thread_id += 1

        # get n datapoints
        infos = list()
        counter = 0
        while i < len(data):
            infos.append({
                data_arg_name: data[i], **thread_kwargs, **func_kwargs
            })
            i += 1
            counter += 1
            if counter == n:
                break

        # Create a thread for genome and run
        thread = GenericThreadN(func, infos, thread_id=thread_id, output=output)
        thread.start()
        thread_id += 1

        active_threads.append(thread)

        # wait until number of active threads is low
        if len(active_threads) >= simultaneous_runs:
            wait_for_any(active_threads)

        # time.sleep(5)

    wait_for_all(active_threads)

    return [
        x for l in output.values() for x in l
    ]


def run_slice_per_thread(data, func, data_arg_name, func_kwargs, **kwargs):
    # type: (Collection[Any], Callable, str, Dict[str, Any], Dict[str, Any]) -> Union[List[Any], None]
    simultaneous_runs = get_value(kwargs, "simultaneous_runs", 8)
    n = get_value(kwargs, "n", math.ceil(len(data) / simultaneous_runs), valid_type=int)
    arg_name_threadid = get_value(kwargs, "arg_name_threadid", None, type=str)

    output = dict()  # type: Dict[Any, List[Any]]

    data = list(data)

    if n * simultaneous_runs > len(data):
        n = math.ceil(len(data) / simultaneous_runs)

    active_threads = list()
    thread_id = 0

    thread_kwargs = {}


    i = 0
    thread_id = 0
    while i < len(data):

        # get slice indices
        left = i
        right_excluded = min(i + n, len(data))

        data_slice = data[left:right_excluded]

        if arg_name_threadid is not None:
            thread_kwargs[arg_name_threadid] = thread_id
            thread_id += 1

        thread = GenericThread(func, {data_arg_name: data_slice, **thread_kwargs, **func_kwargs}, output=output)

        thread.start()
        thread_id += 1

        active_threads.append(thread)

        # wait until number of active threads is low
        if len(active_threads) >= simultaneous_runs:
            wait_for_any(active_threads)

        # time.sleep(5)

    wait_for_all(active_threads)

    return [
        x for l in output.values() for x in l
    ]
