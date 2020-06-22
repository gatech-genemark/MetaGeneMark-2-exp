import logging
import threading
import time
from typing import *

from mg_general.general import get_value

log = logging.getLogger(__name__)


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

        time.sleep(30)


class GenericThread(threading.Thread):
    def __init__(self, func, func_kwargs, **kwargs):
        # type: (Callable, Dict[str, Any], Dict[str, Any]) -> None
        threading.Thread.__init__(self)
        self._func = func
        self._func_kwargs = func_kwargs

    def run(self):
        self._func(**self._func_kwargs)


def run_one_per_thread(data, func, data_arg_name, func_kwargs, **kwargs):
    # type: (Iterable[Any], Callable, str, Dict[str, Any], Dict[str, Any]) -> Any

    simultaneous_runs = get_value(kwargs, "simultaneous_runs", 8)

    active_threads = list()
    thread_id = 0
    for dp in data:

        # Create a thread for genome and run
        thread = GenericThread(func, {data_arg_name: dp, **func_kwargs})
        thread.start()
        thread_id += 1

        # wait until number of active threads is low
        if len(active_threads) >= simultaneous_runs:
            wait_for_any(active_threads)

        time.sleep(5)

    wait_for_all(active_threads)
