#!/usr/bin/env python3
# coding: utf-8

# MCAC
# Copyright (C) 2020 CORIA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Tools related to dask
"""
# import dask.array as da
# import cupy as cp
# import numpy as np
from dask import compute
from dask import persist
# noinspection PyProtectedMember
from dask.base import get_scheduler
from dask.diagnostics import ProgressBar
from dask.distributed import Client
from dask.distributed import progress
from psutil import cpu_count


def progress_compute(*dfs):
    """
    Compute the given dataframe, array, dataarray or datasets showing a progress bar
    """
    scheduler = get_scheduler()
    # distribute = "Client" in type(scheduler).__name__ and hasattr(scheduler, "get")
    distribute = "Client" in str(scheduler)

    if distribute:
        futures = persist(*dfs)
        progress(*futures, notebook=False, multi=False)
        res = compute(*futures)
        print()
    else:
        with ProgressBar():
            res = compute(*dfs)

    # for df in res:
    #     for coordname, coord in df.coords.items():
    #         if isinstance(coord.data, cp.ndarray):
    #             coord.data = coord.data.get()
    #     for varname, var in df.data_vars.items():
    #         if isinstance(var.data, cp.ndarray):
    #             var.data = var.data.get()

    if len(dfs) == 1:
        return res[0]

    return res


class DaskDistribute:
    """
    Context for dask distribute
    """
    __slots__ = ("n_workers", "threads_per_worker", "client")

    def __init__(self, n_workers=None, threads_per_worker=1):
        self.threads_per_worker = threads_per_worker
        self.n_workers = n_workers if n_workers is not None else cpu_count()
        self.client = None

    def start(self):
        """
        Initialize the workers
        """
        self.client = Client(n_workers=self.n_workers, threads_per_worker=self.threads_per_worker)
        return self.client

    def stop(self):
        """
        Shutdown the workers
        """
        self.client.close()

    def __enter__(self):
        return self.start()

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.stop()


# def _to_cpu(arr):
#     if isinstance(arr, np.ndarray):
#         return arr
#
#     assert isinstance(arr, da.Array)
#
#     if isinstance(arr._meta, np.ndarray):
#         return arr
#
#     if isinstance(arr._meta, cp.ndarray):
#         return arr.map_blocks(lambda x: x.get(), dtype=arr.dtype)
#     print(arr)
#     raise ValueError("What is this ?")
#
#
# def to_cpu(*dfs):
#     res = dfs[:]
#
#     for df in res:
#         for coordname, coord in df.coords.items():
#             coord.data = _to_cpu(coord.data)
#         for varname, var in df.data_vars.items():
#             var.data = _to_cpu(var.data)
#     return res
