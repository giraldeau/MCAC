#!/usr/bin/env python3
# coding=utf-8
"""
Set of Python libraries for post-processing MCAC data
"""
import os
from multiprocessing.pool import ThreadPool, Pool

import dask
from dask.cache import Cache
import psutil

from .plot import view_agg, view_proj_agg
from .reader import MCAC
from .tools import autocorrelation3d
from .tools import coverages
from .tools import discretize
from .tools import inertia
from .tools import volume_surface
from .writer import export_ddscat
try:
    # noinspection PyUnresolvedReferences
    from .version import version
except ModuleNotFoundError:
    version = "Unknown"

cache = Cache(2e9)  # Leverage two gigabytes of memory
cache.register()    # Turn cache on globally

ncores = psutil.cpu_count(logical=False)
dask.config.set(pool=ThreadPool(ncores))
# dask.config.set(pool=Pool(ncores))
# dask.config.set(scheduler='processes')
# dask.config.set(scheduler="synchronous")
os.environ['OMP_NUM_THREADS'] = os.environ.get('OMP_NUM_THREADS', str(ncores))


__all__ = ["MCAC",
           "view_agg", "view_proj_agg",
           "autocorrelation3d", "coverages", "discretize", "inertia", "volume_surface",
           "export_ddscat"]
