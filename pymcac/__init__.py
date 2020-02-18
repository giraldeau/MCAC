#!/usr/bin/env python3
# coding=utf-8

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
