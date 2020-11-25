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

import psutil

from .plot import view_agg
from .plot import view_proj_agg
from .reader import MCAC
from .tools import autocorrelation3d
from .tools import coverages
from .tools import discretize
from .tools import inertia
from .tools import volume_surface
from .tools import dask_distribute, JupyterDaskDistribute
from .tools import progress_compute
from .tools import groupby_agg
from .tools import xarray_to_ddframe
from .tools import xarray_to_frame
from .tools import mobility_diameter


from .writer import export_ddscat
try:
    # noinspection PyUnresolvedReferences
    from .version import version
except ModuleNotFoundError:
    version = "Unknown"

ncores = psutil.cpu_count(logical=False)
os.environ['OMP_NUM_THREADS'] = os.environ.get('OMP_NUM_THREADS', str(ncores))

__all__ = ["MCAC",
           "view_agg", "view_proj_agg",
           "autocorrelation3d", "coverages", "discretize", "inertia", "volume_surface",
           "export_ddscat",
           "progress_compute", "dask_distribute", "JupyterDaskDistribute",
           "groupby_agg",
           "xarray_to_frame", "xarray_to_ddframe",
           "mobility_diameter"]
