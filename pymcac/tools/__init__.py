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
Set of post-processing tools for MCAC
"""

from .autocorrelation import autocorrelation3d
from .coverage import coverages
from .discretize import discretize
from .inertia import inertia
from .volume_surface import volume_surface
from .core import DaskDistribute
from .core import progress_compute
from .core import groupby_agg
from .core import groupby_apply
from .core import xarray_to_ddframe
from .core import xarray_to_frame


__all__ = ["autocorrelation3d", "coverages", "discretize", "inertia", "volume_surface",
           "progress_compute", "DaskDistribute",
           "groupby_agg", "groupby_apply",
           "xarray_to_frame", "xarray_to_ddframe"]
