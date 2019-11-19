#!/usr/bin/env python3
# coding=utf-8
"""
Set of Python libraries for post-processing MCAC data
"""

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


__all__ = ["MCAC",
           "view_agg", "view_proj_agg",
           "autocorrelation3d", "coverages", "discretize", "inertia", "volume_surface",
           "export_ddscat"]
