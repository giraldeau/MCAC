#!/usr/bin/env python3
# coding=utf-8
"""
Compute surface and volume of an aggregate
"""
from typing import Optional, Union

import numpy as np
import pandas as pd
from skimage import measure

from pymcac.tools.discretize import discretize


def volume_surface_sbl(spheres: Union[pd.DataFrame, np.ndarray],
                       aggregate: Optional[pd.Series] = None):
    """
        Compute surface and volume using SBL
    """
    cols = ["Posx", "Posy", "Posz", "Radius"]
    if aggregate is not None:
        t, label = aggregate.name

        # get the spheres of the given aggregate
        spheres = spheres[spheres["Label"] == label].loc[t]

    if isinstance(spheres, pd.DataFrame):
        spheres = spheres[cols].values.copy()

    return compute_volume_surface_sbl(spheres)


def volume_surface_disc(spheres: pd.DataFrame,
                        aggregate: Optional[pd.Series] = None,
                        resolution: int = 128,
                        alphagangue: float = 0):
    """
        Compute surface and volume using discretisation

        alphagangue is a parameter allowing some gangue around the aggregate
    """

    data, (x, y, z) = discretize(spheres, aggregate,
                                 resolution=resolution, alphagangue=alphagangue)

    dx = abs(x[0, 0, 0] - x[1, 0, 0])
    dy = abs(y[0, 0, 0] - y[0, 1, 0])
    dz = abs(z[0, 0, 0] - z[0, 0, 1])

    dvol = dx*dy*dy

    volume = data.sum() * dvol

    # compute the 3d surface of the aggregate
    verts, faces, *_ = measure.marching_cubes_lewiner(data, 0.,
                                                      spacing=(dx, dy, dz))

    surface = measure.mesh_surface_area(verts, faces)

    return volume, surface


try:
    from .sbl_wrapper import compute_volume_surface_sbl
    WithSBL = True
    volume_surface = volume_surface_sbl
except ImportError:
    WithSBL = False
    volume_surface = volume_surface_disc
