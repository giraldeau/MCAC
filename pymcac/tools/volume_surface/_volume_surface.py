#!/usr/bin/env python3

# MCAC
# Copyright (C) 2020 CORIA
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Compute surface and volume of an aggregate."""
import multiprocessing
from typing import Optional, Union

import numpy as np
import pandas as pd
from skimage import measure

from pymcac.tools.discretize import discretize_spherelist


def volume_surface_sbl(
    spheres: Union[pd.DataFrame, np.ndarray], aggregate: Optional[pd.Series] = None
):
    """Compute surface and volume using SBL."""
    cols = ["Posx", "Posy", "Posz", "Radius"]
    if aggregate is not None:
        t, label = aggregate.name

        # get the spheres of the given aggregate
        spheres = spheres[spheres["Label"] == label].loc[t]

    if isinstance(spheres, pd.DataFrame):
        spheres = spheres[cols].values.copy()

    def thread_shield(spheres, volume, surface):
        _volume, _surface = compute_volume_surface_sbl(spheres)
        volume.value = _volume
        surface.value = _surface

    with multiprocessing.Manager() as manager:
        volume = manager.Value(float, 0.0)  # type: ignore
        surface = manager.Value(float, 0.0)  # type: ignore
        p = multiprocessing.Process(target=thread_shield, args=(spheres, volume, surface))
        p.start()
        p.join()
        if p.exitcode:
            print("retrying")
            # there is at least one contact point that make sbl crash
            spheres += 1e-10
            p = multiprocessing.Process(target=thread_shield, args=(spheres, volume, surface))
            p.start()
            p.join()
        if p.exitcode:
            raise RuntimeError("SBL is crashing")

        return volume.value, surface.value


def volume_surface_disc(
    spheres: Union[pd.DataFrame, np.ndarray],
    aggregate: Optional[pd.Series] = None,
    resolution: int = 128,
    alphagangue: float = 0,
):
    """Compute surface and volume using discretisation.

    alphagangue is a parameter allowing some gangue around the aggregate
    """
    cols = ["Posx", "Posy", "Posz", "Radius"]
    if aggregate is not None:
        t, label = aggregate.name

        # get the spheres of the given aggregate
        spheres = spheres[spheres["Label"] == label].loc[t]

    if isinstance(spheres, pd.DataFrame):
        spheres = spheres[cols].values.copy()

    data, (x, y, z) = discretize_spherelist(spheres, resolution=resolution, alphagangue=alphagangue)

    dx = abs(x[0, 0, 0] - x[1, 0, 0])
    dy = abs(y[0, 0, 0] - y[0, 1, 0])
    dz = abs(z[0, 0, 0] - z[0, 0, 1])

    dvol = dx * dy * dz

    volume = data.sum() * dvol

    # compute the 3d surface of the aggregate
    verts, faces, *_ = measure.marching_cubes(data, 0.0, spacing=(dx, dy, dz))

    surface = measure.mesh_surface_area(verts, faces)

    return volume, surface


try:
    from .sbl_wrapper import compute_volume_surface_sbl

    WithSBL = True
    volume_surface = volume_surface_sbl
except ImportError:
    print("SBL unavailable")
    WithSBL = False
    volume_surface = volume_surface_disc
