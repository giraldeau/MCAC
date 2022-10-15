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

"""This library discretize a list of sphere."""

from typing import Optional, Tuple

import numpy as np
import pandas as pd
from numba import njit
from scipy.special import erf


def discretize(
    spheres: pd.DataFrame,
    aggregate: Optional[pd.Series] = None,
    resolution: int = 128,
    alphagangue: float = 0.4,
) -> Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Discretize the aggregate in a grid.

    alphagangue is a parameter allowing some gangue around the aggregate
    """
    cols = ["Posx", "Posy", "Posz", "Radius"]

    if aggregate is not None:
        t, label = aggregate.name

        # get the spheres of the given aggregate
        spheres = spheres[spheres["Label"] == label].loc[t]

    spherelist = spheres[cols].values.copy()

    return discretize_spherelist(spherelist, resolution, alphagangue)


def discretize_spherelist(
    spheres: np.ndarray, resolution: int = 128, alphagangue: float = 0.4
) -> Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Discretize a list of spheres in a grid.

    alphagangue is a parameter allowing some gangue around the aggregate
    """
    xbounds, ybounds, zbounds = surrounding_box(spheres)

    # fill the domain : positive value if inside at least a sphere
    if alphagangue > 0:
        x, y, z, data = _gangued_discretize(
            xbounds, ybounds, zbounds, spheres, resolution, alphagangue
        )
    else:
        x, y, z, data = _simple_discretize(xbounds, ybounds, zbounds, spheres, resolution)

    return (data > 0).astype(float), (x, y, z)


@njit(nogil=True, cache=True)
def _gangued_discretize(xbounds, ybounds, zbounds, spheres, resolution, alphagangue):
    xmin, xmax = xbounds
    ymin, ymax = ybounds
    zmin, zmax = zbounds

    x = np.linspace(xmin, xmax, resolution)
    y = np.linspace(ymin, ymax, resolution)
    z = np.linspace(zmin, zmax, resolution)

    data = np.zeros((resolution, resolution, resolution))
    for posx, posy, posz, radius in spheres:
        d = (
            (posx - x.reshape(-1, 1, 1)) ** 2
            + (posy - y.reshape(1, -1, 1)) ** 2
            + (posz - z.reshape(1, 1, -1)) ** 2
        )
        data += 0.5 * (1 + erf(-(d / radius - 1) / alphagangue))

    return x, y, z, data


@njit(nogil=True, cache=True)
def _simple_discretize(xbounds, ybounds, zbounds, spheres, resolution):
    xmin, xmax = xbounds
    ymin, ymax = ybounds
    zmin, zmax = zbounds

    resolution_x = resolution / (xmax - xmin)
    resolution_y = resolution / (ymax - ymin)
    resolution_z = resolution / (zmax - zmin)

    x = np.linspace(xmin, xmax, resolution)
    y = np.linspace(ymin, ymax, resolution)
    z = np.linspace(zmin, zmax, resolution)

    data = np.full((resolution, resolution, resolution), -np.inf)
    for posx, posy, posz, radius in spheres:

        imin = int(np.floor((posx - radius - xmin) * resolution_x))
        imax = int(np.ceil((posx + radius - xmin) * resolution_x)) + 1
        jmin = int(np.floor((posy - radius - ymin) * resolution_y))
        jmax = int(np.ceil((posy + radius - ymin) * resolution_y)) + 1
        kmin = int(np.floor((posz - radius - zmin) * resolution_z))
        kmax = int(np.ceil((posz + radius - zmin) * resolution_z)) + 1

        d = (
            (posx - x[imin:imax].reshape(-1, 1, 1)) ** 2
            + (posy - y[jmin:jmax].reshape(1, -1, 1)) ** 2
            + (posz - z[kmin:kmax].reshape(1, 1, -1)) ** 2
        )
        data[imin:imax, jmin:jmax, kmin:kmax] = np.maximum(
            data[imin:imax, jmin:jmax, kmin:kmax], radius**2 - d
        )

    return x, y, z, data


def surrounding_box(
    spheres: np.ndarray,
) -> Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
    """Compute dimensions of the surrounding box."""
    rmax = spheres[:, 3].max()
    xmin, xmax = spheres[:, 0].min() - rmax, spheres[:, 0].max() + rmax
    ymin, ymax = spheres[:, 1].min() - rmax, spheres[:, 1].max() + rmax
    zmin, zmax = spheres[:, 2].min() - rmax, spheres[:, 2].max() + rmax

    return (xmin, xmax), (ymin, ymax), (zmin, zmax)


if __name__ == "__main__":
    from pymcac import MCAC, validation_data_path

    simu = MCAC(validation_data_path / "pytest_data/")

    # Read all data
    Spheres, Aggregates = simu.spheres, simu.aggregates

    last_agg = Aggregates.iloc[-1]

    discretize(Spheres, last_agg)
