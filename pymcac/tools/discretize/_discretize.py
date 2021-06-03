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
This library discretize a list of sphere
"""

from typing import Optional, Tuple

import numpy as np
import pandas as pd
from scipy.special import erf


def discretize(
    spheres: pd.DataFrame,
    aggregate: Optional[pd.Series] = None,
    resolution: int = 128,
    alphagangue: float = 0.4,
) -> Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Discretize the aggregate in a grid

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
    """
    Discretize a list of spheres in a grid

    alphagangue is a parameter allowing some gangue around the aggregate
    """

    xbounds, ybounds, zbounds = surrounding_box(spheres)

    (x, y, z) = mkgrid(xbounds, ybounds, zbounds, resolution)

    # fill the domain : positive value if inside at least a sphere
    if alphagangue > 0:
        data = np.zeros_like(x)
        # noinspection PyTypeChecker
        for posx, posy, posz, radius in spheres:
            d = np.sqrt((posx - x) ** 2 + (posy - y) ** 2 + (posz - z) ** 2)
            # noinspection PyTypeChecker
            data += 0.5 * (1 + erf(-(d / radius - 1) / alphagangue))
        # centering on 0
        data -= 0.5
    else:
        data = -np.ones_like(x) * np.inf
        # noinspection PyTypeChecker
        for posx, posy, posz, radius in spheres:
            d = (posx - x) ** 2 + (posy - y) ** 2 + (posz - z) ** 2
            data = np.maximum(data, radius ** 2 - d)

    return data > 0, (x, y, z)


def surrounding_box(
    spheres: np.ndarray,
) -> Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
    """
    compute dimensions of the surrounding box
    """
    rmax = spheres[:, 3].max()
    xmin, xmax = spheres[:, 0].min() - rmax, spheres[:, 0].max() + rmax
    ymin, ymax = spheres[:, 1].min() - rmax, spheres[:, 1].max() + rmax
    zmin, zmax = spheres[:, 2].min() - rmax, spheres[:, 2].max() + rmax

    return (xmin, xmax), (ymin, ymax), (zmin, zmax)


def mkgrid(
    xbounds: Tuple[float, float],
    ybounds: Tuple[float, float],
    zbounds: Tuple[float, float],
    resolution: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    compute dimensions of the surrounding box
    """
    xmin, xmax = xbounds
    ymin, ymax = ybounds
    zmin, zmax = zbounds

    # create a grid of the requested resolution
    # noinspection Mypy
    x, y, z = np.mgrid[
        xmin : xmax : resolution * 1j,  # type:ignore
        ymin : ymax : resolution * 1j,  # type:ignore
        zmin : zmax : resolution * 1j,  # type:ignore
    ]

    return x, y, z


if __name__ == "__main__":
    from pathlib import Path

    from pymcac import MCAC

    # The folder with all .h5 and .xmf files
    data_dir = Path("python-analysis/output_dir/")

    # # Read all data
    # Spheres, Aggregates = MCAC(data_dir).read()

    # last_agg = Aggregates.iloc[-1]

    # discretize(Spheres, last_agg)
