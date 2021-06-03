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
Plot a 3d view of the aggregate
"""
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# noinspection PyUnresolvedReferences
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from skimage import measure

from pymcac.tools.discretize import discretize


def view_agg(
    spheres: pd.DataFrame,
    aggregate: Optional[pd.Series] = None,
    resolution: int = 32,
    alphagangue: float = 0.4,
) -> None:
    """
    Plot a 3d view of the aggregate

    alphagangue is a parameter allowing some gangue around the aggregate
    """
    data, grid = discretize(spheres, aggregate, resolution, alphagangue)

    view_agg_disc(data, grid)


def view_agg_disc(data: np.ndarray, grid: Tuple[np.ndarray, np.ndarray, np.ndarray]) -> None:
    """
    Plot a 3d view of the discretized aggregate

    alphagangue is a parameter allowing some gangue around the aggregate
    """
    x, y, z = grid

    # extract space steps
    dx = abs(x[0, 0, 0] - x[1, 0, 0])
    dy = abs(y[0, 0, 0] - y[0, 1, 0])
    dz = abs(z[0, 0, 0] - z[0, 0, 1])

    # compute the 3d surface of the aggregate
    verts, faces, *_ = measure.marching_cubes_lewiner(data, 0.0, spacing=(dx, dy, dz))

    verts[:, 0] += x[0, 0, 0]
    verts[:, 1] += y[0, 0, 0]
    verts[:, 2] += z[0, 0, 0]

    # Plot the result
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_trisurf(verts[:, 0], verts[:, 1], faces, verts[:, 2])

    plt.show()


if __name__ == "__main__":
    from pathlib import Path

    from pymcac import MCAC

    # The folder with all .h5 and .xmf files
    data_dir = Path("python-analysis/output_dir/")

    # # Read all data
    # Spheres, Aggregates = MCAC(data_dir).read()

    # last_agg = Aggregates.iloc[-1]

    # view_agg(Spheres, last_agg)
