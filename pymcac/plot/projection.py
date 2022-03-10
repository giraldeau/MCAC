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

"""Plot a 2d projection of the aggregate."""

from typing import Callable, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pymcac.tools.discretize import discretize


def view_proj_agg(
    spheres: pd.DataFrame,
    aggregate: Optional[pd.Series] = None,
    resolution: int = 128,
    alphagangue: float = 0.4,
    projdir: str = "z",
    reduce: Callable = np.sum,
) -> None:
    """Plot a 2d projection of the aggregate in the direction x, y or z.

    if reduce is np.sum (or np.mean) the aggregate will be "transparent"
    if reduce is np.max  the aggregate will be fully opaque

    alphagangue is a parameter allowing some gangue around the aggregate
    """
    data, grid = discretize(spheres, aggregate, resolution, alphagangue)

    view_proj_disc(data, grid, projdir, reduce)


def view_proj_disc(
    data: np.ndarray,
    grid: Tuple[np.ndarray, np.ndarray, np.ndarray],
    projdir: str = "z",
    reduce: Callable = np.sum,
) -> None:
    """Plot a 2d projection of the discretized aggregate.

    if reduce is np.sum (or np.mean) the aggregate will be "transparent"
    if reduce is np.max  the aggregate will be fully opaque

    alphagangue is a parameter allowing some gangue around the aggregate
    """
    x, y, z = grid

    # binarization
    # noinspection PyTypeChecker,PyUnresolvedReferences
    data = (data > 0.0).astype(float)

    # start the figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # projection
    if projdir == "z":
        data = -reduce(data, axis=2)
        x = x[:, :, 0]
        y = y[:, :, 0]
        ax.set_xlabel("x")
        ax.set_ylabel("y")
    if projdir == "y":
        data = -reduce(data, axis=1)
        x = z[:, 0, :]
        y = x[:, 0, :]
        ax.set_xlabel("z")
        ax.set_ylabel("x")
    if projdir == "x":
        data = -reduce(data, axis=0)
        x = y[0, :, :]
        y = z[0, :, :]
        ax.set_xlabel("y")
        ax.set_ylabel("z")

    ax.contourf(x, y, data, cmap="gray")
    plt.show()


if __name__ == "__main__":
    from pymcac import MCAC, validation_data_path

    simu = MCAC(validation_data_path / "pytest_data/")

    # Read all data
    Spheres, Aggregates = simu.spheres, simu.aggregates

    last_agg = Aggregates.iloc[-1]

    view_proj_agg(Spheres, last_agg)
