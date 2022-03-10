#!/usr/bin/env python3

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

"""Compute the overlapping coefficient COV."""

import os
from typing import Optional

import numpy as np
import pandas as pd

from .coverages_cython import coverages_cython, label_argsort


def coverages(
    spheres: pd.DataFrame, aggregates: pd.DataFrame, nprocs: Optional[int] = None
) -> np.ndarray:
    """Compute the overlapping coefficient COV.

    (which may happen with surface growth)
    """
    if nprocs is None:
        nprocs = len(os.sched_getaffinity(0))

    # Extracting
    data = spheres[["Label", "Posx", "Posy", "Posz", "Radius"]].copy()

    # Counting spheres per timestep
    nspheres = spheres["Posx"].groupby("Time").count().values
    nspheres = np.cumsum(nspheres)

    # Reset index
    data.reset_index(level="Num", drop=True, inplace=True)
    data.reset_index(inplace=True)

    # New index
    data.set_index(["Time", "Label"], inplace=True)

    # Sorting
    label_idx = data.index.get_level_values("Label").values
    idx = label_argsort(nspheres, label_idx, nprocs)

    # extracting
    X = data["Posx"].values[idx]
    Y = data["Posy"].values[idx]
    Z = data["Posz"].values[idx]
    R = data["Radius"].values[idx]

    # Get NP (already sorted)
    npp = aggregates["Np"].values

    # Compute
    return coverages_cython(npp, X, Y, Z, R, nprocs)


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np

    from pymcac import MCAC, validation_data_path

    simu = MCAC(validation_data_path / "pytest_data/")

    # Read all data
    Spheres, Aggregates = simu.spheres, simu.aggregates

    # Filter out aggregates that have only one spheres
    Aggregates = Aggregates[Aggregates["Np"] > 1].copy()

    # Compute the overlapping coefficient for all the aggregates
    Aggregates["cov"] = coverages(Spheres, Aggregates)

    # Compute the mean overlapping coefficient for each time step
    time_cov = Aggregates["cov"].groupby(by="Time").agg(np.mean)

    # Plot the evolution of the overlapping coefficient in time.
    fig, ax = plt.subplots()
    ax.loglog(time_cov)
    ax.set_xlabel("Time (s)", fontsize=9)
    ax.set_ylabel("Overlapping coefficient", fontsize=9)
    plt.suptitle("Time evolution of the overlapping coefficient", fontsize=11)
    plt.show()
