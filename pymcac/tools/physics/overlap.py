#!/usr/bin/env python3
# coding: utf-8

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
Function to compute usual physical quantities
"""
from typing import Tuple, Union

import numpy as np
import xarray as xr
from numba import njit

from pymcac.tools.core.groupby import groupby_aggregate


def overlapping(Spheres: xr.Dataset, Aggregates: xr.Dataset) -> Union[xr.DataArray, xr.Dataset]:
    return groupby_aggregate(
        Spheres,
        Aggregates,
        overlap_in_one_agg,
        ["Posx", "Posy", "Posz", "Radius"],
        None,
        {
            "cij_av": np.float64,
            "cij_median": np.float64,
            "Cij_v30": np.float64,
            "Cij_v20": np.float64,
            "tot_coll": np.int64,
            "n_c_avg": np.float64,
        },
    )


@njit(nogil=True, cache=True)
def overlap_in_one_agg(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, r: np.ndarray
) -> Tuple[np.float64, ...]:
    n = x.size

    tot_coll = np.float64(0.0)
    cij_av = cij_median = Cij_v30 = Cij_v20 = n_c_avg = np.float64(0.0)

    if n <= 1:
        return cij_av, cij_median, Cij_v30, Cij_v20, tot_coll, n_c_avg

    v = r ** 3

    # worst case scenario
    Cij = np.empty((n * (n - 1) // 2, 2))
    n_c_i = np.zeros(n)

    k = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            dij = np.sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2 + (z[i] - z[j]) ** 2)
            rij = r[i] + r[j]

            if 0.0 < dij <= rij * (1.0 + 1e-12):
                cij_k = (rij - dij) / rij
                vij_k = v[i] + v[j]
                Cij[k] = cij_k, vij_k
                k += 1
                n_c_i[i] += 2

    if not k:
        return cij_av, cij_median, Cij_v30, Cij_v20, tot_coll, n_c_avg

    tot_coll = np.float64(2.0 * k)

    cij, vij = Cij[:k].T

    cij_av = np.mean(cij)
    cij_median = np.median(cij)

    vij_sum: np.float64 = np.sum(vij)
    Cij_v20 = (np.sum(vij * cij ** 2) / vij_sum) ** (1 / 2)
    Cij_v30 = (np.sum(vij * cij ** 3) / vij_sum) ** (1 / 3)

    n_c_avg = np.mean(n_c_i)

    return cij_av, cij_median, Cij_v30, Cij_v20, tot_coll, n_c_avg


if __name__ == "__main__":
    import numpy as np

    from pymcac import MCAC, validation_data_path

    simu = MCAC(validation_data_path / "pytest_data/")

    # Read all data
    Spheres, Aggregates = simu.xspheres, simu.xaggregates
    overlap = overlapping(Spheres.compute(), Aggregates.compute()).compute()

    print(overlap)
