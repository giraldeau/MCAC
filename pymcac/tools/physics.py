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

import dask.array as da
import numpy as np
import xarray as xr
from numba import njit
from scipy.optimize import root_scalar


def mobility_diameter(f_agg: xr.DataArray, **kwargs) -> xr.DataArray:
    """
    Compute the mobility diameter based on f_agg
    """
    required = {"A1", "A2", "A3", "lambda_g", "mu_g"}
    if required - set(kwargs.keys()):
        raise ValueError(f"You must provide the following physical parameters: {required}")

    if isinstance(f_agg.data, da.Array):
        return xr.apply_ufunc(_vectorized_compute_dm,
                              f_agg, kwargs=kwargs,
                              dask="parallelized", output_dtypes=[float])
    else:
        return xr.apply_ufunc(_vectorized_compute_dm,
                              f_agg, kwargs=kwargs)


def _vectorized_compute_dm(vect_f_agg: np.ndarray, x0: float = 1e-8, **kwargs) -> np.ndarray:
    """
    Compute the mobility diameter based on f_agg
    """
    res = np.empty_like(vect_f_agg)
    for i, f_agg in enumerate(vect_f_agg):
        res[i] = x0 = _compute_dm(f_agg, x0=x0, **kwargs)
    return res


def _compute_dm(f_agg: float, x0: float = 1e-8, **kwargs) -> float:
    """
    Compute the mobility diameter based on f_agg
    """
    dm_params = ["A1", "A2", "A3", "lambda_g", "mu_g"]
    params = [kwargs[k] for k in dm_params]
    kwargs = {k: kwargs[k] for k in kwargs if k not in dm_params}

    dm_norm = x0
    x0 /= dm_norm
    rootresults = root_scalar(_normalized_dm_equation, args=(f_agg, dm_norm, *params),
                              fprime=True, fprime2=True,
                              x0=x0, x1=10 * x0, **kwargs)
    if not rootresults.converged:
        # raise ValueError("Warning, compute_Dm failed to find a solution ")
        print("Warning, compute_dm failed to find a solution ")
        return np.nan
    return rootresults.root * dm_norm


@njit(nogil=True, cache=True)
def _normalized_dm_equation(Dm: float = 1.0,
                            f_agg: float = 1e-12,
                            dm_norm: float = 1e-8,
                            A1: float = 1.142,
                            A2: float = 0.558,
                            A3: float = 0.999,
                            lambda_g: float = 5e-7,
                            mu_g: float = 6e-5):
    """
    Normalized version of the mobility diameter equation

    The default values are purely indicative
    """
    f_norm = 1 / f_agg / 1_000

    f, df, ddf = _dm_equation(Dm * dm_norm, f_agg, A1, A2, A3, lambda_g, mu_g)

    return (
        f / f_norm,
        df / f_norm * dm_norm,
        ddf / f_norm * dm_norm * dm_norm
        )


@njit(nogil=True, cache=True)
def _dm_equation(dm: float = 1e-8,
                 f_agg: float = 1e-12,
                 A1: float = 1.142,
                 A2: float = 0.558,
                 A3: float = 0.999,
                 lambda_g: float = 5e-7,
                 mu_g: float = 6e-5):
    """
    Mobility diameter equation

    The default values are purely indicative
    """
    kn = 2 * lambda_g / dm

    cunningham = 1 + A1 * kn + A2 * kn * np.exp(-A3 / kn)
    cunningham_prime = A1 * kn + A2 * (A3 + kn) * np.exp(-A3 / kn)
    cunningham_prime2 = 2 * A1 * kn + A2 * (A3 ** 2 / kn + 2 * A3 + 2 * kn) * np.exp(-A3 / kn)

    return (
        3 * np.pi * mu_g * dm - f_agg * cunningham,
        3 * np.pi * mu_g + f_agg * cunningham_prime / dm,
        - f_agg * cunningham_prime2 / dm ** 2
        )
