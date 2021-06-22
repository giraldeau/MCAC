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
Tools to mimic the pandas API
"""
from typing import Any, Callable, Dict, List, Tuple, Union

import dask.array as da
import numpy as np
import xarray as xr

from pymcac.tools.core.dask_tools import aligned_rechunk

from .sorting import sortby


def groupby_aggregate(
    Spheres: xr.Dataset,
    Aggregates: xr.Dataset,
    fn: Union[str, Callable],
    sph_names_in: Union[str, List[str]],
    agg_names_in: Union[str, List[str]] = None,
    meta_out: Dict[str, Any] = None,
    reduction=True,
    dtype=None,
) -> Union[xr.DataArray, xr.Dataset]:
    """TODO use xarray map block and remove intermediate dtype"""
    if not reduction:
        raise NotImplementedError
    if isinstance(sph_names_in, str):
        sph_names_in = [sph_names_in]
    if agg_names_in is None:
        agg_names_in = []
    elif isinstance(agg_names_in, str):
        agg_names_in = [agg_names_in]

    for i, name_in in enumerate(sph_names_in):
        if name_in in ("Time", "Num"):
            sph_names_in[i] = "k" + name_in
    for i, name_in in enumerate(agg_names_in):
        if name_in in ("Time", "Label"):
            agg_names_in[i] = "k" + name_in

    if set(agg_names_in) & set(agg_names_in):
        raise ValueError("You cannot have twice the same variable name")

    if meta_out is None:
        meta_out = {name_in: Spheres[name_in].dtype for name_in in sph_names_in}
    if dtype is None:
        dtypes = set(meta_out.values())
        if len(dtypes) == 1:
            dtype = dtypes.pop()
        else:
            dtype = np.float64

    Spheres = sortby(Spheres, ["Time", "Label"])
    Aggregates = sortby(Aggregates, ["Time", "Label"])

    if reduction:
        Aggregates = aligned_rechunk(Aggregates, "Time", k=None)
        Spheres = aligned_rechunk(Spheres, Time=Aggregates.chunks["Time"])
        chunks = Aggregates.chunks["k"]
    else:
        Spheres = aligned_rechunk(Spheres, "Time", k=None)
        Aggregates = aligned_rechunk(Aggregates, Time=Spheres.chunks["Time"])
        chunks = Spheres.chunks["k"]

    arrs_in = [Spheres[name_in].data for name_in in sph_names_in] + [
        Aggregates[name_in].data for name_in in agg_names_in
    ]

    argin_string = (
        "x" + ", x".join(map(str, range(len(arrs_in) + 1))) + ", fn, nout, intermediate_dtype"
    )
    argout_string = (
        "x0, (x"
        + ", x".join(map(str, range(1, len(arrs_in) + 1)))
        + ",), fn, nout, intermediate_dtype"
    )

    res = da.map_blocks(
        eval(f"lambda {argin_string}: groupby_aggregate_in_block({argout_string})"),
        Aggregates.Np.data,
        *arrs_in,
        fn=fn,
        nout=len(meta_out),
        intermediate_dtype=dtype,
        dtype=dtype,
        new_axis=[0],
        chunks=(len(meta_out), chunks),
    )

    res_ds: Union[xr.DataArray, xr.Dataset] = xr.Dataset()
    for i, (col, dtype) in enumerate(meta_out.items()):
        res_ds[col] = ("k",), res[i, :].astype(dtype)

    if reduction:
        res_ds = res_ds.assign_coords(Aggregates.coords)
    else:
        res_ds = res_ds.assign_coords(Spheres.coords)

    datavars = list(res_ds.data_vars.keys())
    if len(datavars) == 1:
        (data_var,) = datavars
        res_ds = res_ds[data_var]

    return res_ds


def groupby_aggregate_in_block(
    npagg: np.ndarray, arrs: Tuple[np.ndarray, ...], fn, nout, intermediate_dtype
) -> np.ndarray:
    res = np.empty((nout, npagg.size), dtype=intermediate_dtype)
    start = 0
    for iagg, nsph in enumerate(npagg):

        end = start + nsph

        subsets = (arr[start:end] for iarr, arr in enumerate(arrs))
        res[:, iagg] = fn(*subsets)

        start = end

    return res
