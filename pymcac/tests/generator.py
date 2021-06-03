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
#:
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import dask.array as da
import numpy as np
import xarray as xr
from numba import njit

from pymcac.tools.core.dask_tools import aligned_rechunk, not_aligned_rechunk


@njit
def is_sorted(a):
    """Check that a numpy array is sorted"""
    for left, right in zip(a[:-1], a[1:]):
        if left > right:
            return False
    return True


def generate_dummy_aggregates_data(
    tstart=5.0, tend=10.0, nt=3, nagg=5, sort_info=False, dask=0, full=True
):
    """Generate dummy aggregate data with minimal data"""
    assert tstart <= tend
    assert nt > 0
    assert nagg > 0

    Time = np.linspace(tstart, tend, nt)

    Label = np.arange(nagg, dtype=np.int64)

    N = np.random.randint(1, nagg + 1, nt)
    N[np.random.randint(nt)] = nagg

    Nt = np.bincount(N)[:0:-1].cumsum()[::-1]

    kTime = np.concatenate([np.repeat(time, n) for time, n in zip(Time, N)])
    kLabel = np.concatenate([np.arange(n) for n in N])

    ds = xr.Dataset()
    ds.coords["Time"] = ("Time",), Time
    ds.coords["Label"] = ("Label",), Label
    ds.coords["kTime"] = ("k",), kTime
    ds.coords["kLabel"] = ("k",), kLabel
    ds.coords["nLabel"] = ("Time",), N
    ds.coords["nTime"] = ("Label",), Nt

    ds["data"] = ("k",), np.random.rand(ds.sizes["k"])

    if nagg == 1:
        ds = ds.drop_dims("Label").drop_vars({"kLabel", "nLabel"} & set(ds.coords))
    if nt == 1:
        ds = ds.drop_dims("Time").drop_vars({"kTime", "nTime"} & set(ds.coords))

    if nagg == 1 and nt > 1:
        ds["data"] = ds.data.swap_dims({"k": "Time"})
        ds = ds.drop_dims("k")
    elif nagg > 1 and nt == 1:
        ds["data"] = ds.data.swap_dims({"k": "Label"})
        ds = ds.drop_dims("k")
    elif nagg == 1 and nt == 1:
        ds = ds.isel(k=0)

    if dask:
        if nt > 1:
            ds = aligned_rechunk(ds, Time=dask)
        elif nagg > 1:
            ds = aligned_rechunk(ds, Label=dask)
        else:
            ds = not_aligned_rechunk(ds, {ds.data.dims: -1})

    if not full:
        ds = ds.data

    if sort_info:
        ds.attrs["sort"] = []
        if nt > 1:
            ds.attrs["sort"] += ["Time"]
        if nagg > 1:
            ds.attrs["sort"] += ["Label"]

    return ds


def generate_dummy_spheres_data(
    tstart=5.0, tend=10.0, nt=3, nsph=7, aggregates=None, sort_info=False, dask=0, full=True
):
    """Generate dummy spheres data with minimal data"""

    assert tstart <= tend
    assert nt > 0
    assert nsph > 0

    Time = np.linspace(tstart, tend, nt)
    Num = np.arange(nsph, dtype=np.int64)

    if aggregates is None:
        N = np.random.randint(1, nsph + 1, nt)
    else:
        aggregates = aggregates.compute()
        # We need at least one sphere per aggregate
        if "nLabel" in aggregates:
            N = np.concatenate(
                [np.random.randint(nagg_t, nsph + 1, 1) for nagg_t in aggregates.nLabel.data]
            )
        else:
            if "Label" in aggregates.dims:
                assert nt == 1
                nagg_t = aggregates.data.size
            else:
                nagg_t = 1
            N = np.random.randint(nagg_t, nsph + 1, nt)

    N[np.random.randint(nt)] = nsph

    Nt = np.bincount(N)[:0:-1].cumsum()[::-1]

    kTime = np.concatenate([np.repeat(time, n) for time, n in zip(Time, N)])
    kNum = np.concatenate([np.arange(n) for n in N])

    ds = xr.Dataset()
    ds.coords["Time"] = ("Time",), Time
    ds.coords["Num"] = ("Num",), Num
    ds.coords["kTime"] = ("k",), kTime
    ds.coords["kNum"] = ("k",), kNum
    ds.coords["nNum"] = ("Time",), N
    ds.coords["nTime"] = ("Num",), Nt

    ds["data"] = ("k",), np.random.rand(ds.sizes["k"])

    if aggregates:
        if "nLabel" in aggregates:
            Label = np.concatenate(
                [
                    np.random.permutation(
                        np.append(
                            np.arange(nagg_t, dtype=np.int64),
                            np.random.randint(0, nagg_t, nsph - nagg_t),
                        )
                    )
                    for nagg_t, nsph in zip(aggregates.nLabel.values, N)
                ]
            )
        else:
            if "Label" in aggregates.dims:
                assert nt == 1
                nagg_t = aggregates.data.size
            else:
                nagg_t = 1
            Label = np.concatenate(
                [
                    np.random.permutation(
                        np.append(
                            np.arange(nagg_t, dtype=np.int64),
                            np.random.randint(0, nagg_t, nsph - nagg_t),
                        )
                    )
                    for nsph in N
                ]
            )

        ds["Label"] = ("k",), Label

    if nsph == 1:
        ds = ds.drop_dims("Num").drop_vars({"kNum", "nNum"} & set(ds.coords))
    if nt == 1:
        ds = ds.drop_dims("Time").drop_vars({"kTime", "nTime"} & set(ds.coords))

    if nsph == 1 and nt > 1:
        ds["data"] = ds.data.swap_dims({"k": "Time"})
        if aggregates:
            ds["Label"] = ds.Label.swap_dims({"k": "Time"})
        ds = ds.drop_dims("k")
    elif nsph > 1 and nt == 1:
        ds["data"] = ds.data.swap_dims({"k": "Num"})
        if aggregates:
            ds["Label"] = ds.Label.swap_dims({"k": "Num"})
        ds = ds.drop_dims("k")
    elif nsph == 1 and nt == 1:
        ds = ds.isel(k=0)

    if dask:
        if nt > 1:
            ds = aligned_rechunk(ds, Time=dask)
        elif nsph > 1:
            ds = aligned_rechunk(ds, Num=dask)
        else:
            ds = not_aligned_rechunk(ds, k=-1)

    if not full:
        ds = ds.data

    if sort_info:
        ds.attrs["sort"] = []
        if nt > 1:
            ds.attrs["sort"] += ["Time"]
        if nsph > 1:
            ds.attrs["sort"] += ["Num"]

    return ds


def generate_dummy_data(tstart=5.0, tend=10.0, nt=3, nagg=5, nsph=7, sort_info=False, dask=0):
    """Generate dummy data with minimal data"""

    assert nagg <= nsph

    aggregates = generate_dummy_aggregates_data(
        tstart=tstart, tend=tend, nt=nt, nagg=nagg, sort_info=sort_info, dask=dask
    )
    spheres = generate_dummy_spheres_data(
        tstart=tstart,
        tend=tend,
        nt=nt,
        nsph=nsph,
        aggregates=aggregates,
        sort_info=sort_info,
        dask=dask,
    )
    Label = compute_Np_from_Label(spheres)
    if aggregates.data.chunks is not None:
        Label = da.from_array(Label, chunks=aggregates.data.chunks)

    if aggregates.data.dims:
        aggregates["Np"] = aggregates.data.dims, Label
    else:
        aggregates["Np"] = aggregates.data.dims, Label[0]

    return aggregates, spheres


def compute_Np_from_Label(spheres):
    spheres_df = spheres.Label.to_dataframe()
    if "kTime" in spheres:
        Np = spheres_df.groupby(["kTime", "Label"], sort=True).size()
    else:
        Np = spheres_df.groupby("Label", sort=True).size()
    return Np.values


def compute_Label_from_Np(aggregates):
    if not aggregates.Np.dims:
        aggregates = aggregates.expand_dims("None")
    aggregates_df = aggregates.Np.to_dataframe()

    if {"Time", "Label"}.issubset(set(aggregates.dims)):
        Label = aggregates_df.groupby(["kTime", "kLabel"], sort=True).apply(
            lambda df: np.repeat(df.kLabel, df.Np)
        )
        Label = Label.values
    elif "Time" in aggregates.dims:
        gb_dim = "kTime" if "k" in aggregates.dims else "Time"
        Label = aggregates_df.groupby([gb_dim], sort=True).apply(lambda df: np.repeat(0, df.Np))
        Label = np.concatenate(Label.values)
    elif "Label" in aggregates.dims:
        gb_dim = "kLabel" if "k" in aggregates.dims else "Label"
        Label = aggregates_df.groupby([gb_dim], sort=True).apply(
            lambda df: np.repeat(df.index, df.Np)
        )
        Label = np.concatenate(Label.values)
    else:
        Label = np.repeat(0, int(aggregates.Np))
    return Label
