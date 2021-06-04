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
from functools import wraps

import dask.array as da
import numpy as np
import pytest
import xarray as xr
from numba import njit

from pymcac.tools.core.various import get_idx_name

from .generator import compute_Label_from_Np
from .generator import compute_Np_from_Label
from .generator import generate_dummy_aggregates_data
from .generator import generate_dummy_data
from .generator import generate_dummy_spheres_data


@njit
def is_sorted(a):
    """Check that a numpy array is sorted"""
    for left, right in zip(a[:-1], a[1:]):
        if left > right:
            return False
    return True


def repeat(n_repetition):
    def decorator(function):
        @wraps(function)
        def wrapper(*args, **kwargs):
            for i in range(n_repetition):
                function(*args, **kwargs)

        return wrapper

    return decorator


def check_dims(ds):
    # assert "k" in ds.dims, "Missing the dimension k"
    # assert "Time" in ds.dims
    # assert ("Num" in ds.dims) or ("Label" in ds.dims)
    assert set(ds.dims).issubset({"k", "Label", "Num", "Time"}), "At least one dim is unknown"

    for coord in ds.coords.values():
        assert len(coord.dims) <= 1, "All coords should 1D at most"

    if isinstance(ds, xr.Dataset):
        for var in ds.data_vars.values():
            assert len(var.dims) <= 1, "All data_vars should be 1D at most"
    else:
        assert len(ds.dims) <= 1, "All data_vars should be 1D at most"


def check_vars(ds):
    """Check that all the mandatory variables are here"""
    coords = set(ds.coords)
    idx_name = get_idx_name(ds)
    assert idx_name in ("Num", "Label", None), f"{idx_name} is unknown"

    if "Time" in ds.dims:
        assert ds.Time.dims[0] == "Time", "Time dimension should be Time"
        if "k" in ds.dims:
            assert {"Time", "kTime"}.issubset(coords), "Missing Time and/or kTime"
            assert ds.kTime.dims[0] == "k", "kTime dimension should be k"
        if idx_name is not None:
            assert "nTime" in coords, "Missing nTime"
            assert ds.nTime.dims[0] == idx_name, f"nTime dimension should be {idx_name}"

    if idx_name is not None:
        assert ds[idx_name].dims[0] == idx_name, f"{idx_name} dimension should be {idx_name}"
        if "k" in ds.dims:
            assert {idx_name, f"k{idx_name}"}.issubset(coords), f"Missing {idx_name} or k{idx_name}"
            assert ds[f"k{idx_name}"].dims[0] == "k", f"k{idx_name} dimension should be k"
        if "Time" in ds.dims:
            assert f"n{idx_name}" in coords, f"Missing n{idx_name}"
            assert ds[f"n{idx_name}"].dims[0] == "Time", f"n{idx_name} dimension should be Time"

    # if "sort" in ds.attrs:
    #     for sort_idx in ds.attrs["sort"]:
    #         assert not sort_idx.startswith("k")
    #         assert sort_idx in ds.coords


def check_sorted(ds):
    """check that coords are sorted and that sort info is relevant"""

    def sorted_coords(ds):
        for c, coord in ds.coords.items():
            if c not in ("nTime", "nLabel", "nNum"):
                assert is_sorted(coord.values), f"{c} should be sorted"

    def sorted_kcoords(ds):
        idx = np.arange(ds.sizes["k"])
        for sort_idx in ds.attrs["sort"][::-1]:
            if f"k{sort_idx}" in ds.coords:
                data = ds[f"k{sort_idx}"].values[idx]
            else:
                data = ds[sort_idx].values[idx]
            idx = idx[data.argsort(kind="stable")]
        assert is_sorted(idx), "It seems that your data is not sorted"

    if "k" not in ds.dims:
        sorted_coords(ds)
        return

    if isinstance(ds, xr.Dataset):
        sorted_coords(ds.drop_dims(["k"]))

    if "sort" in ds.attrs:
        if set(ds.dims) - {"k"}:
            sorted_kcoords(ds.drop_dims([dim for dim in ds.dims if dim != "k"]))
        else:
            sorted_kcoords(ds)


def check_internal_kconsistency(ds):
    if "k" not in ds.dims:
        return

    ds = ds.compute()
    idx_name = get_idx_name(ds)

    if "Time" in ds.dims:
        assert set(ds.kTime.values) == set(
            ds.Time.values
        ), "Time and kTime contains differents values"

    if idx_name not in (None, "None"):
        assert set(ds[f"k{idx_name}"].values) == set(
            ds[idx_name].values
        ), f"{idx_name} and k{idx_name} contains differents values"

    if "Time" in ds.dims and idx_name not in (None, "None"):
        assert ds.nTime.sum() == ds.sizes["k"], "the sum of nTime should be the size of k"
        assert (
            ds[f"n{idx_name}"].sum() == ds.sizes["k"]
        ), f"the sum of n{idx_name} should be the size of k"
        assert np.array_equal(
            ds.nTime, np.bincount(ds[f"n{idx_name}"])[:0:-1].cumsum()[::-1]
        ), f"nTime and n{idx_name} are not coherent"

    elif "Time" in ds.dims:
        assert ds.sizes["Time"] == ds.sizes["k"], "The size of Time should be the size of k"
    elif idx_name not in (None, "None"):
        assert (
            ds.sizes[idx_name] == ds.sizes["k"]
        ), f"The size of {idx_name} should be the size of k"


def check_consistency(spheres, aggregates):
    spheres_time = spheres.coords.get("Time", None)
    aggregates_time = aggregates.coords.get("Time", None)
    if (spheres_time is not None) or (aggregates_time is not None):
        assert spheres_time is not None, "The Time should be present in spheres"
        assert aggregates_time is not None, "The Time should be present in aggregates"

        assert np.array_equal(
            spheres_time, aggregates_time
        ), "The Time should be the same in spheres and aggregate"
    Np = compute_Np_from_Label(spheres)
    if aggregates.Np.dims:
        assert np.array_equal(
            aggregates.Np, Np
        ), "Np computed from sphere's Label are not equal to the aggregate's Np"
    else:
        assert (
            aggregates.Np == Np[0]
        ), "Np computed from sphere's Label are not equal to the aggregate's Np"

    def sort_by_label(df):
        return df.sort_values("Label")

    df_unsorted_labels = spheres.Label.to_dataframe()
    if spheres_time is not None:
        df_sorted_labels = df_unsorted_labels.groupby("kTime", sort=True).apply(sort_by_label)
    else:
        df_sorted_labels = df_unsorted_labels.sort_values("Label")
    sorted_labels = df_sorted_labels.Label.values

    labels = compute_Label_from_Np(aggregates)
    assert np.array_equal(
        sorted_labels, labels
    ), "labels computed from aggregate's Np are not equal to the sphere's Labels"


def check_dask_consistency(ds):
    if isinstance(ds, xr.DataArray):
        chunk = ds.chunks if "k" in ds.dims else None
        for coord in ds.coords.values():
            assert coord.chunks == chunk, f"{coord.name} is not correctly chunked"
        return

    if "k" not in ds.dims:
        return

    idx_name = get_idx_name(ds)

    for dim in ("k", "Time", idx_name):
        chunk = ds.chunks.get(dim, None)
        if chunk is not None:
            chunk = (chunk,)

        if dim == "k":
            for coord in ds.coords.values():
                if dim in coord.dims:
                    assert coord.chunks == chunk, f"{coord.name} is not correctly chunked"
        elif chunk is not None:
            if "sort" in ds.attrs:
                assert ds.attrs["sort"][0] == dim, f"{dim} should not be chunked if not sorted"

            if dim == "Time":
                coord = da.from_array(ds[dim].data, chunks=ds[f"n{idx_name}"].chunks)
            else:
                coord = da.from_array(ds[dim].data, chunks=ds.nTime.chunks)

            aligned = da.map_blocks(
                lambda coord_block, kcoord_block: np.array([set(coord_block) == set(kcoord_block)]),
                coord,
                ds["k" + dim].data,
                dtype=bool,
            )
            assert np.all(aligned), f"The chunks are not aligned of {dim}"

        for var in ds.data_vars.values():
            if dim in var.dims:
                assert var.chunks == chunk, f"{var.name} is not correctly chunked"


def check_data(ds, aligned=True):
    assert isinstance(ds, (xr.Dataset, xr.DataArray))

    print(ds)

    # print("dims")
    check_dims(ds)
    # print("vars")
    check_vars(ds)
    # print("sorted")
    check_sorted(ds)

    if isinstance(ds, xr.Dataset):
        # print("internal")
        check_internal_kconsistency(ds)
    if aligned:
        # print("dask")
        check_dask_consistency(ds)


@pytest.mark.parametrize("nt", [1, 3])
@pytest.mark.parametrize("nobj", [1, 5])
@pytest.mark.parametrize("data_type", ["aggregates", "spheres"])
@pytest.mark.parametrize("dask", [0, 1, 2])
@pytest.mark.parametrize("full", [True, False])
@repeat(10)
def test_generate(nt, nobj, data_type, dask, full):
    if data_type == "aggregates":
        data = generate_dummy_aggregates_data(
            nt=nt, nagg=nobj, sort_info=True, dask=dask, full=full
        )
    elif data_type == "spheres":
        data = generate_dummy_spheres_data(nt=nt, nsph=nobj, sort_info=True, dask=dask, full=full)
    else:
        raise ValueError(f"data_type {data_type} not understood")

    check_data(data)


@pytest.mark.parametrize("nt", [1, 3])
@pytest.mark.parametrize("nagg", [1, 5])
@pytest.mark.parametrize("dask", [0, 1, 2])
@repeat(10)
def test_generate_both(nt, nagg, dask):
    aggregates, spheres = generate_dummy_data(nt=nt, nagg=nagg, sort_info=True, dask=dask)

    check_data(aggregates)
    check_data(spheres)
    check_consistency(spheres, aggregates)


# TODO hypothesis
