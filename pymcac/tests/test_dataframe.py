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
import pytest
import xarray as xr
from numpy.testing import assert_array_equal

from pymcac.tools.core.dataframe import ddframe_to_xarray
from pymcac.tools.core.dataframe import frame_to_xarray
from pymcac.tools.core.dataframe import xarray_to_ddframe
from pymcac.tools.core.dataframe import xarray_to_frame

from .generator import generate_dummy_aggregates_data, generate_dummy_spheres_data
from .test_dask import raise_if_computed


@pytest.mark.parametrize("nt", [1, 3])
@pytest.mark.parametrize("nobj", [1, 5])
@pytest.mark.parametrize("data_type", ["aggregates", "spheres"])
@pytest.mark.parametrize("dask", [0, 1, 2])
@pytest.mark.parametrize("full", [True, False])
def test_xarray_to_ddframe(nt, nobj, data_type, dask, full):
    if data_type == "aggregates":
        data = generate_dummy_aggregates_data(nt=nt, nagg=nobj, dask=dask, full=full)
    elif data_type == "spheres":
        data = generate_dummy_spheres_data(nt=nt, nsph=nobj, dask=dask, full=full)
    else:
        raise ValueError(f"data_type {data_type} not understood")

    ref = data.compute()
    df = xarray_to_ddframe(data).compute()
    res = xr.Dataset.from_dataframe(df)

    print(ref)
    print(df)
    print(res)

    if full:
        ref = ref.data
    for coord in ref.coords:
        assert_array_equal(ref[coord].values, res[coord])
    assert_array_equal(ref.values, res["data"])


def test_xarray_to_ddframe_no_compute():
    data = generate_dummy_aggregates_data()

    data["trigger"] = ("k",), da.from_delayed(
        raise_if_computed(), dtype=int, shape=(data.sizes["k"],)
    )
    xarray_to_ddframe(data)


@pytest.mark.parametrize("nt", [1, 3])
@pytest.mark.parametrize("nobj", [1, 5])
@pytest.mark.parametrize("data_type", ["aggregates", "spheres"])
@pytest.mark.parametrize("dask", [0, 1, 2])
@pytest.mark.parametrize("full", [True, False])
def test_xarray_to_frame(nt, nobj, data_type, dask, full):
    if data_type == "aggregates":
        data = generate_dummy_aggregates_data(
            nt=nt, nagg=nobj, sort_info=True, dask=dask, full=full
        )
    elif data_type == "spheres":
        data = generate_dummy_spheres_data(nt=nt, nsph=nobj, sort_info=True, dask=dask, full=full)
    else:
        raise ValueError(f"data_type {data_type} not understood")

    ref = data.compute()
    res = xarray_to_frame(data)
    for col in res.index.names:
        res["k" + col] = res.index.get_level_values(col)

    print(ref)
    print(res)

    if not full:
        ref = ref.to_dataset()
    for data in ref.data_vars:
        assert_array_equal(ref[data].values, res[data])
    for coord in ref.coords:
        if coord[0] == "k":
            assert_array_equal(ref[data].values, res[data])
        elif coord[0] == "n":
            pass
        else:
            assert_array_equal(ref[coord].values, res.index.unique(coord))


@pytest.mark.parametrize("nt", [1, 3])
@pytest.mark.parametrize("nobj", [1, 5])
@pytest.mark.parametrize("data_type", ["aggregates", "spheres"])
@pytest.mark.parametrize("dask", [0, 1, 2])
def test_ddframe_to_xarray(nt, nobj, data_type, dask):
    if data_type == "aggregates":
        data = generate_dummy_aggregates_data(nt=nt, nagg=nobj, dask=dask, full=False)
    elif data_type == "spheres":
        data = generate_dummy_spheres_data(nt=nt, nsph=nobj, dask=dask, full=False)
    else:
        raise ValueError(f"data_type {data_type} not understood")

    ref = data
    df = xarray_to_ddframe(data)
    res = ddframe_to_xarray(df)

    print(ref.compute())
    print(df.compute())
    print(res.compute())

    assert ref.identical(res)
    assert ref.compute().identical(res.compute())


@pytest.mark.parametrize("dask", [0, 1, 2])
def test_ddframe_to_xarray_no_compute(dask):
    data = generate_dummy_aggregates_data(dask=dask)

    df = xarray_to_ddframe(data)
    lengths = tuple(df.index.map_partitions(len, enforce_metadata=False).compute())
    chunks = data.data.chunks if dask else (-1,)

    data["trigger"] = ("k",), da.from_delayed(
        raise_if_computed(), dtype=int, shape=(data.sizes["k"],)
    ).rechunk(chunks)

    df = xarray_to_ddframe(data)
    ddframe_to_xarray(df, lengths=lengths)


@pytest.mark.parametrize("nt", [1, 3])
@pytest.mark.parametrize("nobj", [1, 5])
@pytest.mark.parametrize("data_type", ["aggregates", "spheres"])
@pytest.mark.parametrize("dask", [0, 1, 2])
@pytest.mark.parametrize("full", [True, False])
def test_frame_to_xarray(nt, nobj, data_type, dask, full):
    if data_type == "aggregates":
        data = generate_dummy_aggregates_data(
            nt=nt, nagg=nobj, sort_info=True, dask=dask, full=full
        )
    elif data_type == "spheres":
        data = generate_dummy_spheres_data(nt=nt, nsph=nobj, sort_info=True, dask=dask, full=full)
    else:
        raise ValueError(f"data_type {data_type} not understood")

    ref = data.compute()
    df = xarray_to_frame(data)
    res = frame_to_xarray(df, full=full)

    print(ref)
    # print(df)
    print(res)

    assert ref.identical(res)
