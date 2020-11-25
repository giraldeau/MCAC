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

from pymcac.tools.core.dask_tools import not_aligned_rechunk
from pymcac.tools.core.sorting import sortby
from .generator import generate_dummy_aggregates_data
from .test_dask import raise_if_computed
from .test_data import check_data


@pytest.mark.parametrize("dask", [0, 5])
@pytest.mark.parametrize("nagg", [1, 31])
@pytest.mark.parametrize("full", [True, False])
def test_sortby_nothing_to_do(dask, nagg, full):
    data = generate_dummy_aggregates_data(nt=29, nagg=nagg, dask=dask, sort_info=True, full=full)
    sorted_data = sortby(data, "Time")

    if full:
        for v, var in data.data_vars.items():
            assert var.data is sorted_data[v].data
    else:
        assert data.data is sorted_data.data

    for k, v in data.attrs.items():
        assert v is sorted_data.attrs[k]


@pytest.mark.parametrize("dask", [0, 5])
@pytest.mark.parametrize("nagg", [1, 31])
@pytest.mark.parametrize("full", [True, False])
def test_sortby_no_change(dask, nagg, full):
    data = generate_dummy_aggregates_data(nt=29, nagg=nagg, dask=dask, full=full)
    sorted_data = sortby(data, "Time")

    check_data(sorted_data)
    if full:
        if dask:
            assert "Time" in sorted_data.chunks
        else:
            assert "Time" not in sorted_data.chunks
        assert "Label" not in sorted_data.chunks
        if "k" in data.chunks:
            assert len(sorted_data.chunks["k"]) == len(data.chunks["k"])
    else:
        if sorted_data.chunks is not None:
            assert len(sorted_data.chunks) == len(data.chunks)

    data.attrs["sort"] = ["Time"]
    assert sorted_data.identical(data)

    if full:
        for v, var in data.data_vars.items():
            assert var.data is not sorted_data[v].data
    else:
        assert data.data is not sorted_data.data

    for k, v in data.attrs.items():
        assert v is not sorted_data.attrs[k]


@pytest.mark.parametrize("dask", [0, 5])
@pytest.mark.parametrize("nt", [1, 29])
@pytest.mark.parametrize("full", [True, False])
def test_sortby_Label(dask, nt, full):
    data = generate_dummy_aggregates_data(nt=nt, nagg=31, dask=dask, full=full)
    sorted_data = sortby(data, "Label")

    check_data(sorted_data)
    if full:
        assert "Time" not in sorted_data.chunks
        if dask:
            assert "Label" in sorted_data.chunks
        if "k" in data.chunks:
            assert len(sorted_data.chunks["k"]) == len(data.chunks["k"])
    else:
        if sorted_data.chunks is not None:
            assert len(sorted_data.chunks) == len(data.chunks)

    assert "Label" == sorted_data.attrs["sort"][0]

    if full:
        sorted_data = sorted_data.data
        data = data.data

    test = sorted_data.to_dataframe().reset_index(drop=True)
    if "k" in data.dims:
        ref = data.to_dataframe().sort_values(by="kLabel", kind="stable").reset_index(drop=True)
    else:
        ref = data.to_dataframe().sort_values(by="Label", kind="stable").reset_index(drop=True)
    assert test.equals(ref)


@pytest.mark.parametrize("dask", [0, 5])
def test_sortby_two(dask):
    data = generate_dummy_aggregates_data(nt=29, nagg=31, dask=dask)
    sorted_data = sortby(data, ["Time", "Label"])

    check_data(sorted_data)

    if dask:
        assert "Time" in sorted_data.chunks
    else:
        assert "Time" not in sorted_data.chunks
    assert "Label" not in sorted_data.chunks
    if "k" in data.chunks:
        assert len(sorted_data.chunks["k"]) == len(data.chunks["k"])

    data.attrs["sort"] = ["Time", "Label"]
    assert sorted_data.identical(data)


@pytest.mark.parametrize("dask", [0, 5])
def test_sortby_two_seq(dask):
    data = generate_dummy_aggregates_data(nt=29, nagg=31, dask=dask)

    Label_sorted_data = sortby(data, "Label")
    sorted_data = sortby(Label_sorted_data, "Time")

    check_data(sorted_data)
    if dask:
        assert "Time" in sorted_data.chunks
    else:
        assert "Time" not in sorted_data.chunks
    assert "Label" not in sorted_data.chunks
    if "k" in data.chunks:
        assert len(sorted_data.chunks["k"]) == len(data.chunks["k"])

    data.attrs["sort"] = ["Time", "Label"]
    assert sorted_data.identical(data)


@pytest.mark.parametrize("dask", [0, 5])
def test_sortby_other(dask):
    data = generate_dummy_aggregates_data(nt=29, nagg=31, dask=dask)

    sorted_data = sortby(data, "data")
    check_data(sorted_data)

    assert "Time" not in sorted_data.chunks
    assert "Label" not in sorted_data.chunks
    if "k" in data.chunks:
        assert len(sorted_data.chunks["k"]) == len(data.chunks["k"])

    test = sorted_data.data.to_dataframe().reset_index(drop=True)
    ref = data.data.to_dataframe().sort_values(by="data", kind="stable").reset_index(drop=True)
    assert test.equals(ref)
    assert "data" == sorted_data.attrs["sort"][0]


@pytest.mark.parametrize("dask", [0, 5])
def test_sortby_dim_and_other(dask):
    data = generate_dummy_aggregates_data(nt=29, nagg=31, dask=dask)
    sorted_data = sortby(data, ["Time", "data"])

    check_data(sorted_data)

    if dask:
        assert "Time" in sorted_data.chunks
    else:
        assert "Time" not in sorted_data.chunks
    assert "Label" not in sorted_data.chunks
    if "k" in data.chunks:
        assert len(sorted_data.chunks["k"]) == len(data.chunks["k"])

    test = sorted_data.data.to_dataframe().reset_index(drop=True)
    ref = data.data.to_dataframe().sort_values(by=["kTime", "data"], kind="stable").reset_index(drop=True)
    assert test.equals(ref)
    assert sorted_data.attrs["sort"] == ["Time", "data"]


@pytest.mark.parametrize("dask", [0, 5])
def test_sortby_other_and_back(dask):
    data = generate_dummy_aggregates_data(nt=29, nagg=31, dask=dask)

    sorted_data = sortby(data, "data")
    revert_data = sortby(sorted_data, ["Time", "Label"])

    data.attrs["sort"] = ["Time", "Label", "data"]

    assert revert_data.identical(data)


@pytest.mark.parametrize("dask", [0, 5])
def test_sortby_useless_sort(dask):
    data = generate_dummy_aggregates_data(nt=29, nagg=31, dask=dask)

    by = ["Time", "Label", "data"] + ["data", "Label", "Time"]

    sorted_data = sortby(data, by)
    assert sorted_data.attrs["sort"] == ["Time", "Label", "data"]


@pytest.mark.parametrize("nt", [1, 29])
def test_sortby_no_compute(nt):
    data = generate_dummy_aggregates_data(nt=nt, dask=5)
    chunks = data.chunks
    data["trigger"] = data.data.dims, da.from_delayed(raise_if_computed(),
                                                      dtype=int, shape=data.data.shape)
    data = not_aligned_rechunk(data, chunks=chunks)

    sortby(data, "data")
    sortby(data, "trigger")
    sortby(data, ["trigger", "data"])
    sortby(data, ["data", "trigger"])
