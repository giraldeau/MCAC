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

"""Test the groupby functions."""

import dask.array as da
import dask.dataframe as dd
import numpy as np
import pandas as pd
import pytest
import xarray as xr

from pymcac.tools.core.dask_tools import not_aligned_rechunk
from pymcac.tools.core.dataframe import groupby_agg, groupby_apply, xarray_to_frame
from pymcac.tools.core.groupby import groupby_aggregate
from pymcac.tools.core.sorting import sortby

from .generator import generate_dummy_aggregates_data, generate_dummy_data
from .test_dask import raise_if_computed
from .test_data import check_data


def pd_cumsum(df):
    """Cumsum written for pandas."""
    res = df.copy()
    res["data2"] = res["data"].cumsum()
    return res


def pd_cumsum_inplace(df):
    """Inplace Cumsum written for pandas."""
    res = df.copy()
    res["data"] = res["data"].cumsum()
    return res


def pd_duplicate_new_frame(df):
    """Duplicate a column of a pandas dataframe."""
    res = df.copy()
    res["data2"] = res["data"]
    return res


def pd_duplicate_new_series(df):
    """Duplicate a column of a pandas dataframe into a new series."""
    return pd_duplicate_new_frame(df)["data2"]


def custom_mean_pd(df):
    """Mean written for pandas."""
    return df.mean()


def noop_pd(df):
    """Do nothing."""
    return df


custom_mean_dd = dd.Aggregation(
    name="custom_mean",
    chunk=lambda s: (s.count(), s.sum()),
    agg=lambda count, sum: (count.sum(), sum.sum()),
    finalize=lambda count, sum: sum / count,
)

# ##################################
# ####### using dataframes #########
# ##################################

# *********************************
# ********* aggregation ***********
# *********************************


@pytest.mark.parametrize("dask", [0, 5])
@pytest.mark.parametrize("full", [True, False])
@pytest.mark.parametrize("sort_info", [True, False])
def test_groupby_agg_no_change(dask, full, sort_info):
    """Check an aggregation on one element."""
    aggregates = generate_dummy_aggregates_data(
        nt=29, nagg=31, dask=dask, full=full, sort_info=sort_info
    )
    # print(aggregates.compute())

    test = groupby_agg(aggregates, ["Time", "Label"], [("data", "sum", "data")])
    check_data(test)

    sorted_data = sortby(test, ["Time", "Label"])
    # print(sorted_data)
    # print(sorted_data.compute())

    aggregates.attrs["sort"] = ["Time", "Label"]
    assert sorted_data.identical(aggregates)


@pytest.mark.parametrize("dask", [0, 5])
@pytest.mark.parametrize("full", [True, False])
def test_groupby_agg_new_var(dask, full):
    """Check an aggregation on a new var."""
    aggregates = generate_dummy_aggregates_data(
        nt=29, nagg=31, dask=dask, full=full, sort_info=True
    )

    test = groupby_agg(
        aggregates, ["Time", "Label"], [("data", "sum", "data"), ("data2", "mean", "data")]
    )
    check_data(test)
    assert np.allclose(test.data, test.data2)

    if full:
        aggregates = aggregates.data
    aggregates.attrs["sort"] = ["Time", "Label"]

    sorted_data = sortby(test.data2, ["Time", "Label"])
    sorted_data.name = "data"

    assert sorted_data.identical(aggregates)


@pytest.mark.parametrize("sort", [True, False])
@pytest.mark.parametrize("fn", ["mean", np.sum, "custom_mean"])
@pytest.mark.parametrize("dask", [0, 5])
@pytest.mark.parametrize("full", [True, False])
@pytest.mark.parametrize("nt", [1, 29])
@pytest.mark.parametrize("nagg", [1, 31])
def test_groupby_agg(sort, fn, dask, full, nt, nagg):
    """test_groupby_agg."""
    if nt == 1 and nagg == 1:
        return
    aggregates = generate_dummy_aggregates_data(nt=nt, nagg=nagg, dask=dask, full=full)
    if isinstance(aggregates, xr.DataArray):
        aggregates = aggregates.to_dataset()
    chunks = aggregates.chunks
    aggregates["Np"] = aggregates.data.dims, np.random.randint(1, 10, aggregates.data.size)
    aggregates = not_aligned_rechunk(aggregates, chunks=chunks)

    if fn == "custom_mean":
        if aggregates.data.chunks is None:
            fn = custom_mean_pd
        else:
            fn = custom_mean_dd

    test = groupby_agg(
        aggregates, "Np", [("new_var", fn, "data"), ("data", "sum", "data")], sort=sort
    )
    # check_data(test)

    test = test.to_dataframe()[["new_var", "data"]]

    if fn is custom_mean_dd:
        fn = custom_mean_pd

    ref = (
        aggregates[["data", "Np"]]
        .to_dataframe()
        .groupby(by="Np", sort=True)
        .agg(
            new_var=pd.NamedAgg(column="data", aggfunc=fn),
            data=pd.NamedAgg(column="data", aggfunc="sum"),
        )
    )

    sorted_test = test.sort_index()

    if sort:
        assert sorted_test.equals(test)

    assert np.allclose(sorted_test, ref)


@pytest.mark.parametrize("dask_index", [True, False])
def test_groupby_agg_with_index(dask_index):
    """test_groupby_agg_with_index."""
    aggregates = generate_dummy_aggregates_data(nt=29, nagg=31, dask=50)
    chunks = aggregates.chunks
    aggregates["Np"] = ("k",), np.random.randint(1, 10, aggregates.sizes["k"])
    aggregates = not_aligned_rechunk(aggregates, chunks=chunks)

    values = da.unique(aggregates["Np"].values)
    if not dask_index:
        values = values.compute()

    test = groupby_agg(
        aggregates, "Np", [("data", "sum", "data")], index_arrays=values, length=9
    ).compute()
    ref = groupby_agg(aggregates, "Np", [("data", "sum", "data")]).compute()
    # check_data(test)

    assert test.equals(ref)


def test_groupby_agg_no_compute():
    """test_groupby_agg_no_compute."""
    aggregates = generate_dummy_aggregates_data(nt=29, nagg=31, dask=5)
    chunks = aggregates.chunks
    aggregates["Np"] = ("k",), np.random.randint(1, 10, aggregates.sizes["k"])
    aggregates["trigger"] = ("k",), da.from_delayed(
        raise_if_computed(), dtype=int, shape=(aggregates.sizes["k"],)
    )
    aggregates = not_aligned_rechunk(aggregates, chunks=chunks)

    # groupby_agg(aggregates, "Np", [("trigger", "sum", "trigger")])
    # groupby_agg(aggregates, "trigger", [("data", "sum", "data")])
    groupby_agg(aggregates, "trigger", [("data", "sum", "data")], index_arrays=da.arange(10))
    groupby_agg(aggregates, "Np", [("trigger", "sum", "trigger")], index_arrays=da.arange(10))


# *********************************
# ************ apply **************
# *********************************


@pytest.mark.parametrize("dask", [0, 5])
@pytest.mark.parametrize("full", [True, False])
@pytest.mark.parametrize("sort_info", [True, False])
def test_groupby_apply_no_change(dask, full, sort_info):
    """test_groupby_apply_no_change."""
    aggregates = generate_dummy_aggregates_data(
        nt=29, nagg=31, dask=dask, full=full, sort_info=sort_info
    )
    # print(aggregates)

    test = groupby_apply(aggregates, ["Time", "Label"], np.sum, "data", {"data": float})
    check_data(test)
    # print(test)

    sorted_data = sortby(test, ["Time", "Label"])
    check_data(sorted_data)

    # if full:
    #     aggregates = aggregates.data
    aggregates.attrs["sort"] = ["Time", "Label"]
    assert sorted_data.identical(aggregates)


@pytest.mark.parametrize("dask", [0, 5])
@pytest.mark.parametrize("full", [True, False])
def test_groupby_apply_new_frame(dask, full):
    """test_groupby_apply_new_frame."""
    aggregates = generate_dummy_aggregates_data(
        nt=29, nagg=31, dask=dask, full=full, sort_info=True
    )

    test = groupby_apply(
        aggregates,
        ["Time", "Label"],
        pd_duplicate_new_frame,
        "data",
        {"data": float, "data2": float},
    )
    check_data(test)
    assert np.allclose(test.data, test.data2)

    if full:
        aggregates = aggregates.data
    aggregates.attrs["sort"] = ["Time", "Label"]

    sorted_data = sortby(test.data, ["Time", "Label"])

    assert sorted_data.identical(aggregates)


@pytest.mark.parametrize("dask", [0, 5])
@pytest.mark.parametrize("full", [True, False])
def test_groupby_apply_new_series(dask, full):
    """test_groupby_apply_new_series."""
    aggregates = generate_dummy_aggregates_data(
        nt=29, nagg=31, dask=dask, full=full, sort_info=True
    )

    test = groupby_apply(
        aggregates, ["Time", "Label"], pd_duplicate_new_series, "data", {"data2": float}
    )
    check_data(test)

    if full:
        test = test.data2
        aggregates = aggregates.data
    assert np.allclose(test, aggregates.data)

    aggregates.attrs["sort"] = ["Time", "Label"]

    sorted_data = sortby(test, ["Time", "Label"])
    sorted_data.name = "data"

    assert sorted_data.identical(aggregates)


@pytest.mark.parametrize("inplace", [True, False])
@pytest.mark.parametrize("dask", [0, 5])
@pytest.mark.parametrize("full", [True, False])
@pytest.mark.parametrize("nt", [1, 29])
@pytest.mark.parametrize("nagg", [1, 31])
def test_groupby_apply(inplace, dask, full, nt, nagg):
    """test_groupby_apply."""
    if nt == 1 and nagg == 1:
        return
    aggregates = generate_dummy_aggregates_data(nt=nt, nagg=nagg, dask=dask, full=full)
    if isinstance(aggregates, xr.DataArray):
        aggregates = aggregates.to_dataset()
    chunks = aggregates.chunks
    aggregates["Np"] = aggregates.data.dims, np.random.randint(1, 10, aggregates.data.size)
    aggregates = not_aligned_rechunk(aggregates, chunks=chunks)

    if inplace:
        fn = pd_cumsum_inplace
        meta = {"data": float}
    else:
        fn = pd_cumsum
        meta = {"data": float, "data2": float}

    cols = ["Np"] + list(meta.keys())

    print(f"{aggregates=}")

    res_ds = groupby_apply(aggregates, "Np", fn, "data", meta, sort=False)
    print(f"{res_ds=}")
    print(f"{sortby(res_ds, ['Np', 'data'])=}")

    res_df = xarray_to_frame(res_ds, multi=False).reset_index()[cols].sort_values(cols)
    print(f"{res_df=}")

    ref = (
        xarray_to_frame(aggregates, multi=False)
        .reset_index()
        .groupby(by="Np", sort=False)
        .apply(fn)
        .reset_index(drop=True)[cols]
        .sort_values(cols)
    )
    print(f"{ref=}")

    assert np.allclose(res_df.values, ref.values)


@pytest.mark.parametrize("inplace", [True, False])
def test_groupby_apply_no_compute(inplace):
    """test_groupby_apply_no_compute."""
    aggregates = generate_dummy_aggregates_data(nt=29, nagg=31, dask=5)
    chunks = aggregates.chunks
    aggregates["Np"] = ("k",), np.random.randint(1, 10, aggregates.sizes["k"])
    aggregates["trigger"] = ("k",), da.from_delayed(
        raise_if_computed(), dtype=int, shape=(aggregates.sizes["k"],)
    )
    aggregates = not_aligned_rechunk(aggregates, chunks=chunks)

    if inplace:
        fn = pd_cumsum_inplace
        meta = {"data": float}
    else:
        fn = pd_cumsum
        meta = {"data": float, "data2": float}

    groupby_apply(aggregates, "trigger", fn, "data", meta)
    groupby_apply(aggregates, "Np", fn, "trigger", meta)


# ##################################
# ########### manually  ############
# ##################################


@pytest.mark.parametrize("dask", [0, 5])
def test_groupby_aggregate_reduction(dask):
    """test_groupby_aggregate_reduction."""
    aggregates, spheres = generate_dummy_data(nt=29, nagg=31, nsph=37, dask=dask)

    test = groupby_aggregate(spheres, aggregates, np.mean, "data")  # .compute()
    check_data(test)
    # print(test)

    ref = groupby_agg(spheres, ["Time", "Label"], [("data", np.mean, "data")])
    sorted_ref = sortby(ref, ["Time", "Label"]).data
    # print(sorted_ref)

    assert np.allclose(test, sorted_ref)


@pytest.mark.skip(reason="NOT IMPLEMENTED")
@pytest.mark.parametrize("dask", [0, 5])
def test_groupby_aggregate_not_reduction(dask):
    """test_groupby_aggregate_not_reduction."""
    assert False


def test_groupby_aggregate_no_compute():
    """test_groupby_aggregate_no_compute."""
    aggregates, spheres = generate_dummy_data(nt=29, nagg=31, nsph=37, dask=5)
    chunks = spheres.chunks
    spheres["trigger"] = ("k",), da.from_delayed(
        raise_if_computed(), dtype=int, shape=(spheres.sizes["k"],)
    )
    spheres = not_aligned_rechunk(spheres, chunks=chunks)

    groupby_aggregate(spheres, aggregates, np.mean, "trigger")
