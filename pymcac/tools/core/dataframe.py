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
from typing import Any
from typing import Callable
from typing import Dict
from typing import List
from typing import Tuple
from typing import Union

import dask.array as da
import dask.dataframe as dd
import numpy as np
import pandas as pd
import xarray as xr

from .dask_tools import progress_compute


def xarray_to_ddframe(ds: xr.Dataset,
                      set_index: Union[bool, str] = None,
                      dim_order: List[str] = None) -> dd.DataFrame:
    """
    Convert a xarray.Dataset into a dask.dataframe.DataFrame.

    The dimensions, coordinates and data variables in this dataset form
    the columns of the DataFrame.

    The difference with ds.to_dask_dataframe is that
        * it does not produce a column for a dim with no coord.

    Parameters
    ----------
    ds: xr.Dataset
    dim_order : list, optional
        Hierarchical dimension order for the resulting dataframe. All
        arrays are transposed to this order and then written out as flat
        vectors in contiguous order, so the last dimension in this list
        will be contiguous in the resulting DataFrame. This has a major
        influence on which operations are efficient on the resulting dask
        dataframe.

        If provided, must include all dimensions on this dataset. By
        default, dimensions are sorted alphabetically.
    set_index : str, optional
        If set_index=True, the dask DataFrame is indexed by this dataset's
        dimension coordinate.
        If set_index=str, the dask DataFrame is indexed by this coordinate.

    Returns
    -------
    dask.dataframe.DataFrame
    """
    if dim_order is None:
        dim_order = list(ds.dims)
    elif set(dim_order) != set(ds.dims):
        raise ValueError(
                "dim_order {} does not match the set of dimensions on this "
                "Dataset: {}".format(dim_order, list(ds.dims))
                )

    ordered_dims = {k: ds.dims[k] for k in dim_order}

    columns = list(ordered_dims)
    columns.extend(k for k in ds.coords if k not in ds.dims)
    columns.extend(ds.data_vars)

    series_list = []
    for name in columns:
        try:
            var = ds.variables[name]
        except KeyError:
            continue
            # # dimension without a matching coordinate
            # size = ds.dims[name]
            # data = da.arange(size, chunks=size, dtype=np.int64)
            # var = Variable((name,), data)

        # IndexVariable objects have a dummy .chunk() method
        if isinstance(var, xr.IndexVariable):
            var = var.to_base_variable()

        dask_array = var.set_dims(ordered_dims).chunk(ds.chunks).data
        series = dd.from_array(dask_array.reshape(-1), columns=[name])
        series_list.append(series)

    df = dd.concat(series_list, axis=1)

    if set_index:
        if len(dim_order) == 1:
            (dim,) = dim_order
            df = df.set_index(dim)
        else:
            # triggers an error about multi-indexes, even if only one
            # dimension is passed
            df = df.set_index(dim_order)
    return df


def xarray_to_frame(ds: xr.Dataset) -> pd.DataFrame:
    """
    Convert a xarray.Dataset into a pandas.DataFrame.

    Compute the dataset and transform is as a pandas.DataFrame:
     * the dimensions and coordinates forms a multi_index
     * the data variables forms the columns of the DataFrame.

    Returns
    -------
    pandas.DataFrame
    """
    ds = progress_compute(ds)

    index = {**{"Time": None}, **{col: ds[col].data for col in ds.coords}}
    multi_index = pd.MultiIndex.from_arrays(arrays=list(index.values()), names=list(index.keys()))
    df = pd.DataFrame(data={col: ds[col].data for col in ds.data_vars}, index=multi_index)

    return df


def groupby_agg(ds: xr.Dataset,
                by: Union[str, List[str]],
                agg: List[Tuple[str, Union[str, Callable], str]],
                sort: bool = True,
                lengths: Tuple[int] = None,
                set_index: bool = True,
                **kwargs: Any) -> xr.Dataset:
    """
    Mimic the groupby().agg mechanism of pandas

    Parameters
    ----------
    ds : dataset
    by : label or list of labels
        Used to determine the groups for the groupby.
    agg : list of tuple (name_out, fn, name_in)
        name_out/name_in being the name of the resulting/source column
        fn being the name of a function, or the function itself
    sort: bool, default True
        Sort group keys. Get better performance by turning this off.
    lengths: Tuple of int, optional
        size of the resulting chunks if known.
        If not given, will be computed.
    set_index: bool = True
        If grouping by one label, use the resulting index as data dimension
    **kwargs
        Keyword arguments to pass to agg.

    Returns
    -------
    xarray.Dataset
    """
    if isinstance(by, str):
        by = [by]
    split_out = kwargs.get("kwargs", 1)

    content = list(ds.data_vars.keys()) + list(ds.coords.keys())
    keep = [name_in for name_out, fn, name_in in agg] + by
    remove = [var for var in content if var not in keep]

    ds = ds.drop_vars(remove)

    df_gb = xarray_to_ddframe(ds).groupby(by=by, sort=sort)
    res = xr.Dataset()
    index = by[0] if set_index else "k"
    for name_out, fn, name_in in agg:
        serie = df_gb[name_in].agg(fn, **kwargs)
        if not res.coords:
            if len(by) == 1:
                (coord,) = by
                index = coord if set_index else "k"
                if split_out == 1:
                    res.coords[coord] = (index,), serie.index.to_dask_array(lengths)
                    if lengths is None:
                        (size,) = res.coords[coord].shape
                        if np.isnan(size):
                            lengths = tuple(serie.index.map_partitions(len, enforce_metadata=False).compute())
                            res.coords[coord] = (index,), serie.index.to_dask_array(lengths)
                        else:
                            lengths = (size,)
                else:
                    arr = serie.index.to_dask_array(lengths)
                    if lengths is None:
                        lengths = serie.index.map_partitions(len, enforce_metadata=False)
                        lengths, arr = da.compute(lengths, arr)
                        lengths = tuple(lengths)
                    res.coords[coord] = (index,), arr
            else:
                index = "k"
                df_coords = serie.index.to_frame(lengths)
                if lengths is None:
                    lengths = tuple(serie.index.map_partitions(len, enforce_metadata=False).compute())
                for coord in by:
                    res.coords[coord] = (index,), df_coords[coord].to_dask_array(lengths)
        res[name_out] = (index,), serie.to_dask_array(lengths)

    if len(agg) == 1:
        (data_var,) = list(res.data_vars.keys())
        res = res[data_var]

    return res


def groupby_apply(ds: xr.Dataset,
                  by: Union[str, List[str]],
                  fn: Union[str, Callable],
                  name_in: Union[str, List[str]],
                  meta_out: Dict[str, Any] = None,
                  sort: bool = False,
                  lengths=None,
                  **kwargs) -> xr.Dataset:
    """
    Mimic the groupby().apply mechanism of pandas

    Parameters
    ----------
    ds : dataset
    by : label or list of labels
        Used to determine the groups for the groupby.
    name_in : label or list of labels
        Name(s) of the source column(s) passed to fn.
    fn : str or callable
    meta_out : Dict[name, dtype], optionnal
        Description of the new columns added by fn (for dask).
    sort: bool, default False
        Sort group keys. Get better performance by turning this off.
    lengths: Tuple of int, optional
        size of the resulting chunks if known.
        If not given, will be computed.
    **kwargs
        Keyword arguments to pass to apply.

    Returns
    -------
    xarray.Dataset
    """

    if isinstance(by, str):
        by = [by]

    if isinstance(name_in, str):
        name_in = [name_in]

    coords = list(ds.coords.keys())
    # noinspection PyTypeChecker
    keep = list(set(name_in + by + coords))
    remove = [var for var in ds.data_vars.keys() if var not in keep]

    ds = ds.drop_vars(remove)
    with_dask = ([ds[var].chunks is not None for var in ds.data_vars.keys()]
                 + [ds.coords[coord].chunks is not None for coord in ds.coords.keys()])

    if any(with_dask):
        df = xarray_to_ddframe(ds)

        res_df = df.groupby(by=by, sort=sort).apply(fn, meta={**df.dtypes.to_dict(), **meta_out}, **kwargs)

        res = xr.Dataset()

        if lengths is None:
            lengths = tuple(res_df.index.map_partitions(len, enforce_metadata=False).compute())

        for coord in coords:
            res.coords[coord] = ("k",), res_df[coord].to_dask_array(lengths)
        for name_out in meta_out:
            res[name_out] = ("k",), res_df[name_out].to_dask_array(lengths)
    else:
        df = ds.to_dataframe()
        res_df = df.groupby(by=by, **kwargs).apply(fn, **kwargs).reset_index()
        idx = res_df.index
        idx.name = "k"
        res_df.index = idx
        res = xr.Dataset(res_df)
    return res
