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
import warnings
from typing import Any, Callable, Dict, Hashable, List, Mapping, Tuple, Union, cast

import dask.array as da
import dask.dataframe as dd
import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm

from .dask_tools import progress_compute
from .sorting import sortby


def xarray_to_ddframe(
    ds: xr.Dataset, set_index: Union[bool, str] = None, dim_order: List[str] = None
) -> dd.DataFrame:
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
        dim_order = [str(d) for d in ds.dims]
    elif set(dim_order) != set(ds.dims):
        raise ValueError(
            "dim_order {} does not match the set of dimensions on this "
            "Dataset: {}".format(dim_order, list(ds.dims))
        )

    ordered_dims = {k: ds.dims[k] for k in dim_order}

    columns = list(ordered_dims)
    columns.extend(str(k) for k in ds.coords if k not in ds.dims)
    columns.extend(str(d) for d in ds.data_vars)

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
            [dim] = dim_order
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

    index = {**{"Time": None}, **{str(col): ds[col].data for col in ds.coords}}
    multi_index = pd.MultiIndex.from_arrays(arrays=list(index.values()), names=list(index.keys()))
    df = pd.DataFrame(data={col: ds[col].data for col in ds.data_vars}, index=multi_index)

    return df


def try_extracting_index(index_arrays, length):
    known_indexes = True
    for i, index_array in enumerate(index_arrays):
        if isinstance(index_array, xr.DataArray):
            index_array = index_array.data
        if isinstance(index_array, (np.ndarray, da.Array)):
            index_length = index_array.size
            if np.isnan(index_length):
                index_array = index_array.compute_chunk_sizes()
                index_length = index_array.size
            if length is not None and not length == index_length:
                raise ValueError("Your indexes lengths are not coherent")
            length = index_length
            if isinstance(index_array, da.Array):
                index_array = index_array.rechunk(length)
        if index_array is None:
            known_indexes = False
        index_arrays[i] = index_array

    lengths = (length,) if length is not None else None
    return known_indexes, index_arrays, lengths


def groupby_agg(
    ds: Union[xr.DataArray, xr.Dataset],
    by: Union[str, List[str]],
    agg: List[Tuple[str, Union[str, Callable], str]],
    sort: bool = True,
    index_arrays=None,
    length: int = None,
) -> Union[xr.DataArray, xr.Dataset]:
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
    if isinstance(index_arrays, tuple):
        index_arrays = list(index_arrays)
    if not isinstance(index_arrays, list):
        index_arrays = [index_arrays for _ in by]

    if isinstance(ds, xr.DataArray):
        ds = ds.to_dataset(promote_attrs=True)

    known_indexes, index_arrays, lengths = try_extracting_index(index_arrays, length)

    relevent_dims = {str(dim) for var in ds.data_vars.values() for dim in var.dims}
    if len(relevent_dims) > 1:
        raise ValueError("You cannot groupby_agg with mixed shape Dataset")
    if not relevent_dims:
        raise ValueError("Nothing to group, your arrays are 0d")
    [k] = relevent_dims

    ds_only_k = ds.drop_dims([dim for dim in ds.dims if dim != k])
    ds_only_k = ds_only_k.rename(
        {coord: str(coord)[1:] for coord in ds_only_k.coords if str(coord).startswith("k")}
    )

    content = list(ds_only_k.data_vars.keys()) + list(ds_only_k.coords.keys())
    keep = [name_in for name_out, fn, name_in in agg] + by
    remove = [var for var in content if var not in keep]

    ds_only_k = ds_only_k.drop_vars(remove)
    dask = k in ds_only_k.chunks

    res: Union[xr.DataArray, xr.Dataset] = xr.Dataset()
    if dask:
        if sort and not known_indexes:
            index_serie = xarray_to_ddframe(ds_only_k[by]).groupby(by=by, sort=sort).first().index
            index_arrays, lengths = index_serie_to_index_array(index_serie, lengths)

        df_gb = xarray_to_ddframe(ds_only_k).groupby(by=by, sort=sort)
    else:
        df_gb = ds_only_k.to_dataframe().groupby(by=by, sort=sort)

    index = by[0]
    for name_out, fn, name_in in agg:
        serie = df_gb[name_in].agg(fn)

        if not res.coords:
            if not known_indexes:
                index_arrays, lengths = index_serie_to_index_array(serie.index, lengths)
            if len(by) == 1:
                [index] = by
                res.coords[index] = (index,), index_arrays[0]
            else:
                index = k
                for coord, index_array in zip(by, index_arrays):
                    res.coords[coord] = (index,), index_array

        if dask:
            res[name_out] = (index,), serie.to_dask_array(lengths)
        else:
            res[name_out] = (index,), serie.values

    for dim in ("Time", "Label", "Num"):
        if f"k{dim}" in res.coords:
            res = res.assign_coords({dim: ds[dim]})
            if ds[f"n{dim}"].dims[0] in res.dims:
                res = res.assign_coords({f"n{dim}": ds[f"n{dim}"].compute()})

    if "k" in res.dims:
        for coord_ in res.coords:
            coord = str(coord_)
            if f"k{coord}" in ds.coords:
                res = res.rename(cast(Mapping[Hashable, Hashable], {coord: f"k{coord}"}))

    datavars = res.data_vars.keys()
    if len(datavars) == 1:
        [datavar] = datavars
        res = res[datavar]

    return res


def index_serie_to_index_array(index_serie, lengths):
    index_frame = index_serie.to_frame()

    index_arrays = []
    for i, coord in enumerate(index_frame.columns):
        if isinstance(index_frame, dd.DataFrame):
            index_array = index_frame[coord].to_dask_array(lengths)
        else:
            index_array = index_frame[coord].values
        [size] = index_array.shape
        if np.isnan(size):
            index_serie.persist()
            index_array = index_serie.to_dask_array(lengths=True)
            [lengths] = index_array.chunks
        index_arrays += [index_array]
    return index_arrays, lengths


def groupby_apply(
    ds: xr.Dataset,
    by: Union[str, List[str]],
    fn: Union[str, Callable],
    name_in: Union[str, List[str]],
    meta_out: Dict[str, Any],
    sort: bool = False,
    lengths=None,
    **kwargs,
) -> Union[xr.DataArray, xr.Dataset]:
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

    ds = sortby(ds, by)

    if isinstance(name_in, str):
        name_in = [name_in]

    ds_only_k = ds.drop_dims([dim for dim in ds.dims if dim != "k"])
    ds_only_k = ds_only_k.rename(
        {coord: str(coord)[1:] for coord in ds_only_k.coords if str(coord).startswith("k")}
    )

    coords = [str(coord) for coord in ds_only_k.coords.keys()]
    # noinspection PyTypeChecker
    keep = list(set(name_in + by + coords))
    remove = [var for var in ds_only_k.data_vars.keys() if var not in keep]

    ds_only_k = ds_only_k.drop_vars(remove)
    with_dask = [ds_only_k[var].chunks is not None for var in ds_only_k.data_vars.keys()] + [
        ds_only_k.coords[coord].chunks is not None for coord in ds_only_k.coords.keys()
    ]
    res: Union[xr.DataArray, xr.Dataset] = xr.Dataset()
    if any(with_dask):
        df = xarray_to_ddframe(ds_only_k)

        res_df = df.groupby(by=by, sort=sort).apply(
            fn, meta={**df.dtypes.to_dict(), **meta_out}, **kwargs
        )

        if lengths is None:
            lengths = tuple(res_df.index.map_partitions(len, enforce_metadata=False).compute())
        else:
            res_df = res_df.repartition(len(lengths))

        for coord in coords:
            res.coords[coord] = ("k",), res_df[coord].to_dask_array(lengths)
        for name_out in meta_out:
            res[name_out] = ("k",), res_df[name_out].to_dask_array(lengths)
    else:
        df = ds_only_k.to_dataframe()

        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore", category=FutureWarning)
            tqdm.pandas()
        res_df = (
            df.groupby(by=by, sort=sort, **kwargs)
            .progress_apply(fn, **kwargs)
            .reset_index(drop=True)
        )

        for coord in coords:
            res.coords[coord] = ("k",), res_df[coord].to_numpy()
        for name_out in meta_out:
            res[name_out] = ("k",), res_df[name_out].to_numpy()

    res = res.rename({coord: "k" + str(coord) for coord in res.coords})
    res = res.assign_coords(ds.drop_dims("k").coords)

    # if sort:
    #     res = sortby(res, by=by)

    datavars = list(res.data_vars.keys())
    if len(datavars) == 1:
        [data_var] = datavars
        res = res[data_var]
    else:
        [idx_name] = ds.nTime.dims
        for var in ["nTime", f"n{idx_name}"]:
            res[var] = ds[var]
    return res
