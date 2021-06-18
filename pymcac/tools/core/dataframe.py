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
from itertools import chain
from typing import Any, Callable, Dict, Hashable, List, Mapping, Tuple, Union, cast
from functools import wraps

import dask.array as da
import dask.dataframe as dd
import numpy as np
import pandas as pd
import xarray as xr
from tqdm import tqdm
from xarray.core.utils import K

from .dask_tools import progress_compute
from .sorting import sortby


def xarray_to_ddframe(
    ds: Union[xr.DataArray, xr.Dataset],
    set_index: Union[bool, str] = None,
    dim_order: List[Hashable] = None,
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
    if isinstance(ds, xr.DataArray):
        ds = ds.to_dataset()

    if len(ds.coords) > 1:
        no_k = [d for d in ds.dims if d != "k"]
        ds = ds.drop_dims(no_k)

    ordered_dims = ds._normalize_dim_order(dim_order=dim_order)

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
        dim_order = [*ordered_dims]

        if len(dim_order) == 1:
            (dim,) = dim_order
            df = df.set_index(dim)
        else:
            # triggers an error about multi-indexes, even if only one
            # dimension is passed
            df = df.set_index(dim_order)
    else:
        if len(ds.dims) == 1:
            index = df.index
            [index.name] = ds.dims
            df.index = index

    return df


def ddframe_to_xarray(
    df: dd.DataFrame, lengths: Tuple[int] = None
) -> Union[xr.DataArray, xr.Dataset]:
    """Convert a dask.dataframe.DataFrame into an xarray.Dataset"""

    if lengths is None:
        lengths = tuple(
            df.index.map_partitions(len, enforce_metadata=False).compute()
        )  # type:ignore
    lengths = cast(Tuple[int], lengths)

    arrays = [(k, v.to_dask_array(lengths)) for k, v in df.items()]
    index_name = df.index.name if df.index.name is not None else "k"

    res: Union[xr.DataArray, xr.Dataset] = xr.Dataset()
    for name, values in arrays:
        if name[0] == "k":
            res.coords[name] = ("k",), values
        else:
            res[name] = (index_name,), values

    if sum(lengths) == 1:
        res = res.isel(k=0)

    datavars = list(res.data_vars.keys())
    if len(datavars) == 1:
        [data_var] = datavars
        res = res[data_var]

    return res


def xarray_to_frame(ds: Union[xr.DataArray, xr.Dataset], multi=True) -> pd.DataFrame:
    """
    Convert a xarray.Dataset into a pandas.DataFrame.

    Compute the dataset and transform is as a pandas.DataFrame:
     * the dimensions and coordinates forms a multi_index
     * the data variables forms the columns of the DataFrame.

    Returns
    -------
    pandas.DataFrame
    """
    sort = ds.attrs.get("sort", [])

    if isinstance(ds, xr.DataArray):
        ds = ds.to_dataset()

    if len(ds.coords) > 1:
        no_k = [d for d in ds.dims if d != "k"]
        ds = ds.drop_dims(no_k)

    if not multi:
        return ds.to_dataframe()

    if (len(ds.coords) > 1) and len(ds.coords) != len(sort):
        raise ValueError(
            "You must sort the array in order to transform it into a multi_indexed DataFrame"
        )

    ds = progress_compute(ds)

    index = {str(col): ds[col].data for col in ds.coords}
    if len(ds.coords) > 1:
        index = {col: index["k" + col] for col in sort}

    if not index:
        index = {"k": [1]}

    multi_index = pd.MultiIndex.from_arrays(arrays=list(index.values()), names=list(index.keys()))
    df = pd.DataFrame(data={col: ds[col].data for col in ds.data_vars}, index=multi_index)

    return df


def frame_to_xarray(df: pd.DataFrame, full: bool = True) -> Union[xr.DataArray, xr.Dataset]:
    """Convert a pandas.DataFrame into an xarray.Dataset"""

    if df.index.is_monotonic_increasing:
        sort = [c for c in df.index.names if c != "k"]

    df_multi = None
    if df.index.nlevels > 1:
        df_multi = df

        for col in df.index.names:
            df["k" + col] = df.index.get_level_values(col)

        df = df.reset_index(drop=True)
        df.index.name = "k"

    res: Union[xr.DataArray, xr.Dataset] = xr.Dataset.from_dataframe(df)
    if "k" in res.coords:
        res = res.drop_vars("k")

    coords = [coord for coord in res.data_vars if coord[0] == "k"]  # type:ignore
    res = res.set_coords(coords)

    if len(df) == 1:
        res = res.isel(k=0)

    if full and df_multi is not None:
        if df.index.nlevels > 2:
            raise ValueError("You cannot reconstruct a full xarray with more that 2 levels")
        a, b = df_multi.index.names
        df_a = df_multi[df_multi.columns[0]].groupby(a).count()
        df_b = df_multi[df_multi.columns[0]].groupby(b).count()
        res.coords[a] = (a,), df_a.index.values
        res.coords["n" + a] = (b,), df_b.values
        res.coords[b] = (b,), df_b.index.values
        res.coords["n" + b] = (a,), df_a.values

    if not full:
        datavars = list(res.data_vars.keys())
        if len(datavars) == 1:
            [data_var] = datavars
            res = res[data_var]

    if sort is not None:
        res.attrs["sort"] = sort

    return res


def try_extracting_index(index_arrays, length):
    known_indexes = True
    for i, index_array in enumerate(index_arrays):
        if isinstance(index_array, xr.DataArray):
            index_array = index_array.data
        if isinstance(index_array, (np.ndarray, da.Array)):
            index_length = index_array.size
            if np.isnan(index_length):
                if length is not None:
                    # assuming indexes lengths are coherent in order to avoid computation
                    index_length = length
                else:
                    index_array = index_array.compute_chunk_sizes()
                    index_length = index_array.size
            if length is not None and not length == index_length:
                raise ValueError("Your indexes lengths are not coherent")
            length = index_length
            # if isinstance(index_array, da.Array):
            #     index_array = index_array.rechunk(length)
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
    Mimic the groupby().agg mechanism of pandas using dask dataframe

    The result have *one* chunk per variable (`split_out` is not supported)
    This will trigger some computation if `index_arrays` and `length` are not given

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
    length: int, optional
        size of the resulting arrays.
        If not given, will be computed.
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

    known_indexes, index_arrays, lengths = try_extracting_index(index_arrays, length)
    if dask:
        if sort and not known_indexes:
            index_serie = xarray_to_ddframe(ds_only_k[by]).groupby(by=by, sort=sort).first().index
            index_arrays, lengths = index_serie_to_index_array(index_serie, lengths)
            known_indexes = True

        df_gb = xarray_to_ddframe(ds_only_k).groupby(by=by, sort=sort)
    else:
        df_gb = xarray_to_frame(ds_only_k, multi=False).groupby(by=by, sort=sort)

    index = by[0] if len(by) == 1 else k

    res: Union[xr.DataArray, xr.Dataset] = xr.Dataset()
    for name_out, fn, name_in in agg:
        serie = df_gb[name_in].agg(fn)

        if not res.coords:
            if not known_indexes:
                index_arrays, lengths = index_serie_to_index_array(serie.index, lengths)
            if len(by) == 1:
                res.coords[index] = (index,), index_arrays[0]
            else:
                for coord, index_array in zip(by, index_arrays):
                    res.coords[coord] = (index,), index_array

        if dask:
            res[name_out] = (index,), serie.to_dask_array(lengths)
        else:
            res[name_out] = (index,), serie.values

    if "k" in res.dims:
        for coord_ in res.coords:
            coord = str(coord_)
            if f"k{coord}" in ds.coords:
                res = res.rename(cast(Mapping[Hashable, Hashable], {coord: f"k{coord}"}))

    for dim in ("Time", "Label", "Num"):
        if f"k{dim}" in res.coords and dim in ds.coords:
            res = res.assign_coords({dim: ds[dim]})
    for dim in ("Time", "Label", "Num"):
        if f"k{dim}" in res.coords and f"n{dim}" in ds.coords and ds[f"n{dim}"].dims[0] in res.dims:
            res = res.assign_coords({f"n{dim}": ds[f"n{dim}"].compute()})

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
            if lengths is None:
                index_array = index_frame[coord].to_dask_array(True)
                [lengths] = index_array.chunks
            else:
                index_array = index_frame[coord].to_dask_array(lengths)
        else:
            index_array = index_frame[coord].values
            if lengths is None:
                lengths = (index_array.size,)
        index_arrays += [index_array]
    return index_arrays, lengths


def groupby_apply(
    ds: xr.Dataset,
    by: Union[str, List[str]],
    fn: Callable,
    name_in: Union[str, List[str]],
    meta_out: Dict[str, Any] = None,
    sort: bool = True,
    index_arrays=None,
    length: int = None,
    **kwargs,
) -> Union[xr.DataArray, xr.Dataset]:
    """
    Mimic the groupby().apply mechanism of pandas using dask dataframe

    The result have *one* chunk per variable
    This will trigger some computation if `index_arrays` and `length` are not given


    Parameters
    ----------
    ds : dataset
    by : label or list of labels
        Used to determine the groups for the groupby.
    name_in : label or list of labels
        Name(s) of the source column(s) passed to fn.
    fn : callable
    meta_out : Dict[name, dtype]
        Description of the new columns added by fn (required for dask).
    sort: bool, default False
        Sort group keys. Get better performance by turning this off.
    length: int, optional
        size of the resulting arrays.
        If not given, will be computed.
    **kwargs
        Keyword arguments to pass to apply.

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

    if isinstance(name_in, str):
        name_in = [name_in]

    # if sort:
    #     ds = sortby(ds, by)

    if isinstance(ds, xr.DataArray):
        ds = ds.to_dataset(promote_attrs=True)

    relevent_dims = {str(dim) for var in ds.data_vars.values() for dim in var.dims}
    if len(relevent_dims) > 1:
        raise ValueError("You cannot groupby_apply with mixed shape Dataset")
    if not relevent_dims:
        raise ValueError("Nothing to group, your arrays are 0d")
    [k] = relevent_dims

    ds_only_k = ds.drop_dims([dim for dim in ds.dims if dim != k])
    ds_only_k = ds_only_k.rename(
        {coord: str(coord)[1:] for coord in ds_only_k.coords if str(coord).startswith("k")}
    )

    coords = [str(coord) for coord in ds_only_k.coords.keys()]
    content = [str(var) for var in ds_only_k.data_vars.keys()] + coords
    keep = list(set(name_in + by + coords))
    remove = [var for var in content if var not in keep]

    ds_only_k = ds_only_k.drop_vars(remove)
    dask = k in ds_only_k.chunks

    coords_meta = {coord: ds_only_k[coord].dtype for coord in coords}
    meta_by = {b: ds_only_k[b].dtype for b in by if b not in coords}
    if meta_out is None:
        meta_out = {}
    meta_out = {**meta_by, **meta_out}
    meta = {**coords_meta, **meta_out}

    @wraps(fn)
    def fn_(df, *args, **kwargs):
        if k != "k":
            df = df.set_index(k)
        res = fn(df, *args, **kwargs)
        if isinstance(res, pd.Series):
            res = res.to_frame()
        for col in meta_out.keys():
            if col in res.index:
                res = res.T
        cols = list(chain(res.index.names, res.columns))
        for coord in coords:
            if coord not in cols:
                res[coord] = df[coord]
        res = res.reset_index()
        res = res[[key for key in meta.keys() if key in res]]

        return res

    res: Union[xr.DataArray, xr.Dataset] = xr.Dataset()
    if dask:
        df = xarray_to_ddframe(ds_only_k).repartition(npartitions=1)
        # if k != "k":
        #     df = df.reset_index()
        df_gb = df.groupby(by=by, sort=sort)
        res_df = df_gb.apply(fn_, meta=meta, **kwargs)

        lengths = (ds[name_in[0]].size,) if length is None else length

        for name_out in meta_out:
            res[name_out] = (k,), res_df[name_out].to_dask_array(lengths)
        for coord in coords:
            res.coords[coord] = (k,), res_df[coord].to_dask_array(lengths)
    else:
        df = xarray_to_frame(ds_only_k, multi=False)
        if k != "k":
            df = df.reset_index()
        df_gb = df.groupby(by=by, sort=sort, group_keys=False)
        res_df = df_gb.apply(fn_, **kwargs)

        for name_out in meta_out:
            res[name_out] = (k,), res_df[name_out].values
        for coord in coords:
            res.coords[coord] = (k,), res_df[coord].values

    if "k" in res.dims:
        for coord_ in res.coords:
            coord = str(coord_)
            if f"k{coord}" in ds.coords:
                res = res.rename(cast(Mapping[Hashable, Hashable], {coord: f"k{coord}"}))

    for dim in ("Time", "Label", "Num"):
        if f"k{dim}" in res.coords and dim in ds.coords:
            res = res.assign_coords({dim: ds[dim]})
    for dim in ("Time", "Label", "Num"):
        if f"k{dim}" in res.coords and f"n{dim}" in ds.coords and ds[f"n{dim}"].dims[0] in res.dims:
            res = res.assign_coords({f"n{dim}": ds[f"n{dim}"].compute()})

    datavars = list(res.data_vars.keys())
    if len(datavars) == 1:
        [data_var] = datavars
        res = res[data_var]

    return res
