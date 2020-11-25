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
from typing import List
from typing import Union

import dask.array as da
import numpy as np
import xarray as xr
from dask import delayed

from pymcac.tools.core.various import get_idx_name


def sort_on_dask_array(ds, by):
    if ds.chunks is None:
        return False

    for b in by:
        try:
            source_data = ds[f"k{b}"]
        except KeyError:
            try:
                source_data = ds[f"{b}"]
            except KeyError:
                raise ValueError("Sorting index must be present in the dataset")
        if source_data.chunks is not None:
            return True
    return False


def sortby(ds: xr.Dataset, by: Union[str, List[str]], force=False, nchunk=None):
    if isinstance(by, str):
        by = [by]

    if not force:
        if "sort" in ds.attrs:
            uniq_concat = list(dict.fromkeys(by + list(ds.attrs["sort"])))
            k = len(uniq_concat)
            for i, sorti in enumerate(reversed(ds.attrs["sort"])):
                if uniq_concat[k - 1] == sorti:
                    k -= 1
                else:
                    k = len(uniq_concat)
            k = min(k, len(by))
            by = uniq_concat[:k]
        else:
            by = list(dict.fromkeys(by))

    if not by:
        return ds.copy(deep=False)

    if "k" in ds.dims or sort_on_dask_array(ds, by):
        res = custom_sortby(ds, by, nchunk=nchunk)
    else:
        res = ds.sortby(by)

    sortlist = ds.attrs["sort"] if "sort" in ds.attrs else []
    res.attrs["sort"] = list(dict.fromkeys(by + sortlist))

    return res


def custom_sortby(ds: xr.Dataset, by: Union[str, List[str]], nchunk=None):
    res = ds.copy(deep=False)
    was_a_dataarray = False
    if isinstance(res, xr.DataArray):
        was_a_dataarray = True
        res = res.to_dataset(promote_attrs=True)

    idx = from_argsort(res, by)

    if nchunk is None and res.chunks is not None and "k" in res.chunks:
        nchunk = max(1, len(res.chunks["k"]))

    if nchunk is None:
        for var in res.data_vars.values():
            if "k" in var.dims:
                var.data = var.data[idx]
        for coord in res.coords.values():
            if "k" in coord.dims:
                coord.data = coord.data[idx]
    else:
        chunked_slice, by_chunks = get_chunk_idx(res, idx, nchunk, by[0])

        for var in res.data_vars.values():
            if "k" in var.dims:
                var.data = chunk_slice(var.data, chunked_slice)
        for coord in res.coords.values():
            if "k" in coord.dims:
                coord.data = chunk_slice(coord.data, chunked_slice)

        # for dim in ds.dims:
        #     if dim == "k":
        #         continue
        #     chunks = by_chunks if dim == by[0] else -1
        #     for v, var in res.data_vars.items():
        #         if dim in var.dims:
        #             if var.chunks is None:
        #                 var.data = da.from_array(var.data, chunks=chunks)
        #             else:
        #                 var.data = var.data.rechunk(chunks)
        for c, coord in res.coords.items():
            if c in ("nTime", "nLabel", "nNum"):
                if by[0] in coord.dims:
                    if coord.chunks is None:
                        coord.data = da.from_array(coord.data, chunks=by_chunks)
                    else:
                        coord.data = coord.data.rechunk(by_chunks)
                else:
                    res[c] = coord.compute()
                continue

    if was_a_dataarray:
        [name] = res.data_vars
        res = res[name]
    return res


def get_chunk_idx(ds, idx, nchunk, by):
    idx_name = get_idx_name(ds)
    if idx_name is not None:
        if by == "Time" and f"n{idx_name}" in ds.coords:
            return chunk_idx(idx, nchunk, ds[f"n{idx_name}"])
        elif by == idx_name and "nTime" in ds.coords:
            return chunk_idx(idx, nchunk, ds.nTime)

    return chunk_idx(idx, nchunk)


def chunk_idx(idx, nchunk, align_by=None):
    chunks = (tuple([idx.size // nchunk + 1 for _ in range(idx.size % nchunk)])
              + tuple([idx.size // nchunk for _ in range(idx.size % nchunk, nchunk)]),)
    chunked_idx = idx.rechunk(chunks)

    if align_by is not None:
        align_by = align_by.values
        target = np.cumsum((0, *chunked_idx.chunks[0]))
        sizes = align_by.cumsum()
        by_chunks = np.abs(target.reshape(1, -1) - sizes.reshape(-1, 1)).argmin(axis=0) + 1
        by_chunks[0] = 0
        chunks = {0: tuple([align_by[start:end].sum() for start, end in zip(by_chunks[:-1], by_chunks[1:])])}
        chunked_idx = idx.rechunk(chunks)
        by_chunks = tuple(np.diff(by_chunks)),
    else:
        by_chunks = -1

    arg = chunked_idx.map_blocks(np.argsort)
    chunk_idx = arg.map_blocks(lambda x: np.arange(len(x))[x.argsort()])
    tmp = da.map_blocks(lambda x, y: x[y], chunked_idx, arg, dtype=chunked_idx.dtype)

    return (tmp, chunk_idx), by_chunks


def chunk_slice(arr, chunked_slice):
    chunks, chunks_idx = chunked_slice
    chunked_arr = arr.rechunk(-1)[chunks]
    return da.map_blocks(lambda x, y: x[y], chunked_arr, chunks_idx, dtype=arr.dtype)


def from_argsort(ds, by):
    dask = sort_on_dask_array(ds, by)

    idx = None
    for source in by[::-1]:

        try:
            source_data = ds[f"k{source}"]
        except KeyError:
            try:
                source_data = ds[f"{source}"]
            except KeyError:
                raise ValueError("Sorting index must be present in the dataset")

        if "k" in ds.dims and "k" not in source_data.dims:
            raise ValueError("Sorting index must be indexed by k")

        if idx is None or not dask:
            source_data = source_data.data
        else:
            if source_data.chunks is not None:
                source_data = source_data.data.rechunk(-1)
            else:
                source_data = da.from_array(source_data.data, chunks=-1)

        if idx is not None:
            source_data = source_data[idx]

        if dask:
            new_idx = da.from_delayed(argsort(source_data),
                                      dtype="i4", shape=source_data.shape,
                                      meta=np.empty(source_data.shape, "i4"))
        else:
            new_idx = source_data.argsort(kind="stable")

        if idx is None:
            idx = new_idx
        else:
            if dask:
                idx = idx.rechunk(-1)
            idx = idx[new_idx]

    if idx is None:
        if dask:
            idx = da.arange(ds.sizes["k"])
        else:
            idx = np.arange(ds.sizes["k"])
    return idx


@delayed
def argsort(data: np.ndarray) -> np.ndarray:
    return data.argsort(kind="stable")
