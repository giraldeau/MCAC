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
Tools related to dask
"""
from contextlib import contextmanager, nullcontext

import dask.array as da
import numpy as np
import xarray as xr

# import cupy as cp
from dask import compute, persist

# noinspection PyProtectedMember
from dask.base import get_scheduler
from dask.cache import Cache
from dask.diagnostics import ProgressBar
from dask.distributed import Client, performance_report, progress

from pymcac.tools.core.various import get_idx_name


def progress_compute(*dfs):
    """
    Compute the given dataframe, array, dataarray or datasets showing a progress bar
    """
    scheduler = get_scheduler()
    # distribute = "Client" in type(scheduler).__name__ and hasattr(scheduler, "get")
    distribute = "Client" in str(scheduler)

    if distribute:
        futures = persist(*dfs)
        progress(*futures, notebook=False, multi=False)
        res = compute(*futures)
        print()
    else:
        with ProgressBar():
            res = compute(*dfs)

    # for df in res:
    #     for coordname, coord in df.coords.items():
    #         if isinstance(coord.data, cp.ndarray):
    #             coord.data = coord.data.get()
    #     for varname, var in df.data_vars.items():
    #         if isinstance(var.data, cp.ndarray):
    #             var.data = var.data.get()

    if len(dfs) == 1:
        return res[0]

    return res


@contextmanager
def dask_distribute(n_workers=None, threads_per_worker=1, report="", jupyter=False):
    """
    Context for dask distribute
    """
    if jupyter:
        print("Using a IPyParallel cluster seems to be a bad idea...")
        try:
            from ipyparallel import Client as ipyclient

            c = ipyclient()  # connect to IPyParallel cluster
            client = c.become_dask()  # start dask on top of IPyParallel
        except (OSError, ImportError):
            print("IPyParallel cluster not found, falling back to dask distribute")
            client = Client(n_workers=n_workers, threads_per_worker=threads_per_worker)
    else:
        client = Client(n_workers=n_workers, threads_per_worker=threads_per_worker)
    cache = Cache(2e9)
    if report:
        if isinstance(report, str):
            report_ctx = performance_report(filename=report)
        else:
            report_ctx = performance_report()
    else:
        report_ctx = nullcontext()
    with client as c, report_ctx:
        cache.register()
        try:
            yield c
        finally:
            cache.unregister()


class JupyterDaskDistribute:
    def __init__(self, n_workers=None, threads_per_worker=1, report=""):
        self.gen = dask_distribute(
            n_workers=n_workers, threads_per_worker=threads_per_worker, report=report
        ).gen

    def start(self):
        return next(self.gen)

    def stop(self):
        try:
            next(self.gen)
        except StopIteration:
            pass


# def _to_cpu(arr):
#     if isinstance(arr, np.ndarray):
#         return arr
#
#     assert isinstance(arr, da.Array)
#
#     if isinstance(arr._meta, np.ndarray):
#         return arr
#
#     if isinstance(arr._meta, cp.ndarray):
#         return arr.map_blocks(lambda x: x.get(), dtype=arr.dtype)
#     print(arr)
#     raise ValueError("What is this ?")
#
#
# def to_cpu(*dfs):
#     res = dfs[:]
#
#     for df in res:
#         for coordname, coord in df.coords.items():
#             coord.data = _to_cpu(coord.data)
#         for varname, var in df.data_vars.items():
#             var.data = _to_cpu(var.data)
#     return res


def not_aligned_rechunk(ds, chunks=None, **chunks_dict):
    # syntaxic sugar
    if chunks is None:
        chunks = {}
    chunks = {**chunks, **chunks_dict}
    ds = ds.copy()
    for var in ds.values():
        if var.dims:
            [dim] = var.dims
            var_chunks = chunks.get(dim, None)
            if var_chunks is not None:
                if var.chunks is None:
                    var.data = da.from_array(var.data, chunks=var_chunks)
                else:
                    if isinstance(var_chunks, tuple) and isinstance(var_chunks[0], int):
                        var_chunks = (var_chunks,)
                    var.data = var.data.rechunk(var_chunks)
    for c, coord in ds.coords.items():
        if coord.dims:
            [dim] = coord.dims
            if c == dim:
                continue
            coord_chunks = chunks.get(dim, None)
            if coord_chunks is not None:
                if coord.chunks is None:
                    coord.data = da.from_array(coord.data, chunks=coord_chunks)
                else:
                    if isinstance(coord_chunks, tuple) and isinstance(coord_chunks[0], int):
                        coord_chunks = (coord_chunks,)
                    coord.data = coord.data.rechunk(coord_chunks)
    return ds


def _infer_on(chunks, attrs):
    if len(chunks) == 1:
        [on] = chunks  # .keys()
        if on != "k":
            return on

    if "sort" in attrs:
        return attrs["sort"][0]

    raise ValueError("You must give the dimension to align to")


def aligned_rechunk(ds, on=None, chunks=None, **chunks_dict):
    # syntaxic sugar
    if chunks is None:
        chunks = {}
    chunks = {**chunks, **chunks_dict}

    idx_name = get_idx_name(ds)
    if idx_name is None or "Time" not in ds.dims:
        other = "Time" if idx_name is None else idx_name
        kchunks = chunks.get("k", None)
        otherchunks = chunks.get(other, kchunks)
        if kchunks is None:
            kchunks = otherchunks
        return not_aligned_rechunk(ds, {"k": kchunks, other: otherchunks})

    if on is None:
        on = _infer_on(chunks, ds.attrs)

    if "sort" in ds.attrs:
        assert (
            ds.attrs["sort"][0] == on
        ), f"The array should be sorted on on={on}, not on {ds.attrs['sort'][0]}"
    else:
        print("WARNING: you should sort your array before using aligned_rechunk")

    if on not in chunks:
        chunks[on] = None

    if on in ds.dims:
        vals = ds[on].values

        kon = "k" + on
        if on == "Time":
            nums = ds[f"n{idx_name}"].values
        elif on == idx_name:
            nums = ds.nTime.values
        else:
            vals, nums = np.unique(ds[kon].values, return_counts=True)
            raise ValueError("Not tested")
        if chunks[on] is None:
            chunks[on] = ds.chunks.get(on, None)
    else:
        vals, nums = np.unique(ds[on].values, return_counts=True)
        kon = on

    if chunks[on] is None and "k" in chunks:
        if chunks["k"] is None:
            # trusting the caller that k is already well aligned
            chunk = da.map_blocks(
                lambda coord_block: np.array([len(set(coord_block))]), ds[kon].data, dtype=int
            ).compute()
            chunks = {on: tuple(chunk)}
        elif isinstance(chunks["k"], int):
            # The given chunksize is a hint, no more
            nchunk = int(np.round(ds.k.size / chunks["k"]))

            time_chunks = np.searchsorted(
                nums.cumsum() / nums.sum(), np.linspace(0, 1, nchunk, endpoint=False)
            )
            time_chunks = np.append(time_chunks, [time_chunks[-1] + nums.size - time_chunks.max()])
            chunks = {on: tuple([chunksize for chunksize in np.diff(time_chunks) if chunksize > 0])}
        else:
            raise ValueError("I don't understand")

    if "k" in chunks:
        chunks.pop("k")
    assert len(chunks) == 1

    if isinstance(chunks[on], int):
        time_chunks = np.arange(0, vals.size + chunks[on], chunks[on])
        time_chunks[-1] = vals.size
        chunks[on] = tuple([chunksize for chunksize in np.diff(time_chunks) if chunksize > 0])
    else:
        time_chunks = np.cumsum((0, *chunks[on]))

    chunksizes = [nums[start:end].sum() for start, end in zip(time_chunks[:-1], time_chunks[1:])]
    chunks["k"] = (tuple([chunksize for chunksize in chunksizes if chunksize > 0]),)
    return not_aligned_rechunk(ds, **chunks)


# def broadcast_to(source, dest, nums=None, chunks=None):
#     """
#      * Time       -> Time-Num/Label (and transpose)
#      * Num/Label  -> Num/Label-Time  (and transpose)
#      * Label      -> Num-Label
#      * Time-Label -> Time-Num

#     #     T    -> T + N/L
#     #     N/L  -> T + N/L
#     #
#     # L    -> L+N
#     # L    -> T+L+N
#     # T+L  -> T+L+N
#     #
#     # 0 -> T
#     # 0 -> N/L
#     # 0 -> T + N/L
#     """
#     from pymcac.tools.core.sorting import sortby

#     if not source.dims:
#         if isinstance(dest, xr.DataArray):
#             relevent_dims = set(dest.dims)
#         else:
#             relevent_dims = {dim for var in dest.data_vars.values() for dim in var.dims}
#         if len(relevent_dims) > 1:
#             raise ValueError("You cannot broadcast to mixed shape Dataset")
#         if not relevent_dims:
#             return source
#         [k] = relevent_dims

#         res, _ = xr.broadcast(source, dest, exclude=set(dest.dims) - {k})
#         return res

#     missing_dims = set(dest.dims) - set(source.dims) - {"k"}
#     lost_dims = set(source.dims) - set(dest.dims)

#     if not missing_dims and not lost_dims:
#         return source

#     was_a_dataarray = False
#     if isinstance(source, xr.DataArray):
#         was_a_dataarray = True
#         source = source.to_dataset(promote_attrs=True)

#     if isinstance(dest, xr.DataArray):
#         dest = dest.to_dataset(promote_attrs=True)

#     common_dims = set(dest.dims) & set(source.dims) - {"k"}
#     orig_sort = None
#     if nums is None and missing_dims and common_dims and not lost_dims:
#         if len(missing_dims) > 1:
#             raise ValueError("adding more that one dim ?")
#         [missing_dim] = missing_dims

#         if f"n{missing_dim}" in dest.coords and set(dest[f"n{missing_dim}"].dims).issubset(
#             common_dims
#         ):
#             orig_sort = dest.attrs.get("sort", [])
#             common_dims = list(
#                 dict.fromkeys([dim for dim in orig_sort + list(common_dims) if dim in common_dims])
#             )

#             source = sortby(source, common_dims)
#             dest = sortby(dest, common_dims)

#             if any(dest.chunks):
#                 dest = aligned_rechunk(dest, common_dims[0], k=None)
#                 source = aligned_rechunk(
#                     source, chunks={common_dims[0]: dest.chunks[common_dims[0]]}
#                 )
#             elif any(source.chunks):
#                 source = aligned_rechunk(source, common_dims[0], k=None)

#                 dest = aligned_rechunk(dest, chunks={common_dims[0]: source.chunks[common_dims[0]]})

#             if nums is None:
#                 nums = dest[f"n{missing_dim}"]

#     if nums is None:
#         res = __broadcast_with_tag(source, dest)
#     else:
#         res = __broadcast_with_nums(source, dest, nums, chunks=chunks)

#     if orig_sort:
#         res = sortby(res, orig_sort)

#     if was_a_dataarray:
#         [name] = res.data_vars
#         res = res[name]

#     return res


# def __broadcast_with_nums(source, dest, nums, chunks=None):
#     """
#     * Time      -> Time-Num/Label (and transpose)
#     * Num/Label -> Num/Label-Time  (and transpose)
#     * Label     -> Num-Label
#     * Label     -> Num-Label
#     """

#     source_sort = source.attrs.get("sort", [None])[0]
#     dest_sort = dest.attrs.get("sort", [None])[0]
#     if (source_sort is None) or (dest_sort is None):
#         print("WARNING: You should sort you arrays before broadcasting")
#     else:
#         if source_sort != dest_sort:
#             raise ValueError("Your arrays must be aligned for broadcasting")
#         source_chunks = source.chunks.get(source_sort, None)
#         dest_chunks = dest.chunks.get(dest_sort, None)
#         if source_chunks != dest_chunks:
#             # if source_chunks is None:
#             #     source.chunk(dest_chunks)
#             # elif dest_chunks is None:
#             #     dest.chunk(source_chunks)
#             # else:
#             raise ValueError("Your arrays must be aligned for broadcasting")

#     res = dest.drop_vars([v for v in dest.data_vars])

#     if chunks is None:
#         chunks = dest.chunks.get("k", None)
#         chunks = (chunks,) if chunks is not None else None

#     def broadcast_var_block(arr, _nums):
#         return np.concatenate([np.repeat(val, n) for val, n in zip(arr, _nums)])

#     for v, var in source.data_vars.items():
#         if chunks is None:
#             res[v] = ("k",), broadcast_var_block(var.data, nums.data)
#         else:
#             res[v] = ("k",), da.map_blocks(
#                 broadcast_var_block,
#                 var.data,
#                 nums.data,
#                 meta=np.array((), dtype=var.dtype),
#                 chunks=chunks,
#             )

#     return res


# def __broadcast_with_tag(source, dest):
#     """
#     * Time      -> Time-Num/Label (and transpose)
#     * Num/Label -> Num/Label-Time  (and transpose)
#     * Label     -> Num-Label
#     * Label     -> Num-Label
#     """
#     from pymcac.tools.core.sorting import sortby

#     [tag] = set(source.dims) & set(dest.data_vars)

#     dest = dest.rename({tag: f"k{tag}"}).set_coords(f"k{tag}")
#     dest[tag] = source[tag]
#     dest = dest.rename(nTime="nTime_old")
#     dest.coords["nTime_old"] = dest.nTime_old.compute()
#     dest.coords["nTime"] = (tag,), np.unique(dest[f"k{tag}"].values, return_counts=True)[1]

#     orig_sort = dest.attrs.get("sort", None)
#     source = sortby(source, tag)
#     dest = sortby(dest, tag)

#     if "k" in dest.chunks:
#         dest = aligned_rechunk(dest, tag, k=None)
#         source = aligned_rechunk(source, chunks={tag: dest.chunks[tag]})
#     elif "k" in source.chunks:
#         source = aligned_rechunk(source, tag, k=None)
#         dest = aligned_rechunk(dest, chunks={tag: source.chunks[tag]})

#     res = dest.drop_vars([v for v in dest.data_vars])
#     idx = np.abs(source[tag].data.reshape(-1) - dest[f"k{tag}"].data.reshape(-1, 1)).argmin(axis=1)

#     for v, var in source.data_vars.items():
#         res[v] = ("k",), var.data[idx]

#     res = res.drop_dims(tag).reset_coords(f"k{tag}").rename({f"k{tag}": tag, "nTime_old": "nTime"})

#     if orig_sort is not None:
#         res = sortby(res, orig_sort)
#     return res
