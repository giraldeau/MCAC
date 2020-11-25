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
from typing import Hashable
from typing import List
from typing import Sequence
from typing import Tuple
from typing import Union

import dask.array as da
import numpy as np
import xarray as xr
from dask import delayed
from dask.delayed import Delayed
from numba import njit

from pymcac.tools.core.dask_tools import aligned_rechunk
from .sorting import sortby


def groupby_aggregate(Spheres: xr.Dataset,
                      Aggregates: xr.Dataset,
                      fn: Union[str, Callable],
                      sph_names_in: Union[str, List[str]],
                      agg_names_in: Union[str, List[str]] = None,
                      meta_out: Dict[str, Any] = None,
                      reduction=True,
                      dtype=None) -> xr.Dataset:
    """TODO use xarray map block and remove intermediate dtype"""
    if not reduction:
        raise NotImplementedError
    if isinstance(sph_names_in, str):
        sph_names_in = [sph_names_in]
    if agg_names_in is None:
        agg_names_in = []
    elif isinstance(agg_names_in, str):
        agg_names_in = [agg_names_in]

    for i, name_in in enumerate(sph_names_in):
        if name_in in ("Time", "Num"):
            sph_names_in[i] = "k" + name_in
    for i, name_in in enumerate(agg_names_in):
        if name_in in ("Time", "Label"):
            agg_names_in[i] = "k" + name_in

    if set(agg_names_in) & set(agg_names_in):
        raise ValueError("You cannot have twice the same variable name")

    if meta_out is None:
        meta_out = {name_in: Spheres[name_in].dtype for name_in in sph_names_in}
    if dtype is None:
        dtypes = set(meta_out.values())
        if len(dtypes) == 1:
            dtype = dtypes.pop()
        else:
            dtype = np.float64

    Spheres = sortby(Spheres, ["Time", "Label"])
    Aggregates = sortby(Aggregates, ["Time", "Label"])

    if reduction:
        Aggregates = aligned_rechunk(Aggregates, "Time", k=None)
        Spheres = aligned_rechunk(Spheres, Time=Aggregates.chunks["Time"])
        chunks = Aggregates.chunks["k"]
    else:
        Spheres = aligned_rechunk(Spheres, "Time", k=None)
        Aggregates = aligned_rechunk(Aggregates, Time=Spheres.chunks["Time"])
        chunks = Spheres.chunks["k"]

    arrs_in = ([Spheres[name_in].data for name_in in sph_names_in] +
               [Aggregates[name_in].data for name_in in agg_names_in])

    argin_string = "x" + ", x".join(map(str, range(len(arrs_in) + 1))) + ", fn, nout, intermediate_dtype"
    argout_string = "x0, (x" + ", x".join(map(str, range(1, len(arrs_in) + 1))) + ",), fn, nout, intermediate_dtype"

    res = da.map_blocks(eval(f"lambda {argin_string}: groupby_aggregate_in_block({argout_string})"),
                        Aggregates.Np.data,
                        *arrs_in,
                        fn=fn, nout=len(meta_out), intermediate_dtype=dtype,
                        dtype=dtype, new_axis=[0], chunks=(len(meta_out), chunks))

    res_ds = xr.Dataset()
    for i, (col, dtype) in enumerate(meta_out.items()):
        res_ds[col] = ("k",), res[i, :].astype(dtype)

    if reduction:
        res_ds = res_ds.assign_coords(Aggregates.coords)
    else:
        res_ds = res_ds.assign_coords(Spheres.coords)

    datavars = list(res_ds.data_vars.keys())
    if len(datavars) == 1:
        (data_var,) = datavars
        res_ds = res_ds[data_var]

    return res_ds


def groupby_aggregate_in_block(npagg: np.ndarray,
                               arrs: Tuple[np.ndarray, ...],
                               fn, nout, intermediate_dtype) -> np.ndarray:
    res = np.empty((nout, npagg.size), dtype=intermediate_dtype)
    start = 0
    for iagg, nsph in enumerate(npagg):

        end = start + nsph

        subsets = (arr[start:end] for iarr, arr in enumerate(arrs))
        res[:, iagg] = fn(*subsets)

        start = end

    return res


def groupby2(ds: xr.Dataset,
             group: str,
             op: str,
             op_args=(),
             op_kwargs=None,
             squeeze=True,
             restore_coord_dims=None,
             template="reduction",
             keep_coord=None) -> xr.Dataset:
    """TODO avoid xarray groupby"""
    if op_kwargs is None:
        op_kwargs = {}

    kgroup = group if "k" in ds[group].dims else "k" + group

    if isinstance(ds, xr.DataArray):
        ds = ds.to_dataset()
    if isinstance(template, xr.DataArray):
        template = template.to_dataset()

    sorted_ds = sortby(ds, group)
    if "k" in sorted_ds.chunks:
        aligned_ds = aligned_rechunk(sorted_ds, group, k=None)
    else:
        aligned_ds = sorted_ds

    if isinstance(template, str):
        if template == "reduction":
            template = xr.Dataset()
            if group in ds.dims:
                shape = aligned_ds[group].shape
                chunks = aligned_ds.chunks.get(group, None)
            else:
                if aligned_ds[group].chunks is None:
                    chunks = None
                    shape = np.unique(aligned_ds[group].data).sum(),
                else:
                    chunks = tuple(aligned_ds[group].data.map_blocks(
                            lambda ds: np.array([len(np.unique(ds))])
                            ).compute())
                    shape = sum(chunks),

            if chunks is None:
                meta = np.empty(shape)
            else:
                meta = da.empty(shape, chunks=chunks)
            for v, var in aligned_ds.data_vars.items():
                if "k" in var.dims:
                    template[v] = (kgroup,), meta.astype(var.dtype)
            if kgroup in aligned_ds.coords:
                template[kgroup] = (kgroup,), aligned_ds[group].data
        elif template == "transform":
            template = aligned_ds

    ds_only_k = aligned_ds.drop_dims([dim for dim in aligned_ds.dims if dim != "k"])
    if kgroup in template.dims:
        template = template.rename_dims({kgroup: "k"})

    def _groupby_block(ds):
        if keep_coord is not None:
            ds = ds.reset_coords(keep_coord)
        gb = ds.groupby(kgroup, squeeze=squeeze, restore_coord_dims=restore_coord_dims)
        func = getattr(gb, op)
        res = func(*op_args, **op_kwargs)
        if isinstance(res, xr.DataArray):
            res = res.to_dataset()
        if kgroup in res.dims:
            res = res.rename_dims({kgroup: "k"})
        if keep_coord is not None:
            res = res.set_coords(keep_coord)
        return res

    res = ds_only_k.map_blocks(_groupby_block, template=template)

    if len(res.coords) == 1:
        for c, coord in res.coords.items():
            res = res.drop_vars(c)
            if c.startswith("k"):
                c = c[1:]
            res = res.rename_dims(k=c)
            res = res.assign_coords({c: coord.values})

    datavars = list(res.data_vars.keys())
    if len(datavars) == 1:
        [data_var] = datavars
        res = res[data_var]

    return res


def groupby(ds: xr.Dataset,
            group: str,
            func: Callable,
            meta: Union[xr.DataArray, xr.Dataset, Tuple],
            *args,
            nchunk: int = None,
            **kwargs) -> xr.Dataset:
    """WIP"""
    sorted_ds = sortby(ds, group, nchunk=nchunk)
    chunks = None
    if sorted_ds.chunks is not None and "k" in sorted_ds.chunks:
        chunks = sorted_ds.chunks["k"]

    ds_only_k = sorted_ds.drop_dims([dim for dim in sorted_ds.dims if dim != "k"])
    ds_only_k = ds_only_k.rename({coord: coord[1:] for coord in ds_only_k.coords if coord.startswith("k")})
    del ds_only_k.attrs["sort"]

    ds_other = sorted_ds.drop_dims("k")

    if isinstance(meta, xr.DataArray):
        meta = meta.to_dataset()
    if isinstance(meta, xr.Dataset):
        template = meta
    else:
        template = xr.Dataset()
        for name, dim, dtype in meta:
            if dim != group:
                template.coords[dim] = ("k",), ds_only_k.coords[dim].data
                dim = "k"
                if chunks is None:
                    template[name] = (dim,), np.empty(sorted_ds.sizes[dim], dtype=dtype)
                else:
                    template[name] = (dim,), da.empty(sorted_ds.sizes[dim], dtype=dtype, chunks=sorted_ds.chunks[dim])
            else:
                template.coords[dim] = (dim,), ds_other.coords[dim].data
                template[name] = (dim,), da.empty(sorted_ds.sizes[dim], dtype=dtype, chunks=-1)

    kname = "k"
    if "k" in template.dims:
        for c, coord in template.coords.items():
            if "k" in coord.dims:
                kname = c
                break

    label_limits = None
    if group in ds_other.dims:
        if group == "Time":
            limits = np.insert(np.cumsum(ds_other.N.values), 0, 0)
        else:
            limits = np.insert(np.cumsum(ds_other.Nt.values), 0, 0)

        label_limits = ((np.asscalar(label.values), start, end)
                        for label, start, end in zip(ds_other.coords[group], limits[:-1], limits[1:]))

    if ds_only_k.chunks is not None and "k" in ds_only_k.chunks:

        res = ds_only_k.map_blocks(_groupby_block,
                                   args=args,
                                   kwargs={"group": group, "func": func, "kname": kname, "label_limits": label_limits,
                                           **kwargs},
                                   template=template)
    else:
        ds_only_k.coords["k"] = np.arange(ds_only_k.sizes["k"])
        res = _groupby_block(ds_only_k, *args,
                             group=group, func=func, kname=kname, label_limits=label_limits, **kwargs)

    if "k" in res.dims:
        res = res.assign_coords(sorted_ds.coords)
        for var in ["N", "Nt", "itmax"]:
            res[var] = sorted_ds[var]

    return res


def _groupby_block(ds_only_k: xr.Dataset,
                   group: str,
                   func: Callable,
                   kname: str,
                   label_limits,
                   *args,
                   squeeze=True,
                   restore_coord_dims=None,
                   **kwargs) -> xr.Dataset:
    if label_limits is not None:
        gb = ((float(t), ds_only_k.isel(k=slice(start, end)))
              for t, start, end in label_limits)
        gb = filter(lambda x: x[1].sizes["k"] > 0, gb)
    else:

        try:
            kgroup = group
            group_data = ds_only_k[group]
        except KeyError:
            kgroup = "k" + group
            group_data = ds_only_k[kgroup]

        gb = ds_only_k.groupby(group=kgroup, squeeze=squeeze, restore_coord_dims=restore_coord_dims)

    def wrap(ds, label, *args, **kwargs):

        group_ds = ds.swap_dims({"k": kname}).drop_vars([group])
        group_ds.coords[group] = [label]

        res = func(group_ds, *args, **kwargs)

        if isinstance(res, xr.DataArray):
            res = res.to_dataset()

        for var, data in res.data_vars.items():
            if not data.dims:
                res[var] = data.expand_dims({group: [label]})

        if kname in res.dims:
            res = res.rename_dims({kname: "k"})

        return res

    applied = (wrap(ds, label, *args, **kwargs) for label, ds in gb)
    res_ds = xr.concat(applied, dim="k", data_vars="minimal", coords="minimal", compat="no_conflicts")

    return res_ds


def groupby_index(ds: xr.Dataset,
                  indexname: str,
                  fn: Union[str, Callable],
                  meta: Union[xr.DataArray, xr.Dataset],
                  nchunk: int = None) -> xr.Dataset:
    if nchunk is None:
        nchunk = len(ds.chunks)

    if isinstance(meta, xr.DataArray):
        meta = meta.to_dataset()

    chunked_todo = np.array_split(ds.coords[indexname].data, nchunk)

    res = [chunked_select_and_apply(ds, chunk, fn, indexname, meta) for chunk in chunked_todo]
    res_dim = [xr.concat([data.drop_dims([d for d in data.dims if d != dim]) for data in res], dim=str(dim))
               for dim in res[0].dims]
    res = xr.merge(res_dim)
    res.attrs = ds.attrs
    res.attrs["sort"] = indexname
    return res


def chunked_select_and_apply(ds, chunk, fn, indexname,
                             meta: xr.Dataset) -> xr.Dataset:
    ds.coords["k"] = da.arange(len(ds.k))

    idx = np.insert(ds.N.cumsum().values[:-1], 0, 0)
    mask = chunk.min() <= ds.N
    idx_k = (chunk + idx[mask][:, np.newaxis]).flatten()

    ds_k = ds.sel({"k": idx_k, indexname: chunk, "Time": ds.Time[mask]})

    delayed_res: Delayed = delayed_chunked_select_and_apply(ds_k, chunk, fn, indexname)

    res_ds = meta.copy()
    len_k = 1
    if "Time" in meta.dims:
        len_k *= ds_k.Time.size
    if indexname in meta.dims:
        len_k *= ds_k[indexname].size

    for var in meta.data_vars:
        if indexname in meta[var].dims:
            res_ds[var] = (indexname,), da.from_delayed(get_from_delayed(delayed_res, var),
                                                        shape=(ds_k[indexname].size,),
                                                        dtype=meta[var].dtype)
            res_ds[indexname] = ds_k[indexname]
        if "Time" in meta[var].dims:
            res_ds[var] = ("k",), da.from_delayed(get_from_delayed(delayed_res, var),
                                                  shape=(len_k,),
                                                  dtype=meta[var].dtype)
            res_ds["Time"] = ds_k.Time

    if "Time" in res_ds.dims:
        res_ds["N"] = ds_k.N
        missing = [coord for coord in ds_k.coords if coord not in ("Time", indexname, "k")]
        for coord in missing:
            res_ds.coords[coord] = ("k",), da.from_delayed(get_from_delayed(delayed_res, coord),
                                                           shape=(len_k,),
                                                           dtype=ds_k.coords[coord].dtype)
    return res_ds


@delayed
def get_from_delayed(delayed: Delayed, name: Hashable) -> np.ndarray:
    return delayed[name].data


@delayed
def delayed_chunked_select_and_apply(ds: xr.Dataset, chunk: Sequence[int], fn: Callable, indexname: str) -> xr.Dataset:
    idx = np.insert(ds.N.cumsum().values[:-1], 0, 0)

    res = [select_and_apply(ds, k, idx, fn, indexname) for k in chunk]

    res_ds = xr.concat(res, dim="k", data_vars="minimal", coords="minimal", compat="no_conflicts")

    return res_ds


def select_and_apply(ds: xr.Dataset, k: int, idx: np.ndarray, fn: Callable, indexname: str) -> xr.Dataset:
    to_drop = [coord for coord in ds.coords if coord not in ("Time", indexname)]
    mask = k <= ds.N
    ds_k = ds.sel({"k": k + idx[mask], indexname: [k]})

    without_coords = ds_k.drop_vars(["N"] + to_drop).rename({"k": "Time"})

    res = fn(without_coords)

    if isinstance(res, xr.DataArray):
        res = res.to_dataset()
    for var, data in res.data_vars.items():
        if not data.dims:
            res[var] = data.expand_dims({indexname: [k]})

    if "Time" in res.dims:
        res = res.rename_dims({"Time": "k"}).assign_coords(Time=ds_k.Time)
        res["N"] = ds_k.N
        for removed in to_drop:
            res.coords[removed] = ds_k.coords[removed]

    return res


@njit(nogil=True, cache=True)
def groupby_index_in_block(npagg: np.ndarray,
                           time: np.ndarray,
                           naggs: np.ndarray,
                           nsph: np.ndarray,
                           sph_time: np.ndarray,
                           sph_label: np.ndarray,
                           arrs: Tuple[np.ndarray, ...],
                           fn) -> np.ndarray:
    res = np.empty((npagg.size, 6), dtype=np.float64)
    npcum = nsph.cumsum()
    kagg = 0
    for i, (t, nagg) in enumerate(zip(time, naggs)):

        start = npcum[i - 1] if i > 0 else 0
        end = npcum[i]

        l_sph_time = sph_time[start:end]
        l_sph_label = sph_label[start:end]
        l_arrs = [arr[start:end] for arr in arrs]

        for l in np.arange(nagg):
            subset = np.empty((npagg[kagg], len(arrs)), dtype=np.float64)
            ksph = 0

            for j, (sph_t, sph_l) in enumerate(zip(l_sph_time, l_sph_label)):
                if sph_t == t and sph_l == l:
                    subset[ksph] = [l_arr[j] for l_arr in l_arrs]
                    ksph += 1
            res[kagg] = fn(subset)
            kagg += 1
    return res
