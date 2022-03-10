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

"""Tools related to dask."""
from contextlib import contextmanager, nullcontext

import dask.array as da
import numpy as np
from dask import compute, persist

# noinspection PyProtectedMember
from dask.base import get_scheduler
from dask.cache import Cache
from dask.diagnostics import ProgressBar
from dask.distributed import Client, performance_report, progress

from pymcac.tools.core.various import get_idx_name


def progress_compute(*dfs):
    """Show a progress bar while computing the given data.

    Data can be dataframe, array, dataarray or datasets.
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

    if len(dfs) == 1:
        return res[0]

    return res


@contextmanager
def dask_distribute(n_workers=None, threads_per_worker=1, report="", jupyter=False):
    """Context for dask distribute."""
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
    """dask_distribute that can be started/stopped."""

    def __init__(self, n_workers=None, threads_per_worker=1, report=""):
        """Init."""
        self.gen = dask_distribute(
            n_workers=n_workers, threads_per_worker=threads_per_worker, report=report
        ).gen

    def start(self):
        """Start the dask_distribute context."""
        return next(self.gen)

    def stop(self):
        """Stop the dask_distribute context."""
        try:
            next(self.gen)
        except StopIteration:
            pass


def not_aligned_rechunk(ds, chunks=None, **chunks_dict):
    """Rechunck datarray without taking care about aligning the chunks."""
    if chunks is None:
        chunks = {}
    chunks = {**chunks, **chunks_dict}
    ds = ds.copy()
    for var in ds.values():
        if var.dims:
            [dim] = var.dims
            if (var_chunks := chunks.get(dim, None)) is not None:
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
            if (coord_chunks := chunks.get(dim, None)) is not None:
                if coord.chunks is None:
                    coord.data = da.from_array(coord.data, chunks=coord_chunks)
                else:
                    if isinstance(coord_chunks, tuple) and isinstance(coord_chunks[0], int):
                        coord_chunks = (coord_chunks,)
                    coord.data = coord.data.rechunk(coord_chunks)
    return ds


def _infer_on(chunks, attrs):
    """Infer on what align the chunks."""
    if len(chunks) == 1:
        [on] = chunks
        if on != "k":
            return on

    if "sort" in attrs:
        return attrs["sort"][0]

    raise ValueError("You must give the dimension to align to")


def aligned_rechunk(ds, on=None, chunks=None, **chunks_dict):
    """Rechunck datarray taking care that chunks are aligned on an index."""
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
            chunks = {on: tuple(chunksize for chunksize in np.diff(time_chunks) if chunksize > 0)}
        else:
            raise ValueError("I don't understand")

    if "k" in chunks:
        chunks.pop("k")
    assert len(chunks) == 1

    if isinstance(chunks[on], int):
        time_chunks = np.arange(0, vals.size + chunks[on], chunks[on])
        time_chunks[-1] = vals.size
        chunks[on] = tuple(chunksize for chunksize in np.diff(time_chunks) if chunksize > 0)
    else:
        time_chunks = np.cumsum((0, *chunks[on]))

    chunksizes = [nums[start:end].sum() for start, end in zip(time_chunks[:-1], time_chunks[1:])]
    chunks["k"] = (tuple(chunksize for chunksize in chunksizes if chunksize > 0),)
    return not_aligned_rechunk(ds, **chunks)
