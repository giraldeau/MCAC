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
Read the MCAC output files
"""

from pathlib import Path
from typing import Dict
from typing import Optional
from typing import Sequence
from typing import Tuple
from typing import Union

import numpy as np
import pandas as pd
import xarray as xr
from dask import array as da
from dask import dataframe as dd
from dask import delayed

from pymcac.tools.core.dask_tools import aligned_rechunk
from pymcac.tools.core.dataframe import groupby_agg
from pymcac.tools.core.dataframe import xarray_to_ddframe
from pymcac.tools.core.dataframe import xarray_to_frame
from .h5_reader import H5Reader
from .xdmf_reader import XdmfReader
import dask
from numba import njit


CHUNKSIZE = 1  # in Mo


class MCAC:
    """
    This object read the simulation results from MCAC
    """
    __slots__ = ("dir",
                 "_metadata",
                 "times")

    def __init__(self,
                 datadir: Union[str, Path]) -> None:
        self.dir = Path(datadir)
        self._metadata: Optional[Dict[str, Union[bool, float]]] = None
        self.times: Optional[np.ndarray] = None

    @property
    def metadata(self) -> Dict[str, Union[bool, float]]:
        """
        Read the metadata of the simulation from one of the files
        """
        if self._metadata is None:
            # usually spheres are smaller file to read
            files = list(self.dir.glob("Spheres*.xmf"))

            # usually the last one is smaller
            filename = files[-1]

            self._metadata = XdmfReader(filename).extract_metadata()

        return self._metadata

    @property
    def xaggregates(self) -> xr.Dataset:
        """
        Read all the data from the aggregates files

        The result is a large xarray+dask dataset
        """
        files = sorted(list(self.dir.glob("Aggregats*.xmf")))

        _xaggregates = self.read_data(files, indexname='Label', chunksize=500)
        chunks = _xaggregates.Rg.data.rechunk(block_size_limit=CHUNKSIZE * 1024 * 1024).chunks
        return aligned_rechunk(_xaggregates, Time=len(chunks[0]))

    @property
    def sparse_xaggregates(self) -> xr.Dataset:
        """
        Read all the data from the aggregates files

        The result is a large xarray+dask+sparse dataset
        """
        files = sorted(list(self.dir.glob("Aggregats*.xmf")))

        _xaggregates = self.read_data_sparse(files, indexname='Label', chunksize=32)
        chunks = _xaggregates.Rg.data.rechunk(block_size_limit=CHUNKSIZE * 1024 * 1024).chunks
        for var in _xaggregates.values():
            var.data = var.data.rechunk(chunks)
        return _xaggregates

    def get_xaggregates(self,
                        variables: Tuple[str, ...] = tuple()) -> xr.Dataset:
        """
        Read the selected data from the aggregates files

        The result is a large xarray+dask dataset
        """
        if variables:
            return self.xaggregates[variables]
        return self.xaggregates

    # noinspection PyUnusedFunction
    @property
    def ddaggregates(self) -> dd.DataFrame:
        """
        Read all the data from the aggregates files

        The result is a large dask dataframe indexed with time
        """
        return self.get_ddaggregates()

    def get_ddaggregates(self,
                         variables: Tuple[str, ...] = tuple()) -> dd.DataFrame:
        """
        Read the selected data from the aggregates files

        The result is a large dask dataframe indexed with time
        """

        ds = self.get_xaggregates(variables)
        df = xarray_to_ddframe(ds)

        chunk_limits = np.cumsum((0, *ds.chunks[0]), dtype=int)
        chunk_limits[-1] -= 1
        divisions = ds["Time"][chunk_limits].values.tolist()

        return df.set_index("Time", sorted=True, divisions=divisions)

    @property
    def aggregates(self) -> pd.DataFrame:
        """
        Read all the data from the aggregates files

        The result is a large panda multiindex dataframe in the form data.loc[(time, label), attribute]
        """
        return self.get_aggregates()

    def get_aggregates(self,
                       variables: Tuple[str, ...] = tuple()) -> pd.DataFrame:
        """
        Read the selected data from the aggregates files

        The result is a large panda multiindex dataframe in the form data.loc[(time, label), attribute]
        """
        return xarray_to_frame(self.get_xaggregates(variables))

    @property
    def xspheres(self) -> xr.Dataset:
        """
        Read all the data from the spheres files

        The result is a large xarray+dask dataset
        """
        files = sorted(list(self.dir.glob("Spheres*.xmf")))
        _xspheres = self.read_data(files, chunksize=500)
        chunks = _xspheres.Radius.data.rechunk(block_size_limit=CHUNKSIZE * 1024 * 1024).chunks
        return aligned_rechunk(_xspheres, Time=len(chunks[0]))

    @property
    def sparse_xspheres(self) -> xr.Dataset:
        """
        Read all the data from the spheres files

        The result is a large xarray+dask+sparse dataset
        """
        files = sorted(list(self.dir.glob("Spheres*.xmf")))
        _xspheres = self.read_data_sparse(files, chunksize=500)
        chunks = _xspheres.Radius.data.rechunk(block_size_limit=CHUNKSIZE * 1024 * 1024).chunks
        for var in _xspheres.values():
            var.data = var.data.rechunk(chunks)
        return _xspheres

    def get_xspheres(self,
                     variables: Tuple[str, ...] = tuple()) -> xr.Dataset:
        """
        Read the selected data from the spheres files

        The result is a large xarray+dask dataset
        """
        if variables:
            return self.xspheres[variables]
        return self.xspheres

    @property
    def ddspheres(self) -> dd.DataFrame:
        """
        Read all the data from the spheres files

        The result is a large dask dataframe indexed with time
        """
        return self.get_ddspheres()

    def get_ddspheres(self,
                      variables: Tuple[str, ...] = tuple()) -> dd.DataFrame:
        """
        Read the selected data from the spheres files

        The result is a large dask dataframe indexed with time
        """

        ds = self.get_xspheres(variables)
        df = xarray_to_ddframe(ds)

        chunk_limits = np.cumsum((0, *ds.chunks[0]), dtype=int)
        chunk_limits[-1] -= 1
        divisions = ds["Time"][chunk_limits].values.tolist()

        return df.set_index("Time", sorted=True, divisions=divisions)

    @property
    def spheres(self) -> pd.DataFrame:
        """
        Read all the data from the spheres files

        The result is a large panda multiindex dataframe in the form data.loc[(time, label), attribute]
        """
        return self.get_spheres()

    def get_spheres(self,
                    variables: Tuple[str, ...] = tuple()) -> pd.DataFrame:
        """
        Read the selected data from the spheres files

        The result is a large panda multiindex dataframe in the form data.loc[(time, label), attribute]
        """
        return xarray_to_frame(self.get_xspheres(variables))

    def read_data(self,
                  files: Sequence[Union[str, Path]],
                  indexname: str = "Num",
                  chunksize: int = None) -> xr.Dataset:
        """
        Merge the data of all the files into one large dataset
        """

        # files = files[:1]

        datas = {time: data for file in files for time, data in
                 H5Reader.read_file(file, indexname=indexname, chunksize=chunksize).items()}
        self.times = np.fromiter((time for times in datas.keys() for time in times), dtype=float)

        first_datas = next(iter(datas.values()))
        columns = first_datas.keys()

        reversed_datas = {col: [datas[time][col] for time in datas.keys()] for col in columns}
        # i = 0
        # for ichunk, time in enumerate(datas.keys()):
        #     reversed_datas["iTime"][ichunk] = reversed_datas["iTime"][ichunk] + i
        #     reversed_datas["kiTime"][ichunk] = reversed_datas["kiTime"][ichunk] + i
        #     i += len(time)
        del datas

        data_vars = {}
        coords = {}
        for col, data_t in reversed_datas.items():
            if col == "Nt":
                continue
            data = da.concatenate(data_t)

            if data.size == len(self.times):
                dims = "Time",
            else:
                dims = "k",

            if col in ("Time", "kTime", indexname):  # "iTime", "kiTime"
                coords[col] = xr.DataArray(data, dims=dims, name=col)
            else:
                data_vars[col] = xr.DataArray(data, dims=dims, name=col)

        coords["k" + indexname] = coords.pop(indexname)
        coords[indexname] = xr.DataArray(da.arange(coords["k" + indexname].max() + 1), dims=indexname, name=indexname)

        # to_persist = ("Time", "kTime", indexname, "k" + indexname)
        # for k, persisted in zip(to_persist,
        #                         dask.persist(*(coords[k] for k in to_persist))):
        #     coords[k] = persisted
        # data_vars["N"] = data_vars["N"].persist()

        nmax = int(data_vars["N"].max())#.compute())
        # Nt = da.from_delayed(compute_Nt(data_vars["N"]),
        #                      shape=(nmax,), dtype=np.int64, meta=np.empty((nmax,), np.int64)).persist()
        # data_vars["Nt"] = xr.DataArray(Nt, dims=indexname, name="Nt")

        for nt_data in reversed_datas["Nt"]:
            nt_data.resize((nmax,), refcheck=False)
        Nt = sum(reversed_datas["Nt"])
        data_vars["Nt"] = xr.DataArray(Nt, dims=indexname, name="Nt")

        # tmax = da.from_array(coords["Time"][data_vars["Nt"] - 1], chunks=-1)
        # data_vars["tmax"] = xr.DataArray(tmax, dims=indexname, name="tmax")
        # itmax = da.from_array(coords["iTime"][data_vars["Nt"] - 1].values, chunks=-1)
        # data_vars["itmax"] = xr.DataArray(itmax, dims=indexname, name="itmax")

        ds = xr.Dataset(data_vars, coords, attrs={"sort": ["Time", indexname]})

        if "BoxSize" in ds.data_vars and "k" in ds.BoxSize.dims:
            BoxSize = ds.BoxSize
            ds = ds.drop_vars("BoxSize")
            BoxSize = groupby_agg(BoxSize, by="Time", agg=[("BoxSize", "first", "BoxSize")])
            ds["BoxSize"] = BoxSize.chunk({"Time": ds.chunks["Time"]})

        return ds

    def read_data_sparse(self,
                         files: Sequence[Union[str, Path]],
                         indexname: str = "Num",
                         chunksize: int = None) -> xr.Dataset:
        """
        Merge the data of all the files into one large sparse dataset
        """
        import sparse
        # files = files[:1]

        datas = {time: data for file in files for time, data in
                 H5Reader.read_file(file, indexname=indexname, chunksize=chunksize, sparse=True).items()}
        self.times = np.fromiter((time for times in datas.keys() for time in times), dtype=float)
        # test = self.times.size == np.unique(self.times).size

        max_npart = max(s for times_chunk, data in datas.items() for s in data["size"])

        for time, data in datas.items():
            shape = len(time), max_npart
            size = data["size"]
            del data["size"]
            if indexname in data:
                del data[indexname]
            del data["Time"]

            nnz = sum(size)
            it_coord = ((it, ip) for it, t in enumerate(time) for ip in range(size[it]))

            coords = da.from_delayed(compute_coord(it_coord,
                                                   dask_key_name=f"coords for {time} (delayed)"),
                                     dtype="i4", shape=(2, nnz), name=f"coords for {time} (data)",
                                     meta=np.empty((2, nnz), "i4"))

            for col in data:
                data[col] = da.from_delayed(to_sparse(coords, data[col], shape,
                                                      dask_key_name=f"{col} for {time} (delayed)"),
                                            dtype=data[col].dtype, shape=shape, name=f"{col} for {time} (data)",
                                            meta=sparse.COO.from_numpy(np.eye(1)))

        first_datas = next(iter(datas.values()))
        columns = first_datas.keys()

        reversed_datas = {col: [datas[time][col] for time in datas.keys()] for col in columns}
        concatenated_datas = {cols: da.concatenate(data_t) for cols, data_t in reversed_datas.items()}

        data_vars = {}
        coords = {"time": self.times, indexname: range(max_npart)}
        for col, data in concatenated_datas.items():
            if col in ("BoxSize",):
                data_vars[col] = ("time",), data[:, 0].map_blocks(lambda x: x.todense(), dtype=data.dtype)
            else:
                data_vars[col] = ("time", indexname), data
        return xr.Dataset(data_vars, coords)


@delayed
def compute_coord(it_coord) -> np.ndarray:
    """
    Compute coords for sparse dataset
    """
    return np.fromiter(it_coord, dtype="i4,i4").view(np.int32).reshape(-1, 2).T


@delayed
def to_sparse(coords: da.Array,
              data: da.Array,
              shape: Tuple[int, ...]):
    """
    Transform an 1d array into a sparse 2d array
    """
    import sparse
    fill_value = np.nan
    try:
        np.array([fill_value], dtype=data.dtype)
    except ValueError:
        fill_value = 0
    return sparse.COO(coords=coords, data=data, shape=shape, has_duplicates=False, sorted=True, fill_value=fill_value)


@delayed
@njit(nogil=True, cache=True)
def compute_Nt(N):
    return np.bincount(N)[:0:-1].cumsum()[::-1]
