#!/usr/bin/env python3

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

"""Read the MCAC output files."""

from pathlib import Path
from typing import Dict, Iterable, Optional, Sequence, Union, cast

import numpy as np
import pandas as pd
import xarray as xr
from dask import array as da
from dask import dataframe as dd
from dask import delayed
from numba import njit

from pymcac.tools.core.dask_tools import aligned_rechunk
from pymcac.tools.core.dataframe import groupby_agg, xarray_to_ddframe, xarray_to_frame

from .advancement_reader import AdvancementReader
from .h5_reader import H5Reader
from .xdmf_reader import XdmfReader  # type: ignore

CHUNKSIZE = 1  # in Mo


class MCAC:
    """This object read the simulation results from MCAC."""

    __slots__ = ("dir", "_metadata", "_advancement", "_times")

    def __init__(self, datadir: Union[str, Path]) -> None:
        self.dir = Path(datadir)
        self._metadata: Optional[Dict[str, Union[bool, float]]] = None
        self._advancement: pd.DataFrame = None
        self._times: Optional[np.ndarray] = None

    @property
    def metadata(self) -> Dict[str, Union[bool, float]]:
        """Read the metadata of the simulation from one of the files."""
        if self._metadata is None:
            # usually spheres are smaller file to read
            files = list(self.dir.glob("Spheres*.xmf"))

            # usually the last one is smaller
            filename = files[-1]

            self._metadata = XdmfReader(filename).extract_metadata()

        return self._metadata

    @property
    def times(self) -> np.ndarray:
        """Read the metadata of the simulation from one of the files."""
        if self._times is None:
            try:
                self._times = self.advancement.index.values
            except FileNotFoundError:
                files = sorted(list(self.dir.glob("Aggregats*.xmf")))
                self._times = np.unique(
                    np.fromiter(
                        (time for file in files for time in H5Reader(file).read_time()),
                        dtype=float,
                    )
                )
        self._times = cast(np.ndarray, self._times)
        return self._times

    @property
    def advancement(self) -> pd.DataFrame:
        """Read the advancement file of the simulation."""
        if self._advancement is None:
            self._advancement = AdvancementReader.read_advancement(self.dir)

        return self._advancement

    @property
    def xaggregates(self) -> xr.Dataset:
        """Read all the data from the aggregates files.

        The result is a large xarray+dask dataset
        """
        return self.get_xaggregates()

    def get_xaggregates(
        self,
        variables: Iterable[str] = None,
        tmax: float = None,
        nt: int = None,
        time_steps: np.ndarray = None,
    ) -> xr.Dataset:
        """Read the selected data from the aggregates files.

        The result is a large xarray+dask dataset
        """
        files = sorted(list(self.dir.glob("Aggregats*.xmf")))
        xaggregates = self.read_data(
            files, indexname="Label", chunksize=50, tmax=tmax, nt=nt, time_steps=time_steps
        )
        chunks = xaggregates.Rg.data.rechunk(block_size_limit=CHUNKSIZE * 1024 * 1024).chunks
        xaggregates = aligned_rechunk(xaggregates, Time=len(chunks[0]))

        if not variables:
            return xaggregates

        coords = xaggregates.coords
        xaggregates = xaggregates[list(variables)].assign_coords(coords)

        return xaggregates

    # noinspection PyUnusedFunction
    @property
    def ddaggregates(self) -> dd.DataFrame:
        """Read all the data from the aggregates files.

        The result is a large dask dataframe indexed with time
        """
        return self.get_ddaggregates()

    def get_ddaggregates(
        self,
        variables: Iterable[str] = None,
        tmax: float = None,
        nt: int = None,
        time_steps: np.ndarray = None,
    ) -> dd.DataFrame:
        """Read the selected data from the aggregates files.

        The result is a large dask dataframe indexed with time
        """

        ds = self.get_xaggregates(variables, tmax, nt, time_steps)
        df = xarray_to_ddframe(ds)

        chunk_limits = np.cumsum((0, *ds.chunks[0]), dtype=int)
        chunk_limits[-1] -= 1
        divisions = ds["Time"][chunk_limits].values.tolist()

        return df.set_index("Time", sorted=True, divisions=divisions)

    @property
    def aggregates(self) -> pd.DataFrame:
        """Read all the data from the aggregates files.

        The result is a large panda multiindex dataframe in the form
        data.loc[(time, label), attribute]
        """
        return self.get_aggregates()

    def get_aggregates(
        self,
        variables: Iterable[str] = None,
        tmax: float = None,
        nt: int = None,
        time_steps: np.ndarray = None,
    ) -> pd.DataFrame:
        """Read the selected data from the aggregates files.

        The result is a large panda multiindex dataframe in the form
        data.loc[(time, label), attribute]
        """
        return xarray_to_frame(self.get_xaggregates(variables, tmax, nt, time_steps))

    @property
    def xspheres(self) -> xr.Dataset:
        """Read all the data from the spheres files.

        The result is a large xarray+dask dataset
        """
        return self.get_xspheres()

    def get_xspheres(
        self,
        variables: Iterable[str] = None,
        tmax: float = None,
        nt: int = None,
        time_steps: np.ndarray = None,
    ) -> xr.Dataset:
        """Read the selected data from the spheres files.

        The result is a large xarray+dask dataset
        """

        files = sorted(list(self.dir.glob("Spheres*.xmf")))
        _xspheres = self.read_data(files, chunksize=50, tmax=tmax, nt=nt, time_steps=time_steps)
        chunks = _xspheres.Radius.data.rechunk(block_size_limit=CHUNKSIZE * 1024 * 1024).chunks
        xspheres = aligned_rechunk(_xspheres, Time=len(chunks[0]))

        if not variables:
            return xspheres

        coords = xspheres.coords
        xspheres = xspheres[list(variables)].assign_coords(coords)

        return xspheres

    @property
    def ddspheres(self) -> dd.DataFrame:
        """Read all the data from the spheres files.

        The result is a large dask dataframe indexed with time
        """
        return self.get_ddspheres()

    def get_ddspheres(
        self,
        variables: Iterable[str] = None,
        tmax: float = None,
        nt: int = None,
        time_steps: np.ndarray = None,
    ) -> dd.DataFrame:
        """Read the selected data from the spheres files.

        The result is a large dask dataframe indexed with time
        """

        ds = self.get_xspheres(variables, tmax, nt, time_steps)
        df = xarray_to_ddframe(ds)

        chunk_limits = np.cumsum((0, *ds.chunks[0]), dtype=int)
        chunk_limits[-1] -= 1
        divisions = ds["Time"][chunk_limits].values.tolist()

        return df.set_index("Time", sorted=True, divisions=divisions)

    @property
    def spheres(self) -> pd.DataFrame:
        """Read all the data from the spheres files.

        The result is a large panda multiindex dataframe in the form
        data.loc[(time, label), attribute]
        """
        return self.get_spheres()

    def get_spheres(
        self,
        variables: Iterable[str] = None,
        tmax: float = None,
        nt: int = None,
        time_steps: np.ndarray = None,
    ) -> pd.DataFrame:
        """Read the selected data from the spheres files.

        The result is a large panda multiindex dataframe in the form
        data.loc[(time, label), attribute]
        """
        return xarray_to_frame(self.get_xspheres(variables, tmax, nt, time_steps))

    def read_data(
        self,
        files: Sequence[Union[str, Path]],
        indexname: str = "Num",
        chunksize: int = None,
        tmax: float = None,
        nt: int = None,
        time_steps: np.ndarray = None,
    ) -> xr.Dataset:
        """Merge the data of all the files into one large dataset."""
        if time_steps is None:
            if nt is not None:
                if tmax is None:
                    tmax = self.times[-1]
                time_steps = np.linspace(self.times[0], tmax, nt)

        if time_steps is not None:
            idx = np.unique((self.times < time_steps[:, np.newaxis]).argmin(axis=1))
            time_steps = self.times[idx]

        if time_steps is None:
            if tmax is not None:
                nt = (self.times <= tmax).sum()
            else:
                nt = len(self.times)
            time_steps = self.times[:nt]

        data = {
            time: data
            for file in files
            for time, data in H5Reader.read_file(
                file, indexname=indexname, times=time_steps, chunksize=chunksize
            ).items()
        }

        first_data = next(iter(data.values()))
        columns = first_data.keys()

        reversed_data = {col: [data[time][col] for time in data.keys()] for col in columns}
        del data

        data_vars = {}
        coords = {}
        for col, data_t in reversed_data.items():
            if col == "nTime":
                continue
            full_data = da.concatenate(data_t)

            if full_data.size == len(time_steps):
                dims = ("Time",)
            else:
                dims = ("k",)

            if col in ("Time", "nTime", "kTime", indexname, "n" + indexname):
                coords[col] = xr.DataArray(full_data, dims=dims, name=col)
            else:
                data_vars[col] = xr.DataArray(full_data, dims=dims, name=col)

        coords["k" + indexname] = coords.pop(indexname)

        # to_persist = ("Time", "kTime", "k" + indexname, "n" + indexname)
        # for k, persisted in zip(to_persist,
        #                         dask.persist(*(coords[k] for k in to_persist))):
        #     coords[k] = persisted

        coords[indexname] = xr.DataArray(
            da.arange(coords["k" + indexname].max() + 1), dims=indexname, name=indexname
        )

        nmax = int(coords["n" + indexname].max())  # .compute())

        for nt_data in reversed_data["nTime"]:
            nt_data.resize((nmax,), refcheck=False)
        Nt = sum(reversed_data["nTime"])
        coords["nTime"] = xr.DataArray(Nt, dims=indexname, name="nTime")

        # to_persist = ("nTime", indexname)
        # for k, persisted in zip(to_persist,
        #                         dask.persist(*(coords[k] for k in to_persist))):
        #     coords[k] = persisted

        ds = xr.Dataset(data_vars, coords, attrs={"sort": ["Time", indexname]})  # type: ignore

        if "BoxSize" in ds.data_vars and "k" in ds.BoxSize.dims:
            BoxSize = ds.BoxSize
            ds = ds.drop_vars("BoxSize")
            BoxSize = groupby_agg(
                BoxSize,
                by="Time",
                agg=[("BoxSize", "first", "BoxSize")],
                sort=True,
                index_arrays=ds.Time,
            )
            ds["BoxSize"] = BoxSize.chunk({"Time": ds.chunks["Time"]})

        if "BoxSize" in ds.data_vars:
            BoxSize = ds.BoxSize
            ds = ds.drop_vars("BoxSize")
            ds["BoxVolume"] = BoxSize**3

        try:
            for col in self.advancement.columns:
                if col not in ds.data_vars:
                    ds[col] = ("Time",), self.advancement[col][time_steps]
        except FileNotFoundError:
            pass

        if "BoxVolume" not in ds.data_vars:
            print("Warning, the box volume might not be accurate")
            V0 = self.metadata["L"] ** 3
            N0 = self.metadata["N []"]
            if indexname == "Label":
                # aggregates
                N = groupby_agg(
                    ds, "Time", [("N", "sum", "Np")], index_arrays=ds.Time, length=len(time_steps)
                )
            else:
                # spheres
                N = ds["nNum"]
            ds["BoxVolume"] = V0 * N / N0

        return ds


@delayed
@njit(nogil=True, cache=True)
def compute_Nt(N):
    return np.bincount(N)[:0:-1].cumsum()[::-1]
