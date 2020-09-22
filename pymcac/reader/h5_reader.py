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
Read the hdf5 part of MCAC output files
"""

from pathlib import Path
from typing import Dict
from typing import Sequence
from typing import Tuple
from typing import Union

import numpy as np
# import cupy as cp
from dask import array as da
from dask import delayed
from h5py import File as h5File

from .xdmf_reader import XdmfReader


class H5Reader:
    """
    Object containing all functions necessary to read a h5 file
    """
    __slots__ = ("filename",)

    def __init__(self, filename: Union[str, Path]) -> None:
        self.filename = Path(filename).with_suffix(".h5")

    @classmethod
    def read_file(cls,
                  filename: Union[str, Path],
                  chunksize: int = None,
                  indexname: str = None,
                  sparse: bool = False) -> Dict[Tuple[float],
                                                Dict[str, Union[da.Array, np.ndarray]]]:
        """
        Read the h5 file into a mulitple dask arrays (lazy)

        We need info from the xmf file, so this file is read
        """
        return cls(filename).read(indexname=indexname, chunksize=chunksize, sparse=sparse)

    def read(self,
             indexname: str = "Num",
             chunksize: int = None,
             sparse: bool = False) -> Dict[Tuple[float],
                                           Dict[str, Union[da.Array, np.ndarray]]]:
        """
        Read the h5 file into a mulitple dask arrays (lazy)

        We need info from the xmf file, so this file is read
        """

        h5_groups = XdmfReader(self.filename).extract_h5_groups()
        sizes = XdmfReader(self.filename).extract_sizes()

        max_npart = max(sizes.values())

        times = np.sort(np.fromiter(h5_groups.keys(), dtype=float))

        if chunksize is None:
            chunksize = len(times)
        nchunk = min(max(len(times) // chunksize, 1), len(times))

        times_chunks = np.array_split(times, nchunk)

        data: Dict[Tuple[float], Dict[str, Union[da.Array, np.ndarray]]] = {}
        for times_chunk in times_chunks:
            times_chunk = tuple(times_chunk)

            chunk_size = sum(sizes[time] for time in times_chunk)
            chunk_group = {}
            for time in times_chunk:
                step = h5_groups[time]
                for attrib, group in step.items():
                    chunk_group[attrib] = *chunk_group.get(attrib, tuple()), group

            data[times_chunk] = {col: data
                                 for attrib, groups in chunk_group.items()
                                 for col, data in self.read_multiple_array(groups, attrib, chunk_size).items()}

            if "BoxSize" in h5_groups[times_chunk[0]]:
                with h5File(str(self.filename), 'r') as file_h5:
                    BoxSizes = [da.full(sizes[time], file_h5[h5_groups[time]["BoxSize"]][0]) for time in times_chunk]

                data[times_chunk]["BoxSize"] = da.concatenate(BoxSizes)
            elif "ll_box" in data[times_chunk]:
                data[times_chunk]["BoxSize"] = data[times_chunk]["ll_box"]
                del data[times_chunk]["ll_box"]

            if sparse:
                data[times_chunk]["Time"] = np.fromiter((time for time in times_chunk), dtype=float)
                data[times_chunk][indexname] = np.arange(max_npart)
                data[times_chunk]["size"] = np.fromiter((sizes[time] for time in times_chunk), dtype=int)
            else:
                data[times_chunk]["Time"] = da.concatenate([da.full(sizes[time], time) for time in times_chunk])
                data[times_chunk][indexname] = da.concatenate([da.arange(sizes[time]) for time in times_chunk])

        return data

    def read_multiple_array(self,
                            groups: Sequence[str],
                            attrib: str,
                            size: int) -> Dict[str, da.Array]:
        """
        Read multiple blocks of the h5 file into a dask array (lazy)
        """
        dtype = h5File(self.filename, 'r')[groups[0]].dtype
        shape = (size,)
        columns = (attrib,)

        if attrib == "Positions":
            columns = ("Posx", "Posy", "Posz")
            shape = 3, size

        data = da.from_delayed(read_h5_arrays(self.filename, groups, shape, attrib,
                                              dask_key_name=f"{groups} ({attrib})\n from {self.filename} (delayed)"),
                               dtype=dtype, shape=shape, name=f"{groups} ({attrib})\n from {self.filename} (data)",
                               meta=np.empty(shape, dtype))
        # meta=cp.empty(shape, dtype))

        if len(columns) == 1:
            return {columns[0]: data}
        return {col: data[i] for i, col in enumerate(columns)}


@delayed(pure=True, traverse=False)
def read_h5_arrays(filename: Path,
                   datasets: Sequence[str],
                   final_shape: Tuple[int, ...] = (-1,),
                   varname: str = "Unknown") -> np.ndarray:
    """Read multiple blocks of data"""

    # print(f"reading {varname} in {filename} ({datasets})")

    nvars = np.product(final_shape[:-1], dtype=int)
    part_shape = -1, *final_shape[:-1]

    with h5File(filename, 'r') as h5file:

        dtype = h5file[datasets[0]].dtype
        res = np.empty(final_shape, dtype=dtype)
        end = 0
        for dataset in datasets:
            data = h5file[dataset]
            start, end = end, end + data.size // nvars

            res[..., start:end] = np.reshape(data, part_shape).T

    return res
    # return cp.array(res)
