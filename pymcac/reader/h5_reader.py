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
from typing import Optional
from typing import Sequence
from typing import Tuple
from typing import Union

import numpy as np
import pandas as pd
from dask import dataframe as dd
from dask import delayed
from h5py import File as h5File

from .xdmf_reader import XdmfReader


class H5Reader:
    """
    Object containing all functions necessary to read an h5 file
    """
    __slots__ = ("filename",)

    def __init__(self, filename: Union[str, Path]) -> None:
        self.filename = Path(filename).with_suffix(".h5")

    @classmethod
    def read_file(cls,
                  filename: Union[str, Path],
                  index: Optional[str] = None) -> Tuple[np.ndarray, dd.DataFrame]:
        """
        Read the h5 file into a dask dataframe (lazy)

        We need info from the xmf file, so this file is also read
        """
        return cls(filename).read(index=index)

    def read(self, index: str = "Num") -> Tuple[np.ndarray, dd.DataFrame]:
        """
        Read the h5 file into a dask dataframe (lazy)

        We need info from the xmf file, so this file is also read
        """

        h5_groups = XdmfReader(self.filename).extract_h5_groups()
        sizes = XdmfReader(self.filename).extract_sizes()

        data: Dict[float, pd.DataFrame] = {}
        for time, step in h5_groups.items():
            datas = [self.read_one(time, group, attrib, index, sizes)
                     for attrib, group in step.items()
                     if attrib not in ("BoxSize", "Time", index)]
            from functools import reduce
            data[time] = reduce(lambda x, y: x.join(y), datas)

            if "BoxSize" in step:
                with h5File(str(self.filename), 'r') as file_h5:
                    npdata = np.asarray(file_h5[step["BoxSize"]])[0]
                data[time]["BoxSize"] = npdata

        times = np.sort(np.fromiter(data.keys(), dtype=float))

        for time in times:
            data[time]["Time"] = time

        meta = {index: np.int64, **data[times[0]].dtypes.to_dict()}
        del meta["Time"]

        all_data = dd.from_delayed([set_index(data[time]) for time in times],
                                   meta=meta, divisions=(*times, times[-1]))

        return times, all_data

    def read_one(self,
                 time: float,
                 group: str,
                 attrib: str,
                 index: str,
                 sizes: Dict[float, int]) -> dd.DataFrame:
        """
        Read one block of the h5 file into a dask dataframe (lazy)
        """

        if attrib == "Positions":
            meta = {"Posx": np.float64,
                    "Posy": np.float64,
                    "Posz": np.float64}
            shape = (-1, 3)
        elif attrib in ("Label", "Np"):
            meta = {attrib: np.int64}
            shape = (-1,)
        else:
            meta = {attrib: np.float64}
            shape = (-1,)
        h5data = read_h5_frame(self.filename, group, list(meta.keys()), shape=shape, index_name=index)
        return dd.from_delayed([h5data], meta=meta, divisions=(0, sizes[time] - 1))


@delayed
def read_h5_frame(filename: Path,
                  dataset: str,
                  columns: Sequence[str],
                  shape: Tuple[int] = (-1,),
                  index_name: str = "Num") -> pd.DataFrame:
    """Read one block of data"""

    # print("reading", dataset, "in", filename)

    df = pd.DataFrame(np.asarray(h5File(filename, 'r')[dataset]).reshape(shape),
                      columns=columns)
    df.index.name = index_name
    return df


@delayed
def set_index(data: pd.DataFrame) -> pd.DataFrame:
    """set time as index"""

    res = data.reset_index().set_index("Time")
    return res
