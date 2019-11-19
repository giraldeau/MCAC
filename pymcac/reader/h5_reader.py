#!/usr/bin/env python3
# coding: utf-8
"""
Read the hdf5 part of MCAC output files

TODO use dask
"""

from pathlib import Path
from typing import Union, Dict, Optional

import numpy as np
import pandas as pd
from dask import dataframe as dd
from h5py import File as h5File

from .xdmf_reader import XdmfReader


class H5Reader:
    """
    Object containing all functions necessary to read an h5 file
    """
    __slots__ = ("filename", )

    def __init__(self, filename: Union[str, Path]) -> None:
        self.filename = Path(filename).with_suffix(".h5")

    @classmethod
    def read_file(cls,
                  filename: Union[str, Path],
                  index: Optional[str] = None) -> Dict[float, pd.DataFrame]:
        """
        Read the h5 file

        We need info from the xmf file, so this file is also read
        """
        return cls(filename).read(index=index)

    def read(self, index: Optional[str] = None) -> Dict[float, pd.DataFrame]:
        """
        Read the h5 file

        We need info from the xmf file, so this file is also read
        """

        h5_groups = XdmfReader(self.filename).extract_h5_groups()

        data: Dict[float, pd.DataFrame] = {}
        with h5File(str(self.filename), 'r') as file_h5:
            for time, step in h5_groups.items():
                datat = {}
                for attrib, group in step.items():
                    npdata = np.asarray(file_h5[group])
                    if attrib == "Positions":
                        pos = npdata.reshape((-1, 3))
                        datat['Posx'] = pos[:, 0]
                        datat['Posy'] = pos[:, 1]
                        datat['Posz'] = pos[:, 2]
                    else:
                        datat[attrib] = npdata
#                datat["Aggkey"] = list(map(lambda x: hash((time,x)),datat["Label"]))

                n = datat['Posx'].size
                datat["Time"] = np.repeat(datat["Time"], n)
                datat["BoxSize"] = np.repeat(datat["BoxSize"], n)

                del datat["Time"]

                data[time] = pd.DataFrame.from_dict(datat)

                if index is not None:
                    data[time].set_index(index, inplace=True)
        return data
