#!/usr/bin/env python3
# coding: utf-8
"""
Read the MCAC output files
"""

from pathlib import Path
from typing import Dict
from typing import Optional
from typing import Sequence
from typing import Union

import numpy as np
import pandas as pd
from dask import dataframe as dd

from .h5_reader import H5Reader
from .xdmf_reader import XdmfReader


class MCAC:
    """
    This object read and contains the simulation result from MCAC
    """
    __slots__ = ("dir",
                 "_metadata",
                 "_ddaggregates",
                 "_aggregates",
                 "_ddspheres",
                 "_spheres",
                 "times",)

    def __init__(self,
                 datadir: Union[str, Path]) -> None:
        self.dir = Path(datadir)
        self._metadata: Optional[Dict[str, Union[bool, float]]] = None
        self._ddaggregates: Optional[pd.DataFrame] = None
        self._aggregates: Optional[pd.DataFrame] = None
        self._ddspheres: Optional[pd.DataFrame] = None
        self._spheres: Optional[pd.DataFrame] = None
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
    def ddaggregates(self) -> dd.DataFrame:
        """
        Read all the data from the aggregates files

        The result is a large dask dataframe indexed with time
        """
        if self._ddaggregates is None:
            files = sorted(list(self.dir.glob("Aggregats*.xmf")))
            self._ddaggregates = self.read_data(files, index='Label')

        return self._ddaggregates

    @property
    def aggregates(self) -> pd.DataFrame:
        """
        Read all the data from the aggregates files

        The result is a large panda multiindex dataframe in the form data[attribut][time][label]
        """
        if self._aggregates is None:
            self._aggregates = self.ddaggregates.compute().reset_index().set_index(['Time', 'Label'])
        return self._aggregates

    @property
    def ddspheres(self) -> dd.DataFrame:
        """
        Read all the data from the spheres files

        The result is a large dask dataframe indexed with time
        """
        if self._ddspheres is None:
            files = sorted(list(self.dir.glob("Spheres*.xmf")))

            # Read data
            self._ddspheres = self.read_data(files)
        return self._ddspheres

    @property
    def spheres(self) -> pd.DataFrame:
        """
        Read all the data from the spheres files

        The result is a large panda multiindex dataframe in the form data[attribut][time][label]
        """
        if self._spheres is None:
            self._spheres = self.ddspheres.compute().reset_index().set_index(['Time', 'Num'])
        return self._spheres

    def read_data(self,
                  files: Sequence[Union[str, Path]],
                  index: str = "Num") -> dd.DataFrame:
        """
        Read all the data into a dask DataFrame
        """

        datas = [H5Reader.read_file(file, index) for file in files]

        self.times = np.concatenate([t for t, d in datas])
        return dd.concat([d for t, d in datas])
