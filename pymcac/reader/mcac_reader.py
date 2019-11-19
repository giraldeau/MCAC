#!/usr/bin/env python3
# coding: utf-8
"""
Read the MCAC output files

TODO use dask
"""

from pathlib import Path
from typing import Union, Dict, Optional, Sequence, Tuple
from threading import Thread, RLock
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
from dask import dataframe as dd

from .xdmf_reader import XdmfReader
from .h5_reader import H5Reader


class MCAC:
    """
    This object read and contains the simulation result from MCAC
    """
    __slots__ = ("dir",
                 "metadata",
                 "aggregates",
                 "spheres",
                 "times",
                 "seq",
                 "read_data",
                 "have_data")

    def __init__(self,
                 datadir: Union[str, Path],
                 seq: bool = False) -> None:
        self.dir = Path(datadir)
        self.metadata: Optional[Union[bool, float]] = None
        self.aggregates: Optional[pd.DataFrame] = None
        self.spheres: Optional[pd.DataFrame] = None
        self.times: Optional[np.ndarray] = None

        self.seq = seq
        if seq:
            self.read_data = self.read_data_seq
        else:
            self.read_data = self.read_data_par
        self.have_data = False

    def read_metadata(self) -> None:
        """
        Read the metadata of the simulation from one of the files
        """
        # usually spheres are smaller file to read
        files = list(self.dir.glob("Spheres*.xmf"))

        # usually the last one is smaller
        filename = files[-1]

        self.metadata = XdmfReader(filename).extract_metadata()

    def extract_time(self,
                     data: Dict[float, pd.DataFrame],
                     lock: Optional[RLock] = None) -> np.ndarray:
        """
        Extract time steps from data.
        Save them sorted and return them unsorted (for later reordering)
        """
        times = np.fromiter(data.keys(), dtype=float)

        # save sorted time steps
        if lock is not None:
            with lock:
                if self.times is None:
                    self.times = np.sort(times)
        else:
            if self.times is None:
                self.times = np.sort(times)
        return times

    # @profile
    @staticmethod
    def read_data_seq(files: Sequence[Union[str, Path]],
                      index: Optional[str] = None) -> Dict[float, pd.DataFrame]:
        """
        Read all the data into a dict(dict(ndarray))

        the result is in the form data[time][attribut][label]
        """
        args = [(file, index) for file in files]

        datas = [H5Reader.read_file(*arg) for arg in args]

        data = {k: v for d in datas for k, v in d.items()}

        return data

    @staticmethod
    def read_data_par(files: Sequence[Union[str, Path]],
                      index: Optional[str] = None) -> Dict[float, pd.DataFrame]:
        """
        Read all the data into a dict(dict(ndarray))

        the result is in the form data[time][attribut][label]
        """

        args = [(file, index) for file in files]

        with Pool(processes=cpu_count()) as pool:
            datas = pool.starmap(H5Reader.read_file, args)

            data = {k: v for d in datas for k, v in d.items()}

        return data

    # @profile
    def read_all_aggregates(self, lock: Optional[RLock] = None) -> None:
        """
        Read all the data from the aggregates files

        If not already stored, store the time steps

        The result is a large panda multiindex dataframe in the form data[attribut][time][label]
        """
        files = list(self.dir.glob("Aggregats*.xmf"))

        # Read data
        data = self.read_data(files, index='Label')
        self.extract_time(data, lock)

        # turn data into a large panda multiindex dataframe
        self.aggregates = pd.concat(data, names=['Time', 'Label']).sort_index(0)

    # @profile
    def read_all_spheres(self, lock: Optional[RLock] = None) -> None:
        """
        Read all the data from the spheres files

        If not already stored, store the time steps

        The result is a large panda multiindex dataframe in the form data[attribut][time][label]
        """
        files = list(self.dir.glob("Spheres*.xmf"))

        # Read data
        data = self.read_data(files)
        self.extract_time(data, lock)

        # turn data into a large panda multiindex dataframe
        self.spheres = pd.concat(data, names=['Time', 'Num']).sort_index(0)

    def read_seq(self) -> Tuple[Union[bool, float],
                                np.ndarray,
                                pd.DataFrame,
                                pd.DataFrame]:
        """
        Sequential version of the full reader
        """
        self.read_metadata()
        self.read_all_spheres()
        self.read_all_aggregates()

        return self.metadata, self.times, self.spheres, self.aggregates

    def read_par(self) -> Tuple[Union[bool, float],
                                np.ndarray,
                                pd.DataFrame,
                                pd.DataFrame]:
        """
        Parallel version of the full reader
        """
        self.read_metadata()

        lock = RLock()
        mcac = self

        class AggregatesJob(Thread):
            """Read all aggregates"""

            def run(self):
                """Read all aggregates"""
                mcac.read_all_aggregates(lock)

        aggregates_job = AggregatesJob()
        aggregates_job.start()

        class SphereJob(Thread):
            """Read all spheres"""

            def run(self):
                """Read all spheres"""
                mcac.read_all_spheres(lock)

        sphere_job = SphereJob()
        sphere_job.start()
        sphere_job.join()
        aggregates_job.join()

        return self.metadata, self.times, self.spheres, self.aggregates

    def read(self, dask: bool = False) -> Tuple[pd.DataFrame,
                                                pd.DataFrame]:
        """
        Read paraview output or preprocessed files if availables
        """
        data = self.dir / "data.h5"
        if not data.exists():
            if self.seq:
                self.read_seq()
            else:
                self.read_par()

            # self.Spheres.to_hdf(data, key="spheres", mode='w', format="table")
            # self.Aggregates.to_hdf(data, key="aggregats", mode='a', format="table")
            self.have_data = True

        if dask:
            # self.Spheres = pd.utils.from_pandas(self.Spheres, len(self.Spheres) / 2**16 + 1)
            # self.Aggregates = pd.utils.from_pandas(self.Aggregates, len(self.Aggregates) / 2**16 + 1)

            # self.Spheres = dd.from_pandas(self.Spheres.reset_index(0), len(self.Spheres) / 2**16 + 1)
            # self.Aggregates = dd.from_pandas(self.Aggregates.reset_index(0), len(self.Aggregates) / 2**16 + 1)

            self.spheres = dd.read_hdf(str(data), key="spheres")
            self.aggregates = dd.read_hdf(str(data), key="aggregats")
            self.have_data = True

        if not self.have_data:
            self.spheres = pd.read_hdf(data, key="spheres")
            self.aggregates = pd.read_hdf(data, key="aggregats")

        return self.spheres, self.aggregates
