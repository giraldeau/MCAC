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

"""Compute the normalized 3d autocorrelation graph of an aggregate."""

import ctypes
import multiprocessing as mp
import os
from queue import Empty
from typing import Any, Callable, Optional, Tuple

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from pymcac.tools.volume_surface import volume_surface


def autocorrelation3d(
    spheres: pd.DataFrame,
    aggregate: Optional[pd.Series] = None,
    nradius: int = 32,
    nsamples: int = 16,
    repart: Callable = np.geomspace,
    start: Optional[float] = None,
    end: Optional[float] = None,
    nprocs: Optional[int] = None,
    *args,
    **kwargs,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the normalized 3d autocorrelation graph of an aggregate."""
    cols = ["Posx", "Posy", "Posz", "Radius"]
    if aggregate is not None:
        # Get the spheres of the given aggregate
        t, label = aggregate.name

        # get the spheres of the given aggregate
        spheres = spheres[spheres["Label"] == label].loc[t]

    spheres = spheres[cols]

    # Define range of radius
    if start is None:
        start = spheres.Radius.min() / 2
    if end is None:
        if aggregate is not None:
            end = aggregate.Rmax * 2
        else:
            end = spheres.Radius.max() * 2

    radius = np.array([0] + list(repart(start, end, nradius)))

    # convert pandas to numpy
    spheres = spheres.values

    # compute volumes
    if nprocs is None or nprocs > 0:
        volumes = parallel_volumes_autoco(spheres, radius, nsamples, nprocs, *args, **kwargs)
    else:
        volumes = seq_volumes_autoco(spheres, radius, nsamples, *args, **kwargs)

    # intersection volume
    # using A∩B = A + B - A∪b
    # and normalization
    # noinspection PyTypeChecker
    volumes = 2 - volumes / (nsamples * volumes[0])

    # no more contact
    volumes[volumes < 0] = 0

    return radius[1:], volumes[1:]


def parallel_volumes_autoco(
    spheres: np.ndarray,
    lradius: np.ndarray,
    nsamples: int,
    nprocs: Optional[int] = None,
    *args,
    **kwargs,
) -> np.ndarray:
    """Compute all the volumes needed for the autocorrelation.

    including the volume of the initial aggregate
    """
    # prepare output array
    volumes = np.zeros(len(lradius))

    # useless to make a contiguous copy, Process will do it
    with Parallel(worker_job, spheres, nprocs) as p:

        # Put work inside the queue
        njobs = 0
        for i, radius in enumerate(lradius):
            for _ in range(nsamples):
                p.qin.put((i, (radius, args, kwargs)))
                njobs += 1

                # Only compute the volume of the initial aggregate once
                if i == 0:
                    break

        # Receive result
        for _ in tqdm(range(njobs)):
            i, volume = p.qout.get()
            volumes[i] += volume

    return volumes


def worker_job(q: mp.Queue, r: mp.Queue, spheres: np.ndarray) -> None:
    """Pick one job and do it."""
    first = True

    # Each process must have a different random seed
    np.random.seed(None)

    while True:
        try:
            # get work
            i, (radius, args, kwargs) = q.get(block=first)

            # work
            volume = single_volume_autoco(radius, spheres, *args, **kwargs)

            # put result
            r.put((i, volume))

            first = False
        except Empty:
            if q.empty() and q.qsize() == 0:
                # nothing else to do, shutting
                break


def seq_volumes_autoco(
    spheres: np.ndarray, lradius: np.ndarray, nsamples: int, *args, **kwargs
) -> np.ndarray:
    """Compute all the volumes needed for the autocorrelation.

    including the volume of the initial aggregate
    """
    # prepare output array
    volumes = np.zeros(len(lradius))

    # make a contiguous copy
    spheres = spheres.copy()

    # Put work inside the queue
    njobs = (len(lradius) - 1) * nsamples + 1

    it = 0
    with tqdm(total=njobs) as progress_bar:
        for i, radius in enumerate(lradius):
            for _ in range(nsamples):

                volumes[i] += single_volume_autoco(radius, spheres, *args, **kwargs)
                it += 1
                progress_bar.update(it)

                # Only compute the volume of the initial aggregate once
                if i == 0:
                    break

    return volumes


def random_vector() -> np.ndarray:
    """Return a random vector on the unit sphere."""
    theta = np.random.uniform(0, 2 * np.pi)
    phi = np.arccos(1 - 2 * np.random.uniform(0, 1))

    return np.array([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)])


def translated_union(vec: np.ndarray, spheres: np.ndarray) -> np.ndarray:
    """Compute the union of an aggregate and a translated copy of itself."""
    nspheres = spheres.shape[0]

    # duplication
    duplicate = np.vstack((spheres, spheres))

    # translation of one version
    duplicate[:nspheres, :3] += vec

    return duplicate


def single_volume_autoco(radius: float, spheres: np.ndarray, *args, **kwargs) -> float:
    """Compute the volume needed for the autocorrelation.

    Except for the initial volume (i == 0), compute the union of an
    aggregate and a translated copy of itself in a random direction
    """
    if abs(radius) == 0:
        aggregate = spheres
    else:
        # noinspection PyTypeChecker
        vec: np.ndarray = radius * random_vector()
        aggregate = translated_union(vec, spheres)

    # compute volume with volume_surface
    return volume_surface(aggregate, *args, **kwargs)[0]


class Parallel:
    """A context to do parallel computation."""

    __slots__ = ("target", "nprocs", "_shared_array", "shared_array", "qin", "qout", "processes")

    def __init__(self, target: Callable, array: np.ndarray, nprocs: Optional[int] = None) -> None:
        """Init."""
        self.target = target
        self.nprocs = nprocs

        if self.nprocs is None:
            self.nprocs = len(os.sched_getaffinity(0))

        # input data to shared numpy array
        self._shared_array = mp.Array(ctypes.c_double, int(array.size))
        self.shared_array = np.frombuffer(self._shared_array.get_obj())
        self.shared_array = self.shared_array.reshape(array.shape)
        self.shared_array[:] = array

        self.qin: mp.Queue = mp.Queue()
        self.qout: mp.Queue = mp.Queue()

        self.processes = [
            mp.Process(target=target, args=(self.qin, self.qout, self.shared_array))
            for _ in range(self.nprocs)
        ]

    def __enter__(self) -> "Parallel":
        """Start the processes."""
        for p in self.processes:
            p.start()
        return self

    def __exit__(self, exception_type: Any, exception_value: Any, traceback: Any) -> None:
        """End the processes."""
        # first kill them just in case
        for p in self.processes:
            p.terminate()

        for p in self.processes:
            p.join()


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    from pymcac import MCAC, validation_data_path

    simu = MCAC(validation_data_path / "pytest_data/")

    # Read all data
    Spheres, Aggregates = simu.spheres, simu.aggregates

    last_agg = Aggregates.iloc[-1]

    _radius, _volumes = autocorrelation3d(Spheres, last_agg)

    plt.figure()
    plt.loglog(_radius, _volumes, "-o")
    plt.show()
