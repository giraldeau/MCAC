#!/usr/bin/env python3
# coding=utf-8
"""
Plot a 2d projection of the aggregate
"""

from typing import Optional, Callable, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pymcac.tools.discretize import discretize


def view_proj_agg(spheres: pd.DataFrame,
                  aggregate: Optional[pd.Series] = None,
                  resolution: int = 128,
                  alphagangue: float = 0.4,
                  projdir: str = 'z',
                  reduce: Callable = np.sum) -> None:
    """
        Plot a 2d projection of the aggregate in the direction x, y or z

        if reduce is np.sum (or np.mean) the aggregate will be "transparent"
        if reduce is np.max  the aggregate will be fully opaque

        alphagangue is a parameter allowing some gangue around the aggregate
    """
    data, grid = discretize(spheres, aggregate, resolution, alphagangue)

    view_proj_disc(data, grid, projdir, reduce)


def view_proj_disc(data: np.ndarray,
                   grid: Tuple[np.ndarray, np.ndarray, np.ndarray],
                   projdir: str = 'z',
                   reduce: Callable = np.sum) -> None:
    """
        Plot a 2d projection of the discretized aggregate in the direction x, y or z

        if reduce is np.sum (or np.mean) the aggregate will be "transparent"
        if reduce is np.max  the aggregate will be fully opaque

        alphagangue is a parameter allowing some gangue around the aggregate
    """
    x, y, z = grid

    # binarization
    # noinspection PyTypeChecker,PyUnresolvedReferences
    data = (data > 0.).astype(float)

    # start the figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # projection
    if projdir == 'z':
        data = - reduce(data, axis=2)
        x = x[:, :, 0]
        y = y[:, :, 0]
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    if projdir == 'y':
        data = - reduce(data, axis=1)
        x = z[:, 0, :]
        y = x[:, 0, :]
        ax.set_xlabel('z')
        ax.set_ylabel('x')
    if projdir == 'x':
        data = - reduce(data, axis=0)
        x = y[0, :, :]
        y = z[0, :, :]
        ax.set_xlabel('y')
        ax.set_ylabel('z')

    ax.contourf(x, y, data, cmap='gray')
    plt.show()


if __name__ == "__main__":
    from pathlib import Path
    import numpy as np
    from pymcac import MCAC

    # The folder with all .h5 and .xmf files
    data_dir = Path("python-analysis/output_dir/")

    # Read all data
    Spheres, Aggregates = MCAC(data_dir).read()

    last_agg = Aggregates.iloc[-1]

    view_proj_agg(Spheres, last_agg)

