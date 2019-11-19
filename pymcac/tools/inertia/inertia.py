#!/usr/bin/env python3
# coding=utf-8
"""
Compute gyration radius using discretisation
"""
from typing import Optional, Tuple

import numpy as np
import pandas as pd

from pymcac.tools.discretize import discretize


def inertia(spheres: pd.DataFrame,
            aggregate: Optional[pd.Series] = None,
            resolution: int = 128,
            alphagangue: float = 0.4) -> np.ndarray:
    """
        Compute gyration radius using discretisation

        alphagangue is a parameter allowing some gangue around the aggregate
    """
    data, grid = discretize(spheres, aggregate, resolution, alphagangue)

    return inertia_disc(data, grid)


def inertia_disc(data: np.ndarray,
                 grid: Tuple[np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    """
        Compute gyration radius using discretisation
    """
    X, Y, Z = grid

    # noinspection PyTypeChecker
    coords = np.where(data > 0)

    # centering
    x = X[coords] - X[coords].mean()
    y = Y[coords] - Y[coords].mean()
    z = Z[coords] - Z[coords].mean()

    tensor = np.zeros((3, 3))
    tensor[0, 0] = (y ** 2 + z ** 2).mean()
    tensor[1, 1] = (x ** 2 + z ** 2).mean()
    tensor[2, 2] = (x ** 2 + y ** 2).mean()

    tensor[0, 1] = tensor[1, 0] = -(x * y).mean()
    tensor[0, 2] = tensor[2, 0] = -(x * z).mean()
    tensor[1, 2] = tensor[2, 1] = -(y * z).mean()

    eig = np.linalg.eigvals(tensor)

    return eig


if __name__ == "__main__":
    from pathlib import Path
    import numpy as np
    from pymcac import MCAC

    # The folder with all .h5 and .xmf files
    data_dir = Path("python-analysis/output_dir/")

    # Read all data
    Spheres, Aggregates = MCAC(data_dir).read()

    last_agg = Aggregates.iloc[-1]

    inertia_agg = inertia(Spheres, last_agg)
    print("Gyration radius (code): ", last_agg.Rg)
    print("Gyration radius:        ", np.sqrt(np.sum(inertia_agg)/2))
    print("Anisotropy:             ", inertia_agg.max(initial=-np.inf)/inertia_agg.min(initial=-np.inf))
