#!/usr/bin/env python3
# coding=utf-8

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
Ecriture du fichier shape.dat pour DDSCAT
"""

from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd

from pymcac.tools.coverage import coverages
from pymcac.tools.discretize import discretize
from pymcac.tools.inertia import inertia_disc
from pymcac.tools.volume_surface import volume_surface


def export_ddscat(
    spheres: pd.DataFrame,
    aggregate: pd.Series,
    filename: Union[str, Path] = "ddscat_shape.dat",
    resolution: int = 64,
    alphagangue: float = 0.4,
    aggregates: Optional[pd.DataFrame] = None,
    type_limits: Sequence[float] = (0.0,),
) -> None:
    """
    Ecriture du fichier shape.dat pour DDSCAT
    """
    if "cov" not in aggregate:
        if aggregates is None:
            raise ValueError(
                """In export_ddscat
I need the cov attribute
You can compute it before using the `coverages` function
Or give me the aggregates array so I can compute it for you"""
            )

        aggregates["cov"] = coverages(spheres, aggregates)
        aggregate = aggregates.loc[aggregate.name]

    data, grid = discretize(spheres, aggregate, resolution, alphagangue)
    x, y, z = grid
    nx, ny, nz = data.shape

    inertia_agg = inertia_disc(data, grid)
    anisotropy = inertia_agg.max(initial=-np.inf) / inertia_agg.min(initial=-np.inf)

    typed_data = np.zeros_like(data, dtype=int)
    for limit in type_limits:
        # noinspection PyTypeChecker
        typed_data += data > limit

    coords = np.where(typed_data > 0)
    nb_dipoles = len(coords[0])

    volume, _ = volume_surface(spheres, aggregate)

    aeff = (3 * volume / (4 * np.pi)) ** (1 / 3)

    with open(filename, "w") as f:
        f.write(f"Param alpha : {alphagangue:9.2e},")
        f.write(f" Cov : {aggregate['cov']:9.2e},")
        f.write(f" Rg (nm): {aggregate.Rg:9.2e},")
        f.write(f" Anisotropie : {anisotropy:9.2e},")
        f.write(f" aeff (nm) : {aeff:9.2e}\n")

        f.write(f"{nb_dipoles}\t= NAT\n")
        f.write("1.000000 0.000000 0.000000 = targer vector a1 (in TF)\n")
        f.write("0.000000 1.000000 0.000000 = targer vector a2 (in TF)\n")

        f.write(f"{1.}\t{1.}\t{1.}\t")
        f.write(" = lattice spacings (d_x,d_y,d_z)/d\n")

        f.write(f"{nx/2}\t{ny/2}\t{nz/2}\t")
        f.write(" = lattice offset x0(1-3)")
        f.write(" = (x_TF,y_TF,z_TF)/d for dipole 0 0 0\n")

        f.write("JA IX IY IZ ICOMP(x,y,z)\n")
        for i, coord in enumerate(zip(*coords)):
            f.write(f"{i}\t{x[coord]}\t{y[coord]}\t{z[coord]}")
            f.write(f"\t{typed_data[coord]}" * 3)
            f.write("\n")


if __name__ == "__main__":
    from pathlib import Path

    from pymcac import MCAC

    # The folder with all .h5 and .xmf files
    simu = MCAC(Path("python-analysis/output_dir/"))

    # Read all data
    Spheres, Aggregates = simu.spheres, simu.aggregates

    last_agg = Aggregates.iloc[-1]

    export_ddscat(Spheres, last_agg, aggregates=Aggregates, type_limits=[-0.2, 0.0])
