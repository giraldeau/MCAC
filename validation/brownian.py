#!/usr/bin/env python
# coding: utf-8
from pathlib import Path
from time import time

import matplotlib.pyplot as plt
import pandas as pd

from pymcac import MCAC
from pymcac import groupby_apply, groupby_agg
from pymcac import progress_compute, DaskDistribute


if __name__ == "__main__":

    # The folder with all .h5 and .xmf files
    result_path = Path("/Data/WORK/Projets/SRC/MCAC/validation/brownian_data")
    MCACSimulation = MCAC(result_path)

    start = time()
    # small data -> faster to avoid distribute
    # but it should work with this
    # with DaskDistribute():
    if True:
        # Read all data
        Aggregates = MCACSimulation.xaggregates
        Spheres = MCACSimulation.xspheres

        # small data -> faster to avoid dask
        # but it should work without this
        Aggregates, Spheres = progress_compute(Aggregates, Spheres)

        print(f"Reading time : {time() - start}")

        print(MCACSimulation.times.shape)

        k_B = 1.38066e-23  # Boltzmann constant in J/K
        T = 293.15  # Temperature

        D = k_B * T / float(Aggregates["f_agg"][0])
        BoxSize = float(Aggregates["BoxSize"][0])

        print(BoxSize, D)

        def distance(df: pd.DataFrame):
            Posx, Posy, Posz = df.iloc[0][["Posx", "Posy", "Posz"]]
            dx = abs(df["Posx"] - Posx) % BoxSize
            dy = abs(df["Posy"] - Posy) % BoxSize
            dz = abs(df["Posz"] - Posz) % BoxSize
            df["distance"] = dx ** 2 + dy ** 2 + dz ** 2
            return df

        tmp = groupby_apply(Spheres, by="Num",
                            meta_out={"distance": float}, fn=distance, name_in=["Posx", "Posy", "Posz"])
        distances = groupby_agg(tmp, by="Time", agg=[("distance", "mean", "distance")]).to_dataset()

        distances["theorical"] = 6 * D * distances.Time

        proper_time = groupby_agg(Aggregates, by="Time", agg=[
            ("min", "min", "proper_time"),
            ("mean", "mean", "proper_time"),
            ("max", "max", "proper_time"),
            ])

        # useless if done before
        # print("compute")
        # proper_time, distances = progress_compute(proper_time, distances)

    print(f"Total Compute time : {time() - start}")

    # because plotting a dataframe is easyer
    distances = distances.to_dataframe()
    proper_time = proper_time.to_dataframe()

    distances.plot(style=["-", "--"], figsize=(10, 6))
    plt.ylabel(r"$Distance^2\ (m^2)$", fontsize=20)
    plt.xlabel("Time (s)", fontsize=20)
    plt.legend(fontsize=16, loc=0)
    plt.show()

    proper_time.plot(style=["-", "--"], figsize=(10, 6))
    plt.ylabel("Proper time (s)", fontsize=20)
    plt.xlabel("Time (s)", fontsize=20)
    plt.legend(fontsize=16, loc=0)
    plt.show()
