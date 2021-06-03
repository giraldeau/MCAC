#!/usr/bin/env python
# coding: utf-8
from pathlib import Path
from time import time

import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

from pymcac import MCAC, dask_distribute, progress_compute
from pymcac.tools.core.dataframe import groupby_agg, groupby_apply
from pymcac.tools.core.groupby import groupby2, groupby_aggregate
from pymcac.tools.physics.overlap import overlapping

if __name__ == "__main__":

    # The folder with all .h5 and .xmf files
    result_path = Path(
        # "/stockage/samba/Partages/public/MCAC_validation/03_VARYING_PP_polyd/03p1_SIGMAp_1p25/run1"
        "/Data/WORK/Projets/SRC/MCAC/validation/brownian_data"
    )

    MCACSimulation = MCAC(result_path)

    # small data -> faster to avoid distribute
    # but it should work with this
    with dask_distribute(report=True) as c:
        print(c)
        # import dask
        # with dask.config.set(scheduler='single-threaded'):
        # if True:
        start = time()

        # Read all data
        Aggregates = MCACSimulation.xaggregates
        # print(Aggregates)
        # print(Aggregates.compute())
        # print(Aggregates.kTime)
        # print(Aggregates.kTime.values)
        # ts = Aggregates.Time[:5]
        # mask_k = Aggregates.kTime.isin(ts)
        # print(mask_k)
        # print(mask_k.values)
        # exit()
        Spheres = MCACSimulation.xspheres
        # print(Spheres)
        # exit()
        # small data -> faster to avoid dask
        # but it should work without this
        # Aggregates, Spheres = progress_compute(Aggregates, Spheres)

        print(f"Reading time : {time() - start}")
        print(MCACSimulation.times.shape)

        k_B = 1.38066e-23  # Boltzmann constant in J/K
        T = 293.15  # Temperature

        D = k_B * T / float(Aggregates["f_agg"][0])
        BoxSize = float(Aggregates["BoxSize"][0])

        print(BoxSize, D)

        def fdistance(df: pd.DataFrame):
            Posx, Posy, Posz = df.iloc[0][["Posx", "Posy", "Posz"]]
            dx = abs(df["Posx"] - Posx) % BoxSize
            dy = abs(df["Posy"] - Posy) % BoxSize
            dz = abs(df["Posz"] - Posz) % BoxSize
            df["distance"] = dx ** 2 + dy ** 2 + dz ** 2
            return df

        # overlap = overlapping(Spheres, Aggregates)
        # print(overlap)
        # print(overlap.compute())

        # tmp = groupby_apply(
        #     Spheres,
        #     by="Num",
        #     meta_out={"distance": float},
        #     fn=fdistance,
        #     name_in=["Posx", "Posy", "Posz"],
        # )
        # print(tmp)
        # distances = groupby_agg(
        #     tmp, by="Time", agg=[("distance", "mean", "distance")], index_arrays=tmp.Time
        # ).to_dataset()

        BoxSize = 300.0

        def xdistance(ds: xr.DataArray):
            init = ds.isel(k=ds.kTime.argmin())
            dx = abs(ds.Posx - init.Posx) % BoxSize
            dy = abs(ds.Posy - init.Posy) % BoxSize
            dz = abs(ds.Posz - init.Posz) % BoxSize
            res = dx ** 2 + dy ** 2 + dz ** 2
            res = res.rename("distance")
            res["kTime"] = ds.kTime
            return res

        tmp = Spheres[["kTime", "kNum", "Posx", "Posy", "Posz", "Nt", "N"]]
        template = Spheres.Posx.rename("distance")
        distance = groupby2(
            tmp, "Num", op="map", op_args=(xdistance,), template=template
        ).to_dataset()
        distance["N"] = Spheres.N
        distance["Nt"] = Spheres.Nt
        distances = groupby2(distance, "Time", "mean")

        print(distances)
        print(distances.compute())

        distances["theorical"] = 6 * D * distances.Time

        proper_time = groupby_agg(
            Aggregates,
            by="Time",
            agg=[
                ("min", "min", "proper_time"),
                ("mean", "mean", "proper_time"),
                ("max", "max", "proper_time"),
            ],
        )

        # useless if done before
        # print("compute")
        # proper_time, distances = progress_compute(proper_time, distances)

        print(f"Total Compute time : {time() - start}")

    # because plotting a dataframe is easier
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
