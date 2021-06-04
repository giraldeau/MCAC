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
#:
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from pymcac import validation_data_path
from pymcac.reader.mcac_reader import MCAC

from .test_data import check_consistency, check_data


def test_read_metadata():
    metadata = MCAC(validation_data_path / "pytest_data").metadata

    ref = {
        "flux_surfgrowth": 0.0001,
        "u_sg": 5.55556e-08,
        "dfe": 1.78,
        "kfe": 1.3,
        "lambda": 4.98113e-07,
        "rpeqmass": 5.25563e-09,
        "gamma_": 1.378,
        "P [Pa]": 101300.0,
        "T [K]": 1700.0,
        "Mu": 5.662e-05,
        "Rho [kg/m3]": 1800.0,
        "Dpm [nm]": 10.0,
        "sigmaDpm [nm]": 1.2,
        "FV [ppt]": 1e-05,
        "L": 1.06741e-06,
        "N []": 20.0,
    }
    print(ref)
    print(metadata)

    assert ref == metadata


def test_read_xaggregate():
    xaggregates = MCAC(validation_data_path / "pytest_data").xaggregates
    check_data(xaggregates)


def test_read_xspheres():
    xspheres = MCAC(validation_data_path / "pytest_data").xspheres
    check_data(xspheres)


def test_read_xdata():
    sim = MCAC(validation_data_path / "pytest_data")
    xaggregates = sim.xaggregates
    xspheres = sim.xspheres
    check_consistency(xspheres, xaggregates)


# def test_read_xaggregate(nt, nobj, data_type, dask, full):
#     if data_type == "aggregates":
#         data = generate_dummy_aggregates_data(
#             nt=nt, nagg=nobj, sort_info=True, dask=dask, full=full
#         )
#     elif data_type == "spheres":
#         data = generate_dummy_spheres_data(nt=nt, nsph=nobj, sort_info=True, dask=dask, full=full)
#     else:
#         raise ValueError(f"data_type {data_type} not understood")

#     check_data(data)


# @pytest.mark.parametrize("nt", [1, 3])
# @pytest.mark.parametrize("nagg", [1, 5])
# @pytest.mark.parametrize("dask", [0, 1, 2])
# @repeat(10)
# def test_generate_both(nt, nagg, dask):
#     aggregates, spheres = generate_dummy_data(nt=nt, nagg=nagg, sort_info=True, dask=dask)

#     check_data(aggregates)
#     check_data(spheres)
#     check_consistency(spheres, aggregates)


# # TODO hypothesis
