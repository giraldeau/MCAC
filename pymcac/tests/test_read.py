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
import numpy as np
import pytest

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


@pytest.mark.parametrize("variables", [None, ["Rg", "Volume"]])
@pytest.mark.parametrize("tmax", [False, True])
@pytest.mark.parametrize("nt", [False, True])
@pytest.mark.parametrize("time_steps", [False, True])
def test_read_xaggregate(variables, tmax, nt, time_steps):
    sim = MCAC(validation_data_path / "pytest_data")

    _tmax = sim.times[len(sim.times) // 2]
    _nt = 3
    _time_steps = np.linspace(sim.times[0], _tmax, _nt)

    kwargs = {}
    if variables:
        kwargs["variables"] = variables
    if tmax:
        kwargs["tmax"] = _tmax
    if nt:
        kwargs["nt"] = _nt
    if time_steps:
        kwargs["time_steps"] = _time_steps
    print(kwargs)

    xaggregate = sim.get_xaggregates(**kwargs)
    check_data(xaggregate)

    if variables:
        assert list(xaggregate.data_vars.keys()) == variables
    if nt:
        assert len(xaggregate.Time) == _nt
    if tmax:
        assert xaggregate.Time.values[-1] == _tmax
    if time_steps:
        idx = np.unique((sim.times < _time_steps[:, np.newaxis]).argmin(axis=1))
        _time_steps = sim.times[idx]
        assert np.allclose(xaggregate.Time.values, _time_steps)


@pytest.mark.parametrize("variables", [None, ["Radius", "Label"]])
@pytest.mark.parametrize("tmax", [False, True])
@pytest.mark.parametrize("nt", [False, True])
@pytest.mark.parametrize("time_steps", [False, True])
def test_read_xspheres(variables, tmax, nt, time_steps):
    sim = MCAC(validation_data_path / "pytest_data")

    _tmax = sim.times[len(sim.times) // 2]
    _nt = 3
    _time_steps = np.linspace(sim.times[0], _tmax, _nt)

    kwargs = {}
    if variables:
        kwargs["variables"] = variables
    if tmax:
        kwargs["tmax"] = _tmax
    if nt:
        kwargs["nt"] = _nt
    if time_steps:
        kwargs["time_steps"] = _time_steps
    print(kwargs)

    xspheres = sim.get_xspheres(**kwargs)
    check_data(xspheres)

    if variables:
        assert list(xspheres.data_vars.keys()) == variables
    if nt:
        assert len(xspheres.Time) == _nt
    if tmax:
        assert xspheres.Time.values[-1] == _tmax
    if time_steps:
        idx = np.unique((sim.times < _time_steps[:, np.newaxis]).argmin(axis=1))
        _time_steps = sim.times[idx]
        assert np.allclose(xspheres.Time.values, _time_steps)


@pytest.mark.parametrize("variables_sph", [None, ["Radius", "Label"]])
@pytest.mark.parametrize("variables_agg", [None, ["Rg", "Np"]])
@pytest.mark.parametrize("tmax", [False, True])
@pytest.mark.parametrize("nt", [False, True])
@pytest.mark.parametrize("time_steps", [False, True])
def test_read_xdata(variables_sph, variables_agg, tmax, nt, time_steps):
    sim = MCAC(validation_data_path / "pytest_data")

    _tmax = sim.times[len(sim.times) // 2]
    _nt = 3
    _time_steps = np.linspace(sim.times[0], _tmax, _nt)

    kwargs = {}
    if tmax:
        kwargs["tmax"] = _tmax
    if nt:
        kwargs["nt"] = _nt
    if time_steps:
        kwargs["time_steps"] = _time_steps
    print(kwargs)

    xspheres = sim.get_xspheres(variables=variables_sph, **kwargs)
    xaggregates = sim.get_xaggregates(variables=variables_agg, **kwargs)
    check_data(xspheres)
    check_data(xaggregates)
    check_consistency(xspheres, xaggregates)

    if variables_sph:
        assert list(xspheres.data_vars.keys()) == variables_sph
    if variables_agg:
        assert list(xaggregates.data_vars.keys()) == variables_agg
    if nt:
        assert len(xspheres.Time) == _nt
        assert len(xaggregates.Time) == _nt
    if tmax:
        assert xspheres.Time.values[-1] == _tmax
        assert xaggregates.Time.values[-1] == _tmax
    if time_steps:
        idx = np.unique((sim.times < _time_steps[:, np.newaxis]).argmin(axis=1))
        _time_steps = sim.times[idx]
        assert np.allclose(xspheres.Time.values, _time_steps)
        assert np.allclose(xaggregates.Time.values, _time_steps)
