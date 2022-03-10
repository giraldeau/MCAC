#!/usr/bin/env python3
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
from itertools import chain

import dask.array as da
import numpy as np
import pytest
from dask import delayed

from pymcac.tools.core.dask_tools import aligned_rechunk, not_aligned_rechunk
from pymcac.tools.core.sorting import sortby

from .generator import generate_dummy_aggregates_data
from .test_data import check_dask_consistency, check_data


def identical(ds1, ds2):
    if not ds1.identical(ds2):
        return False
    if not ds1.chunks == ds2.chunks:
        return False
    return True


def check_chunks(chunks, target):
    assert max(chunks) == target
    assert min(chunks) > 0
    assert target in chunks


@pytest.mark.parametrize("dask", [0, 1, 2])
def test_not_aligned_rechunk_no_change(dask):
    data = generate_dummy_aggregates_data(dask=dask)

    no_change = not_aligned_rechunk(data, chunks=data.chunks)
    assert identical(data, no_change)


@pytest.mark.parametrize("dask", [0, 1, 2])
@pytest.mark.parametrize("dim", ["k", "Time", "Label"])
def test_not_aligned_rechunk_one_chunk(dask, dim):
    data = generate_dummy_aggregates_data(dask=dask)

    one_chunked = not_aligned_rechunk(data, chunks={dim: 3})

    assert identical(data.compute(), one_chunked.compute())
    assert identical(one_chunked, one_chunked.unify_chunks())
    assert set(one_chunked.chunks) == {dim} | set(data.chunks)
    check_chunks(one_chunked.chunks[dim], 3)
    check_data(one_chunked, aligned=False)


@pytest.mark.parametrize("dask", [0, 1, 2])
@pytest.mark.parametrize("dims", [("k", "Time"), ("Time", "Label")])
def test_not_aligned_rechunk_two_chunk(dask, dims):
    data = generate_dummy_aggregates_data(nt=5, dask=dask)

    chunks = {dim: n for dim, n in zip(dims, (3, 4))}
    two_chunked = not_aligned_rechunk(data, chunks=chunks)

    assert identical(data.compute(), two_chunked.compute())
    assert identical(two_chunked, two_chunked.unify_chunks())
    assert set(two_chunked.chunks) == set(dims) | set(data.chunks)
    for dim, n in chunks.items():
        check_chunks(two_chunked.chunks[dim], n)
    check_data(two_chunked, aligned=False)


@pytest.mark.parametrize("dask", [0, 1, 2])
@pytest.mark.parametrize("dims", [("k", "Time"), ("Time", "Label")])
def test_not_aligned_rechunk_two_chunk_seq(dask, dims):
    data = generate_dummy_aggregates_data(nt=5, dask=dask)

    chunks = {dim: n for dim, n in zip(dims, (3, 4))}
    two_chunked = not_aligned_rechunk(data, chunks=chunks)
    chunk0, chunk1 = ({dim: n} for dim, n in chunks.items())
    two_chunked_seq = not_aligned_rechunk(not_aligned_rechunk(data, chunks=chunk0), chunks=chunk1)
    assert identical(two_chunked, two_chunked_seq)
    two_chunked_seq = not_aligned_rechunk(not_aligned_rechunk(data, chunks=chunk1), chunks=chunk0)
    assert identical(two_chunked, two_chunked_seq)
    check_data(two_chunked_seq, aligned=False)


@pytest.mark.parametrize("dims", [("k", "Time"), ("Time", "Label")])
def test_not_aligned_rechunk_multi(dims):
    data = generate_dummy_aggregates_data(nt=5)

    multi_chunk = data.copy()
    for i, var in enumerate(chain(multi_chunk.values(), multi_chunk.coords.values())):
        if i == 0:
            continue
        if "k" in var.dims:
            var.data = da.from_array(var.data, chunks=i + 1)
    assert identical(data, multi_chunk.compute())

    rechunked = not_aligned_rechunk(multi_chunk, k=3)

    assert identical(data, rechunked.compute())
    assert set(rechunked.chunks) == {"k"}
    check_chunks(rechunked.chunks["k"], 3)
    assert identical(rechunked, rechunked.unify_chunks())
    check_data(rechunked, aligned=False)


def test_not_aligned_rechunk_no_compute():
    data = generate_dummy_aggregates_data()

    data["trigger"] = ("k",), da.from_delayed(
        raise_if_computed(), dtype=int, shape=(data.sizes["k"],)
    )
    not_aligned_rechunk(data, chunks={"k": 3})


@pytest.mark.parametrize("dask", [1, 2])
def test_aligned_rechunk_idempotent(dask):
    data = generate_dummy_aggregates_data(sort_info=True, dask=dask)
    rechunked = aligned_rechunk(data, Time=max(data.chunks["Time"]))
    assert identical(data, rechunked)
    check_data(rechunked)


@pytest.mark.parametrize("dask", [1, 2])
def test_aligned_rechunk_infer_on(dask):
    data = generate_dummy_aggregates_data(sort_info=True, dask=dask)
    rechunked = aligned_rechunk(data)
    assert identical(data, rechunked)
    check_data(rechunked)


@pytest.mark.parametrize("dask", [0, 5])
def test_aligned_rechunk_time(dask):
    data_unchunked = generate_dummy_aggregates_data(sort_info=True, dask=dask)
    data_chunked = aligned_rechunk(data_unchunked, Time=2)
    assert identical(data_unchunked.compute(), data_chunked.compute())
    assert set(data_chunked.chunks) == {"k", "Time"}
    check_chunks(data_chunked.chunks["Time"], 2)
    check_data(data_chunked)


@pytest.mark.parametrize("dask", [0, 5])
def test_aligned_rechunk_label(dask):
    data_unchunked = sortby(generate_dummy_aggregates_data(dask=dask), "Label")
    check_data(data_unchunked)
    data_chunked = aligned_rechunk(data_unchunked, Label=2)
    assert identical(data_unchunked.compute(), data_chunked.compute())
    assert set(data_chunked.chunks) == {"k", "Label"}
    check_chunks(data_chunked.chunks["Label"], 2)
    check_data(data_chunked)


@pytest.mark.parametrize("dask", [0, 5])
def test_aligned_rechunk_other(dask):
    aggregates = generate_dummy_aggregates_data(dask=dask)
    if dask:
        chunks = aggregates.chunks
    aggregates["Np"] = ("k",), np.random.randint(1, 10, aggregates.sizes["k"])
    if dask:
        aggregates = not_aligned_rechunk(aggregates, chunks=chunks)

    data_unchunked = sortby(aggregates, "Np")

    data_chunked = aligned_rechunk(data_unchunked, Np=2)

    assert identical(data_unchunked.compute(), data_chunked.compute())
    assert set(data_chunked.chunks) == {"k"}

    limits_inf = da.map_blocks(
        lambda Np: np.array([Np.min()]), data_chunked["Np"].data, dtype=int
    ).compute()
    limits_sup = da.map_blocks(
        lambda Np: np.array([Np.max()]), data_chunked["Np"].data, dtype=int
    ).compute()
    assert np.all(limits_inf[1:] > limits_inf[:-1])
    assert np.all(limits_sup[1:] > limits_sup[:-1])
    assert np.all(limits_sup >= limits_inf)

    check_data(data_chunked)


def test_aligned_rechunk_args():
    data_unchunked = generate_dummy_aggregates_data()
    data_chunked = aligned_rechunk(data_unchunked, Time=2)
    data_chunked2 = aligned_rechunk(data_unchunked, "Time", Time=2)
    assert identical(data_chunked, data_chunked2)

    data_chunked2 = aligned_rechunk(data_unchunked, Time=data_chunked.chunks["Time"])
    assert identical(data_chunked, data_chunked2)

    data_chunked2 = aligned_rechunk(not_aligned_rechunk(data_unchunked, Time=2), Time=None)
    assert identical(data_chunked, data_chunked2)


def test_aligned_rechunk_mixed():
    data_unchunked = generate_dummy_aggregates_data()

    k_chunked = aligned_rechunk(data_unchunked, "Time", k=3)
    check_dask_consistency(k_chunked)

    mixed = not_aligned_rechunk(data_unchunked, k=k_chunked.chunks["k"])
    k_chunked2 = aligned_rechunk(mixed, "Time", k=None)
    assert identical(k_chunked2, k_chunked)


def test_aligned_rechunk_no_compute():
    data = generate_dummy_aggregates_data()

    data["trigger"] = ("k",), da.from_delayed(
        raise_if_computed(), dtype=int, shape=(data.sizes["k"],)
    )
    aligned_rechunk(data, Time=2)


@delayed
def raise_if_computed():
    raise ValueError("I should never be computed")
