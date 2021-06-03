# cython: language_level=3
# cython: initializedcheck=False
# cython: binding=True
# cython: nonecheck=False
# cython: boundscheck=False
# cython: wraparound=False
# coding=utf-8
"""
Compute the overlapping coefficient COV
"""

import numpy as np

cimport numpy as np
from libc.math cimport sqrt

# noinspection PyUnresolvedReferences

from cython.parallel import prange


cdef double _coverage(const double* X,
                      const double* Y,
                      const double* Z,
                      const double* R,
                      const int points) nogil:
    """
    Compute the coverage coefficient of an aggregate
    C version (data must be contiguous)
    """

    cdef double res, d
    cdef long n, end

    cdef Py_ssize_t i, j

    n = 0
    res = 0.
    for i in range(points):
        for j in range(i + 1, points):
            d = 1. - sqrt((X[i] - X[j]) ** 2. +
                          (Y[i] - Y[j]) ** 2. +
                          (Z[i] - Z[j]) ** 2.) / (R[i] + R[j])
            if d > 0:
                res += d
                n += 1
    if n > 0:
        return res / n
    else:
        return 0


def coverages_cython(const long[::1]& Aggs,
                     const double[::1]& X,
                     const double[::1]& Y,
                     const double[::1]& Z,
                     const double[::1]& R,
                     const long openmp):
    """
    Compute the coverage coefficient of all the aggregates
    Data *MUST* be sorted by time and by label
    """

    cdef Py_ssize_t nsph = len(X)
    cdef Py_ssize_t nAgg = len(Aggs)

    cdef double[::1] res = np.empty(nAgg, dtype=np.double)
    cdef long[::1] isph = np.empty(nAgg, dtype=np.int)

    cdef Py_ssize_t iisph = 0
    cdef Py_ssize_t iagg

    for iagg in range(nAgg):
        isph[iagg] = iisph
        iisph += Aggs[iagg]

    cdef Py_ssize_t istart
    cdef Py_ssize_t iend

    if openmp > 1:
        for iagg in prange(nAgg - 1, -1, -1,
                           nogil=True, schedule="dynamic", chunksize=1, num_threads=openmp):
            istart = isph[iagg]
            res[iagg] = _coverage(&X[istart],
                               &Y[istart],
                               &Z[istart],
                               &R[istart],
                               Aggs[iagg])
    else:
        for iagg in range(nAgg):
            istart = isph[iagg]
            res[iagg] = _coverage(&X[istart],
                                  &Y[istart],
                                  &Z[istart],
                                  &R[istart],
                                  Aggs[iagg])

    return np.asarray(res)


def label_argsort(const long[::1]& nspheres,
                  const long[::1]& values,
                  const long openmp):
    """
    Compute the indexes that would sort values
    for each block of nspheres values (time step)
    """
    cdef Py_ssize_t ntime = nspheres.shape[0]
    cdef Py_ssize_t ndata = values.shape[0]
    cdef long[::1] idx = np.empty(ndata, dtype=np.int)
    cdef Py_ssize_t itime, i
    cdef long start, end


    if openmp > 1:
        for itime in prange(0, ntime, nogil=True, num_threads=openmp):
            start = 0
            if itime > 0:
                start = nspheres[itime-1]
            end = nspheres[itime]
            mergesort(&idx[0], &values[0], start, end - 1)
    else:
        for itime in range(0, ntime):
            start = 0
            if itime > 0:
                start = nspheres[itime-1]
            end = nspheres[itime]
            mergesort(&idx[0], &values[0], start, end - 1)
            start = end

    return np.asarray(idx)


cdef void mergesort(long v[], const long a[], const long low, const long high) nogil:
    """
    Merge sort algorithm
    """
    if low == high:
        v[low] = low
        return

    cdef long m
    if low < high:
        m = (high + low) // 2
        mergesort(v, a, low, m)
        mergesort(v, a, m + 1, high)

        merge(v, a, low, m, high)
        return


cdef void merge(long v[], const long a[], const long low, const long mid, const long high) nogil:
    """
    Merge part of the merge sort algorithm
    """
    cdef long b[high - low + 1]

    cdef long i = low
    cdef long j = mid + 1
    cdef long k = 0

    while i <= mid and j <= high:
        if a[v[i]] <= a[v[j]]:
            b[k] = v[i]
            i += 1
        else:
            b[k] = v[j]
            j += 1
        k += 1

    while i <= mid:
        b[k] = v[i]
        k += 1
        i += 1

    while j <= high:
        b[k] = v[j]
        k += 1
        j += 1

    k -= 1
    while k >= 0:
        v[low + k] = b[k]
        k -= 1

