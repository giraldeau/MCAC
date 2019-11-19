# coding=utf-8
# cython: language_level=3
# cython: initializedcheck=False
# cython: binding=True
# cython: nonecheck=False
# cython: boundscheck=False
# cython: wraparound=False


cdef extern from "SBLVolumeSurface.hpp":
    void compute_volume_surface(const double spherelist[],
                                const int nspheres,
                                double& volume,
                                double& area) nogil

def compute_volume_surface_sbl(const double[:, ::1]& spherelist):
    """Compute surface and volume of an aggregate using SBL"""
    cdef int nspheres = spherelist.shape[0]
    cdef long dim = spherelist.shape[1]

    if dim != 4:
        raise ValueError("We need 4 values: x, y, z and r")

    cdef double volume = 0
    cdef double area = 0

    with nogil:
        compute_volume_surface(&spherelist[0, 0], nspheres, volume, area)

    return volume, area
