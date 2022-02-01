#!/usr/bin/env python3
# coding=utf-8
"""Compile and install pymcac"""
from Cython.Build import cythonize

from numpy.distutils.core import Extension, setup

coverage = Extension(
    name="pymcac.tools.coverage.coverages_cython",
    sources=[
        "pymcac/tools/coverage/coverages_cython.pyx",
    ],
    include_dirs=["pyarcher/libraries/surface_operators/src"],
    extra_compile_args=["-fopenmp", "-O3"],
    extra_link_args=["-fopenmp"],
)
(coverage,) = cythonize(coverage, build_dir="build")

sbl = Extension(
    name="pymcac.tools.volume_surface.sbl_wrapper",
    sources=[
        "pymcac/tools/volume_surface/sbl_wrapper.pyx",
        "pymcac/tools/volume_surface/SBLVolumeSurface.cpp",
    ],
    include_dirs=["pymcac/tools/volume_surface", "ext_bin/sbl/include", "/opt/cgal/include"],
    library_dirs=["/opt/cgal/lib64"],
    extra_compile_args=["-fopenmp", "-O3", "-frounding-math", "-DNDEBUG"],
    libraries=["mpfr", "gmp"],
    extra_link_args=["-fopenmp"],
    language="c++",
)
(sbl,) = cythonize(sbl, build_dir="build")

setup(ext_modules=[coverage, sbl])
