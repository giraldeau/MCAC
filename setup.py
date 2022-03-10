#!/usr/bin/env python3
# Copyright (C) 2021 CORIA
"""Compile and install pymcac."""
from pathlib import Path

from Cython.Build import cythonize
from numpy.distutils.core import Extension, setup

ext_modules = []
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
ext_modules.append(coverage)

if Path("ext_bin/sbl/include").exists():
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
    ext_modules.append(sbl)

setup(ext_modules=ext_modules)
