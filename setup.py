#!/usr/bin/env python3
# coding=utf-8
"""Compile and install pymcac"""
import sys
from subprocess import CalledProcessError, run

from setuptools import Extension as Extension
from setuptools import setup as setup

setup()

# noinspection PyPep8
import numpy as np

coverage = Extension(
    name="pymcac.tools.coverage.coverages_cython",
    sources=["pymcac/tools/coverage/coverages_cython.pyx"],
    include_dirs=[np.get_include()],
    extra_compile_args=["-fopenmp", "-O3"],
    extra_link_args=["-fopenmp"],
)
# noinspection PyUnusedName
coverage.cython_c_in_temp = True

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
# noinspection PyUnusedName
sbl.cython_c_in_temp = True

try:
    # Build the cython extension
    # --------------------------
    setup(ext_modules=[coverage])
except BaseException as e:
    print("\n################### ERROR ##################################", file=sys.stderr)
    print("Unable to compile the coverage library", file=sys.stderr)
    print(e, file=sys.stderr)
    print("################### ERROR ##################################\n", file=sys.stderr)

try:
    # Build the cython extension
    # --------------------------
    setup(ext_modules=[sbl])
except BaseException as e:
    print("\n################### ERROR ##################################", file=sys.stderr)
    print("Unable to compile the sbl wrapper", file=sys.stderr)
    print(e, file=sys.stderr)
    print("################### ERROR ##################################\n", file=sys.stderr)

try:
    run(["git", "branch"], check=True)
except CalledProcessError:
    setup()
else:
    # Add the version from git
    # --------------------------------
    setup(use_scm_version={"write_to": "pymcac/version.py", "fallback_version": "UNKNOWN"})
