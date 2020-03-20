cmake_minimum_required(VERSION 3.13)

project(MCAC-FULL)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

if (NOT SBL_GIT_REP)
    set(SBL_GIT_REP https://gitlab.coria-cfd.fr/MCAC/sbl.git)
endif ()
option(WITH_PYTHON "Compile PyMCAC" ON)

ExternalProject_Add(GSL_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/GSL
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/GSL
    GIT_REPOSITORY https://github.com/microsoft/GSL.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/GSL -DGSL_TEST=OFF
    EXCLUDE_FROM_ALL ON
)

ExternalProject_Add(inipp_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/inipp
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/inipp
    GIT_REPOSITORY https://gitlab.coria-cfd.fr/MCAC/inipp.git
    GIT_TAG no_tests
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/inipp
    EXCLUDE_FROM_ALL ON
)

ExternalProject_Add(sbl_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/sbl
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/sbl
    GIT_REPOSITORY ${SBL_GIT_REP}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/sbl
    EXCLUDE_FROM_ALL ON
)

ExternalProject_Add(xdmf_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/xdmf
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/xdmf
    GIT_REPOSITORY https://gitlab.kitware.com/xdmf/xdmf.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/xdmf
    EXCLUDE_FROM_ALL ON
)

ExternalProject_Add(cgal_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/cgal
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/cgal
    GIT_TAG releases/CGAL-4.9-branch-tip
    GIT_REPOSITORY https://github.com/CGAL/cgal.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/cgal -DWITH_CGAL_Qt5=OFF -DWITH_CGAL_ImageIO=OFF
    EXCLUDE_FROM_ALL ON
)


add_subdirectory(src)

if (WITH_PYTHON)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/PyMCAC.cmake)
endif ()

