cmake_minimum_required(VERSION 3.10)

project(MCAC-FULL)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

if(NOT SBL_GIT_REP)
    set(SBL_GIT_REP https://gitlab.coria-cfd.fr/MCAC/sbl.git)
endif()
option(WITH_PYTHON "Compile PyMCAC" ON)

# Create a list with an alternate separator
string(REPLACE ";" "|" CMAKE_PREFIX_PATH_ALT_SEP "${CMAKE_PREFIX_PATH}")

externalproject_add(
    GSL_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/GSL
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/GSL
    GIT_REPOSITORY https://github.com/microsoft/GSL.git
    GIT_TAG main
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/GSL -DGSL_TEST=OFF
               -DCMAKE_BUILD_TYPE=Release
    EXCLUDE_FROM_ALL ON
)

externalproject_add(
    inipp_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/inipp
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/inipp
    GIT_REPOSITORY https://gitlab.coria-cfd.fr/MCAC/inipp.git
    GIT_TAG no_tests
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/inipp
               -DCMAKE_BUILD_TYPE=Release
    EXCLUDE_FROM_ALL ON
)

externalproject_add(
    sbl_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/sbl
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/sbl
    GIT_REPOSITORY ${SBL_GIT_REP}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/sbl
               -DCMAKE_BUILD_TYPE=Release
    EXCLUDE_FROM_ALL ON
)

externalproject_add(
    xdmf_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/xdmf
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/xdmf
    GIT_REPOSITORY https://gitlab.kitware.com/xdmf/xdmf.git
    LIST_SEPARATOR | # Use the alternate list separator
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/xdmf
               -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH_ALT_SEP} -DBUILD_SHARED_LIBS=ON
               -DCMAKE_BUILD_TYPE=Release
    EXCLUDE_FROM_ALL ON
)

externalproject_add(
    cgal_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/cgal
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/cgal
    GIT_TAG v5.2
    GIT_REPOSITORY https://github.com/CGAL/cgal.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/cgal -DWITH_CGAL_Qt5=OFF
               -DWITH_CGAL_ImageIO=OFF -DCMAKE_BUILD_TYPE=Release
    EXCLUDE_FROM_ALL ON
)

externalproject_add(
    libunwind_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/libunwind
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/libunwind
    GIT_TAG v1.5
    GIT_REPOSITORY https://github.com/libunwind/libunwind.git
    CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ext/libunwind/autogen.sh
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/ext/libunwind/configure
            --prefix=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/libunwind
    EXCLUDE_FROM_ALL ON
)

add_subdirectory(src)

if(WITH_PYTHON)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/PyMCAC.cmake)
endif()
