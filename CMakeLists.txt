cmake_minimum_required(VERSION 3.15)

project(MCAC-FULL)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

ExternalProject_Add(GSL_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/GSL
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/GSL
    GIT_REPOSITORY https://github.com/microsoft/GSL.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/GSL -DGSL_TEST=OFF
    EXCLUDE_FROM_ALL ON)

ExternalProject_Add(inipp_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/inipp
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/inipp
    GIT_REPOSITORY https://github.com/mcmtroffaes/inipp
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/inipp
    EXCLUDE_FROM_ALL ON)

ExternalProject_Add(sbl_project
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/sbl
    INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/sbl
    GIT_REPOSITORY https://gitlab.coria-cfd.fr/MCAC/sbl.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_SOURCE_DIR}/ext_bin/sbl
    EXCLUDE_FROM_ALL ON)

add_subdirectory(src)


# ADD python
