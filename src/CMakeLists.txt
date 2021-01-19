cmake_minimum_required(VERSION 3.10)
project(MCAC CXX)

option(WITH_SBL "Use SBL for aggregate volume measure" OFF)
option(WITH_ARVO "Use ARVO for aggregate volume measure" OFF)
option(WITH_IO "Write output" ON)
option(DEBUG_MCAC "Usefull for debugging MCAC" OFF)

set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

### definitions des dependences  ###
include_directories(../include)

FILE(GLOB_RECURSE SRC_FILES *.cpp)
if (NOT WITH_ARVO)
    LIST(FILTER SRC_FILES EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/arvo/.*" )
endif()
# TODO temporary fix for old cmake 
link_directories(${PROJECT_SOURCE_DIR}/../ext_bin/xdmf/lib)
link_directories(${PROJECT_SOURCE_DIR}/../ext_bin/xdmf/lib64)
# TODO end temporary fix for old cmake 
add_executable(${PROJECT_NAME} ${SRC_FILES})
    
### Definition des options de compilation ###

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../bin/${CMAKE_BUILD_TYPE}/)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)
if (WITH_IO)
    find_package(LibXml2 REQUIRED)
    find_package(HDF5 REQUIRED CXX)
    find_package(Xdmf QUIET)
    if (NOT Xdmf_FOUND)
        add_dependencies(${PROJECT_NAME} xdmf_project)
    endif ()
endif ()
if (WITH_SBL)
    find_package(CGAL QUIET)
    if (CGAL_FOUND)
        include(${CGAL_MODULES_DIR}/CGAL_Macros.cmake)
        cgal_setup_module_path()
    else ()
        add_dependencies(${PROJECT_NAME} cgal_project)
        list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/cgal)
    endif ()
    find_package(MPFR REQUIRED)
    find_package(GMP REQUIRED)
    find_package(Boost 1.60.0 REQUIRED)
    add_definitions(-DCGAL_VERSION_GEQ_4_10)
endif ()
find_package(OpenMP QUIET COMPONENTS CXX)
find_package(Unwind QUIET)
if (NOT Unwind_FOUND)
        add_dependencies(${PROJECT_NAME} libunwind_project)
endif()

message(STATUS "")
message(STATUS "=======================================================================================")
message(STATUS "")
site_name(hostname)
message(STATUS "Compilation de MCAC sur : " ${hostname})
message(STATUS "")


### Options  ###

message(STATUS "Flags:")

target_link_libraries(${PROJECT_NAME} stdc++fs)
target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
if (Unwind_FOUND)
	target_link_libraries(${PROJECT_NAME} unwind::unwind)
else()
	target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../ext_bin/libunwind/include)
	target_link_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../ext_bin/libunwind/lib)
	target_link_libraries(${PROJECT_NAME} unwind)
endif()
add_dependencies(${PROJECT_NAME} GSL_project)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../ext_bin/GSL/include)
add_dependencies(${PROJECT_NAME} inipp_project)
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../ext_bin/inipp/include)

if (DEBUG_MCAC)
    add_definitions(-DDEBUG_MCAC)
endif()
if (WITH_SBL)
    message(STATUS "    With SBL")
    add_definitions(-DWITH_SBL)
    add_dependencies(${PROJECT_NAME} sbl_project)
    target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../ext_bin/sbl/include)
    if (NOT CGAL_FOUND)
        target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../ext_bin/cgal/include)
        target_link_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../ext_bin/cgal/lib64/)
    endif ()
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${MPFR_LIBRARIES} ${GMP_LIBRARIES})
else ()
    message(STATUS "    Without SBL")
endif ()
if (WITH_ARVO)
    message(STATUS "    With Arvo")
    add_definitions(-DWITH_SBL)
else ()
    message(STATUS "    Without Arvo")
endif ()

if (WITH_IO)
    message(STATUS "    With I/O")
    add_definitions(-DWITH_HDF5)
    target_include_directories(${PROJECT_NAME} PRIVATE ${LIBXML2_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${LIBXML2_LIBRARIES})

    target_include_directories(${PROJECT_NAME} PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES})

    if (Xdmf_FOUND)
        target_link_libraries(${PROJECT_NAME} Xdmf::Xdmf Xdmf::XdmfCore)
    else ()
        target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../ext_bin/xdmf/include)
        #target_link_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../ext_bin/xdmf/lib)
        #target_link_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/../ext_bin/xdmf/lib64)
        target_link_libraries(${PROJECT_NAME} XdmfCore Xdmf)
    endif ()

else ()
    message(STATUS "    Without I/O")
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES "Release")
    message(STATUS "    Optimized  build")
    #  set(CMAKE_EXE_LINK_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
elseif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    message(STATUS "    Debugging  build")
    #  set(CMAKE_EXE_LINK_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
endif ()

message(STATUS "${CMAKE_COMPILE_OPTIONS}")
message(STATUS ${CMAKE_CXX_COMPILER_ID} " compiler: "
        ${CMAKE_CXX_COMPILER} " (version: " ${CMAKE_CXX_COMPILER_VERSION} ")")

### Definition des flags de compilation suivant le compilateur ###

# icpc compiler (LINUX)
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "-g -traceback -xHost")
    set(CMAKE_CXX_FLAGS_RELEASE " -ansi-alias -xHost -fp-model fast=2  ") #   -no-prec-div -fp-model fast=2 -init=zero -prof-use -prof-dir=prof_dir/  -prof-gen
    set(CMAKE_CXX_FLAGS_DEBUG "-fp-trap=common -ansi-alias-check -check=conversions,stack,uninit -fp-stack-check -fstack-security-check -ftrapuv -par-runtime-control=3 -vec-guard-write  -Wcheck -O0")

    # g++ compiler (LINUX)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-g -march=native -frounding-math")
    set(CMAKE_CXX_FLAGS_RELEASE " -O2")  # O3 and Ofast does not work with sbl
    set(CMAKE_CXX_FLAGS_DEBUG "-fcheck-new -ggdb3 -gpubnames -g3  -O0 \
      -Wpedantic -Wall -Wextra -Wformat=2 -Wmissing-include-dirs -Wswitch-default -Wswitch-enum -Wuninitialized -Wstrict-overflow=5 \
      -Warray-bounds -Wtrampolines -Wfloat-equal -Wundef -Wshadow -Wunsafe-loop-optimizations -Wcast-qual -Wcast-align -Wconversion \
      -Wuseless-cast -Wlogical-op -Wno-aggressive-loop-optimizations -Wnormalized=nfkc -Wpacked -Wredundant-decls -Winline -Winvalid-pch \
      -Wvector-operation-performance -Wdisabled-optimization -Wnoexcept  -fext-numeric-literals -Wstrict-null-sentinel -Wold-style-cast -Woverloaded-virtual -Wsign-promo \
      -Wzero-as-null-pointer-constant -Wmissing-declarations -Weffc++ -Wpadded -Waggregate-return \
      -Wsuggest-attribute=pure -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wsuggest-attribute=format \
      -Weffc++ -Wunused-variable -Wunused-result -Wunused-parameter\
      -Wzero-as-null-pointer-constant -Wlong-long -Wuseless-cast -Wconversion -Wfloat-equal -Wpacked -Wswitch-default -Wstrict-overflow \
      -Wmissing-declarations \
      -Wno-aggregate-return  -Wno-padded")
    # Problems in SBL
#    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-catch-value -Wno-unused-parameter -Wno-effc++ \
#                                -Wno-zero-as-null-pointer-constant -Wno-conversion -Wno-old-style-cast\
#                                -Wno-switch-enum -Wno-type-limits -Wno-unused-but-set-variable -Wno-pedantic")

elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "-g -march=native")
    set(CMAKE_CXX_FLAGS_RELEASE " -Ofast") #
    set(CMAKE_CXX_FLAGS_DEBUG "-fcheck-new -ggdb3 -gpubnames -g3  -O0  -ftrapv \
                                -Wall -Wextra -Wpedantic -Weverything \
                                -Wno-c++98-compat -Wno-padded -Wno-global-constructors -Wno-unused-macros \
                                -Wno-unreachable-code -Wno-header-hygiene -Wno-c++98-compat-pedantic -Wno-weak-vtables \
                                -Wno-exit-time-destructors")
    # Problems in SBL
#    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-reserved-id-macro -Wno-documentation-unknown-command \
#                                 -Wno-extra-semi -Wno-documentation -Wno-sign-conversion -Wno-old-style-cast\
#                                 -Wno-unused-parameter -Wno-zero-as-null-pointer-constant -Wno-switch-enum \
#                                 -Wno-shorten-64-to-32 -Wno-deprecated")
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES "Release")
    message(STATUS "${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
elseif (${CMAKE_BUILD_TYPE} MATCHES "Debug")
    message(STATUS "${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
endif ()
message(STATUS "")
message(STATUS "=======================================================================================")
message(STATUS "")

set_target_properties(MCAC PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(MCAC PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)

