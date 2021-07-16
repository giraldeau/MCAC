# ~~~
# - Try to find Xdmf
# Once done this will define
#  Xdmf_FOUND - System has Xdmf
#  Xdmf_INCLUDE_DIRS - The Xdmf include directories
#  Xdmf_LIBRARIES - The libraries needed to use Xdmf
#  Xdmf_DEFINITIONS - Compiler switches required for using Xdmf
# ~~~

find_package(PkgConfig QUIET)
pkg_check_modules(PC_Xdmf QUIET Xdmf)
set(Xdmf_DEFINITIONS ${PC_Xdmf_CFLAGS_OTHER})

find_path(
    Xdmf_INCLUDE_DIR
    NAMES XdmfDomain.hpp
    HINTS ${PC_Xdmf_INCLUDEDIR} ${PC_Xdmf_INCLUDE_DIRS}
    PATHS /usr/include
)

find_library(
    Xdmf_LIBRARY
    NAMES Xdmf
    HINTS ${PC_Xdmf_LIBDIR} ${PC_Xdmf_LIBRARY_DIRS}
    PATHS /usr/lib/x86_64-linux-gnu/xdmf/serial/ /usr/lib/x86_64-linux-gnu/xdmf/openmpi/
)

find_library(
    XdmfCore_LIBRARY
    NAMES XdmfCore HINT ${PC_Xdmf_LIBDIR} ${PC_Xdmf_LIBRARY_DIRS}
    PATHS /usr/lib/x86_64-linux-gnu/xdmf/serial/ /usr/lib/x86_64-linux-gnu/xdmf/openmpi/
)

set(Xdmf_INCLUDE_DIRS ${Xdmf_INCLUDE_DIR} ${PC_Xdmf_INCLUDE_DIRS})
set(Xdmf_LIBRARIES ${Xdmf_LIBRARY} ${XdmfCore_LIBRARY})

include(FindPackageHandleStandardArgs)
# ~~~
# handle the QUIETLY and REQUIRED arguments and set Xdmf_FOUND to TRUE
# if all listed variables are TRUE
# ~~~
find_package_handle_standard_args(Xdmf REQUIRED_VARS Xdmf_LIBRARY XdmfCore_LIBRARY Xdmf_INCLUDE_DIR)

if(Xdmf_FOUND AND NOT TARGET Xdmf::Xdmf)
    add_library(Xdmf::Xdmf UNKNOWN IMPORTED)
    set_target_properties(
        Xdmf::Xdmf PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Xdmf_INCLUDE_DIRS}"
    )
    set_target_properties(Xdmf::Xdmf PROPERTIES IMPORTED_LOCATION "${Xdmf_LIBRARY}")
endif()
if(Xdmf_FOUND AND NOT TARGET Xdmf::XdmfCore)
    add_library(Xdmf::XdmfCore UNKNOWN IMPORTED)
    set_target_properties(
        Xdmf::XdmfCore PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Xdmf_INCLUDE_DIRS}"
    )
    set_target_properties(Xdmf::XdmfCore PROPERTIES IMPORTED_LOCATION "${XdmfCore_LIBRARY}")
endif()

mark_as_advanced(Xdmf_INCLUDE_DIR Xdmf_LIBRARY XdmfCore_LIBRARY)
