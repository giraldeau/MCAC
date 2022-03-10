if(NOT CGAL_GENERATOR_SPECIFIC_SETTINGS_FILE_INCLUDED)
    set(CGAL_GENERATOR_SPECIFIC_SETTINGS_FILE_INCLUDED 1)

    message(STATUS "Targeting ${CMAKE_GENERATOR}")

    if(MSVC)
        message(STATUS "Target build environment supports auto-linking")
        set(CGAL_AUTO_LINK_ENABLED TRUE)
    endif()

    if(MSVC15)
        set(CGAL_TOOLSET "vc150")
        message(STATUS "Using VC15 compiler.")
    elseif(MSVC14)
        set(CGAL_TOOLSET "vc140")
        message(STATUS "Using VC14 compiler.")
    elseif(MSVC12)
        set(CGAL_TOOLSET "vc120")
        message(STATUS "Using VC12 compiler.")
    elseif(MSVC11)
        set(CGAL_TOOLSET "vc110")
        message(STATUS "Using VC11 compiler.")
    elseif(MSVC10)
        set(CGAL_TOOLSET "vc100")
        message(STATUS "Using VC10 compiler.")
    elseif(MSVC90)
        set(CGAL_TOOLSET "vc90")
        message(STATUS "Using VC90 compiler.")
    elseif(MSVC80)
        set(CGAL_TOOLSET "vc80")
        message(STATUS "Using VC80 compiler.")
    elseif(MSVC71)
        set(CGAL_TOOLSET "vc71")
        message(STATUS "Using VC71 compiler.")
    else()
        message(STATUS "Using ${CMAKE_CXX_COMPILER} compiler.")
    endif()

    # From james Bigler, in the cmake users list.
    if(APPLE)
        exec_program(
            uname ARGS
            -v
            OUTPUT_VARIABLE DARWIN_VERSION
        )
        string(REGEX MATCH "[0-9]+" DARWIN_VERSION ${DARWIN_VERSION})
        message(STATUS "DARWIN_VERSION=${DARWIN_VERSION}")
        if(DARWIN_VERSION GREATER 8)
            message(STATUS "Mac Leopard detected")
            set(CGAL_APPLE_LEOPARD 1)
        endif()
    endif()

    if(NOT "${CMAKE_CFG_INTDIR}" STREQUAL ".")
        set(HAS_CFG_INTDIR
            TRUE
            CACHE INTERNAL "Generator uses intermediate configuration directory"
        )
        message(STATUS "Generator uses intermediate configuration directory: ${CMAKE_CFG_INTDIR}")
    endif()

endif()
