#-------------------------------------------------
#
# Project created by QtCreator 2017-06-15T16:25:02
#
#-------------------------------------------------

COMPILATOR = "GNU" # GNU or INTEL or CLANG
WITH_IO = "1"           # 1 or 0
WITH_SBL = "0"          # 1 or 0
STATIC = "0"            # 0
PROFILING = "0"         # 0 or 1

TARGET = DLCA
TEMPLATE = app
QT = ""

SOURCES +=\
        storage.hpp \
        aggregat.hpp \
        aggregatList.hpp \
        IO.hpp \
        calcul.hpp \
        physical_model.hpp \
        Sphere.hpp \
        Spherelist.hpp \
        statistics.hpp \
        storagelist.hpp \
        verlet.hpp \
        sblvolumesurface.hpp \
        verlet.cpp \
        physical_model.cpp \
        Sphere_storage.cpp \
        Spherelist_storage.cpp \
        Sphere_physics.cpp \
        Spherelist_physics.cpp \
        aggregat.cpp \
        aggregatList.cpp \
        calcul.cpp\
        main.cpp \
        IO.cpp \
        statistics.cpp \
        sblvolumesurface.cpp

#SOURCES +=\
#        verlet.cpp \
#        physical_model.cpp \
#        Sphere_storage.cpp \
#        Spherelist_storage.cpp \
#        Sphere_physics.cpp \
#        Spherelist_physics.cpp \
#        aggregat.cpp \
#        aggregatList.cpp \
#        mainwindow.cpp\
#        main.cpp \
#        sblvolumesurface.cpp

#HEADERS += \
#        mainwindow.hpp \
#        Sphere.hpp \
#        Spherelist.hpp \
#        physical_model.hpp \
#        aggregat.hpp \
#        aggregatList.hpp \
#        verlet.hpp \
#        storage.hpp \
#        storagelist.hpp
#        sblvolumesurface.hpp \


LIBS += -lstdc++fs

CONFIG += warn_on debug_and_release debug_and_release_target build_all

release: DESTDIR = $$_PRO_FILE_PWD_/../build-Generateur_DLCA/release
debug:   DESTDIR = $$_PRO_FILE_PWD_/../build-Generateur_DLCA/debug

equals(WITH_SBL, "1"){
    INCLUDEPATH +=  /home/pouxa/packages/sbl/install/include
    LIBS += -lCGAL -lmpfr -lgmp
}

equals(WITH_IO, "1"){

    DEFINES += WITH_HDF5
    INCLUDEPATH += /usr/include/libxml2/

    equals(STATIC, "0"){
        LIBS += -lXdmf -lXdmfCore
        LIBS += -lhdf5
    } else {
        INCLUDEPATH += /opt/local/xdmf/include/
        LIBS += /opt/local/xdmf/lib/libXdmf.a /opt/local/xdmf/lib/libXdmfCore.a
        LIBS += /usr/lib64/libhdf5.a  -lz -lsz
    }

    LIBS += -lxml2 -ltiff -ldl
}

equals(PROFILING, "1"){
    LIBS += -lprofiler -ltcmalloc
}

equals(COMPILATOR, "INTEL"){
    QMAKE_CXXFLAGS += -std=c++14 -g -traceback -xHost -qopenmp -static -fPIC -DNDEBUG
    QMAKE_LFLAGS   += -std=c++14 -g -traceback -xHost -qopenmp -DNDEBUG

    ### WARNINGS ###
    QMAKE_CXXFLAGS += -Wall -w3 -diag-enable=3 -Wremarks -Wtrigraphs -Wcomment -Wdeprecated -Wno-effc++ -Wextra-tokens -Wformat -Wformat-security -Wic-pointer -Winline -Wmain
    QMAKE_CXXFLAGS += -Wnon-virtual-dtor -Wpointer-arith -Wreorder -Wreturn-type -Wshadow -Wsign-compare -Wuninitialized -Wunknown-pragmas -Wunused-function
    QMAKE_CXXFLAGS += -Wwrite-strings

    QMAKE_CXXFLAGS_DEBUG += -fp-trap=common -ansi-alias-check -check=conversions,stack,uninit -fp-stack-check -fstack-security-check -ftrapuv -par-runtime-control=3 -vec-guard-write  -Wcheck -O0
    #QMAKE_CXXFLAGS_DEBUG += -check-pointers=rw -check-pointers-dangling=all -check-pointers-narrowing -check-pointers-undimensioned
    #LIBS += -lchkp -lchkpwrap

    QMAKE_CXXFLAGS_RELEASE += -fast -ansi-alias -DNDEBUG -fp-model fast=2
    QMAKE_CXXFLAGS_RELEASE += -ip -ipo # -prof-gen:srcpos
    QMAKE_CXXFLAGS_RELEASE += -qopt-report=5
    QMAKE_LFLAGS_RELEASE   += -qopt-report=5
}

equals(COMPILATOR, "GNU"){

    QMAKE_CXXFLAGS += -g -march=native -fopenmp  -std=c++14 -frounding-math
    QMAKE_LFLAGS   += -g  -fopenmp  -std=c++14

    ### WARNINGS ###
    QMAKE_CXXFLAGS += -Wpedantic -Wall -Wextra -Wformat=2 -Wmissing-include-dirs -Wswitch-default -Wswitch-enum -Wuninitialized -Wstrict-overflow=5
    QMAKE_CXXFLAGS += -Warray-bounds -Wtrampolines -Wfloat-equal -Wundef -Wshadow -Wunsafe-loop-optimizations -Wcast-qual -Wcast-align -Wconversion
    QMAKE_CXXFLAGS += -Wuseless-cast -Wlogical-op -Wno-aggressive-loop-optimizations -Wnormalized=nfkc -Wpacked -Wredundant-decls -Winline -Winvalid-pch
    QMAKE_CXXFLAGS += -Wvector-operation-performance -Wdisabled-optimization -Wnoexcept  -fext-numeric-literals -Wstrict-null-sentinel -Wold-style-cast -Woverloaded-virtual -Wsign-promo
    QMAKE_CXXFLAGS += -Wzero-as-null-pointer-constant -Wmissing-declarations -Weffc++ -Wpadded -Waggregate-return
    QMAKE_CXXFLAGS += -Wsuggest-attribute=pure -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wsuggest-attribute=format
    QMAKE_CXXFLAGS += -Weffc++ -Wunused-variable -Wunused-result -Wunused-parameter
    # QT problems
    QMAKE_CXXFLAGS += -Wzero-as-null-pointer-constant -Wlong-long -Wuseless-cast -Wconversion -Wfloat-equal -Wpacked -Wswitch-default -Wstrict-overflow
    QMAKE_CXXFLAGS += -Wmissing-declarations
    # Not a problems
    QMAKE_CXXFLAGS += -Wno-aggregate-return  -Wno-padded

    #CHECKS
    QMAKE_CXXFLAGS_DEBUG += -fcheck-new -ggdb3 -gpubnames -g3  -O0
    #QMAKE_CXXFLAGS_DEBUG += -fsanitize=address
    #QMAKE_LFLAGS_DEBUG   += -fsanitize=address
    #QMAKE_CXXFLAGS_DEBUG += -fsanitize=leak
    #QMAKE_LFLAGS_DEBUG   += -fsanitize=leak

    # OPTIMISATION
    QMAKE_CXXFLAGS_RELEASE += -Ofast -DNDEBUG -fPIC
    QMAKE_CXXFLAGS_RELEASE += -flto -fwhole-program
    QMAKE_LFLAGS_RELEASE   += -flto

    equals(STATIC, "1"){
        QMAKE_CXXFLAGS += -Wl,--no-allow-shlib-undefined,--no-undefined,--as-needed -static-libgcc -static-libstdc++
        QMAKE_LFLAGS   += -Wl,--no-allow-shlib-undefined,--no-undefined,--as-needed -static-libgcc -static-libstdc++
    }

}


equals(COMPILATOR, "CLANG"){

    QMAKE_CXXFLAGS += -g -march=native -static  -std=c++14 -fopenmp
    QMAKE_LFLAGS   += -g -std=c++14

    ### WARNINGS ###
    QMAKE_CXXFLAGS += -Wpedantic -Wall -Wextra -Wformat=2 -Wmissing-include-dirs -Wswitch-default -Wswitch-enum -Wuninitialized -Wstrict-overflow=5
    QMAKE_CXXFLAGS += -Warray-bounds -Wfloat-equal -Wundef -Wshadow -Wcast-qual -Wcast-align -Wconversion
    QMAKE_CXXFLAGS += -Wpacked -Wredundant-decls -Winline -Winvalid-pch
    QMAKE_CXXFLAGS += -Wdisabled-optimization -Wold-style-cast -Woverloaded-virtual -Wsign-promo
    QMAKE_CXXFLAGS += -Wzero-as-null-pointer-constant -Wmissing-declarations -Weffc++ -Wpadded -Waggregate-return
    QMAKE_CXXFLAGS += -Weffc++ -Wpadded -Wunused-variable -Wunused-result -Wunused-parameter
    # QT problems
    QMAKE_CXXFLAGS += -Wzero-as-null-pointer-constant -Wlong-long -Wconversion -Wfloat-equal -Wpacked -Wswitch-default -Wstrict-overflow
    QMAKE_CXXFLAGS += -Wmissing-declarations
    # Not a problems
    QMAKE_CXXFLAGS += -Wno-aggregate-return -Wno-padded



    #CHECKS
    QMAKE_CXXFLAGS_DEBUG += -fcheck-new -ggdb3 -g3 -O0 -ftrapv
    #QMAKE_CXXFLAGS_DEBUG += -fsanitize=address
    #QMAKE_LFLAGS_DEBUG   += -fsanitize=address
    #QMAKE_CXXFLAGS_DEBUG += -fsanitize=leak
    #QMAKE_LFLAGS_DEBUG   += -fsanitize=leak

    # OPTIMISATION
    QMAKE_CXXFLAGS_RELEASE += -Ofast -DNDEBUG -fPIC -fomit-frame-pointer
    #QMAKE_CXXFLAGS_RELEASE += -flto
    #QMAKE_LFLAGS_RELEASE   += -flto
}

