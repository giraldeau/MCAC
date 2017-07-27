#-------------------------------------------------
#
# Project created by QtCreator 2017-06-15T16:25:02
#
#-------------------------------------------------

COMPILATOR = "INTEL"

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = DLCA
TEMPLATE = app

SOURCES += main.cpp\
        mainwindow.cpp\
    physical_model.cpp \
    Sphere_storage.cpp \
    Sphere_physics.cpp \
    aggregat.cpp

HEADERS += mainwindow.h\
        Sphere.h \
    physical_model.h \
    aggregat.h

FORMS += mainwindow.ui

RESOURCES += \
    mainwindow.qrc

CONFIG += warn_on debug_and_release debug_and_release_target build_all

release: DESTDIR = $$_PRO_FILE_PWD_/../build-Generateur_DLCA/release
debug:   DESTDIR = $$_PRO_FILE_PWD_/../build-Generateur_DLCA/debug

OBJECTS_DIR = $$DESTDIR/.obj
MOC_DIR = $$DESTDIR/.moc
RCC_DIR = $$DESTDIR/.qrc
UI_DIR = $$DESTDIR/.ui

equals(COMPILATOR, "INTEL"){
    QMAKE_CXXFLAGS += -g -traceback -xHost -qopenmp -static -std=c++11
    QMAKE_LFLAGS   += -g -traceback -qopenmp -std=c++11

    ### WARNINGS ###
    QMAKE_CXXFLAGS += -w3 -Wall


    QMAKE_CXXFLAGS_DEBUG += -fp-trap=common -ansi-alias-check -check=conversions,stack,uninit -fp-stack-check -fstack-security-check -ftrapuv -par-runtime-control=3 -vec-guard-write  -Wcheck
    #QMAKE_CXXFLAGS_DEBUG += -check-pointers=rw -check-pointers-dangling=all -check-pointers-narrowing -check-pointers-undimensioned
    #LIBS += -lchkp -lchkpwrap

    QMAKE_CXXFLAGS_RELEASE += -fast -ansi-alias -xHost -DNDEBUG -fp-model fast=2 -fPIC
    QMAKE_CXXFLAGS_RELEASE += -ip -ipo # -prof-gen:srcpos
    QMAKE_CXXFLAGS_RELEASE += -qopt-report=5
    QMAKE_LFLAGS_RELEASE   += -qopt-report=5
}

equals(COMPILATOR, "GNU"){

    QMAKE_CXXFLAGS += -g -march=native -fopenmp -static  -std=c++11
    QMAKE_LFLAGS   += -g  -fopenmp  -std=c++11

    ### WARNINGS ###
    QMAKE_CXXFLAGS += -Wpedantic -Wall -Wextra -Wformat=2 -Wmissing-include-dirs -Wswitch-default -Wswitch-enum -Wuninitialized -Wstrict-overflow=5
    QMAKE_CXXFLAGS += -Warray-bounds -Wtrampolines -Wfloat-equal -Wundef -Wshadow -Wunsafe-loop-optimizations -Wcast-qual -Wcast-align -Wconversion
    QMAKE_CXXFLAGS += -Wuseless-cast -Wlogical-op -Wno-aggressive-loop-optimizations -Wnormalized=nfkc -Wpacked -Wredundant-decls -Winline -Winvalid-pch
    QMAKE_CXXFLAGS += -Wvector-operation-performance -Wdisabled-optimization -Wnoexcept  -fext-numeric-literals -Wstrict-null-sentinel -Wold-style-cast -Woverloaded-virtual -Wsign-promo
    QMAKE_CXXFLAGS += -Wzero-as-null-pointer-constant -Wmissing-declarations -Weffc++ -Wpadded -Waggregate-return
    # Not interesting
    #QMAKE_CXXFLAGS += -Wsuggest-attribute=pure -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wsuggest-attribute=format
    QMAKE_CXXFLAGS += -Wno-missing-declarations -Wno-effc++ -Wno-padded -Wno-aggregate-return -Wno-unused-variable -Wno-unused-result -Wno-unused-parameter
    # QT problems
    QMAKE_CXXFLAGS += -Wno-zero-as-null-pointer-constant -Wno-long-long -Wno-useless-cast -Wno-conversion -Wno-float-equal -Wno-packed -Wno-switch-default -Wno-strict-overflow


    #CHECKS
    QMAKE_CXXFLAGS_DEBUG += -fcheck-new -ggdb3 -gpubnames -g3
    #QMAKE_CXXFLAGS_DEBUG += -fsanitize=address
    #QMAKE_LFLAGS_DEBUG   += -fsanitize=address
    #QMAKE_CXXFLAGS_DEBUG += -fsanitize=leak
    #QMAKE_LFLAGS_DEBUG   += -fsanitize=leak

    # OPTIMISATION
    QMAKE_CXXFLAGS_RELEASE += -Ofast -DNDEBUG -fPIC -fomit-frame-pointer
    QMAKE_CXXFLAGS_RELEASE += -flto -fwhole-program
    QMAKE_LFLAGS_RELEASE   += -flto
}


equals(COMPILATOR, "CLANG"){

    QMAKE_CXXFLAGS += -g -march=native -static  -std=c++11
    QMAKE_LFLAGS   += -g -std=c++11

    ### WARNINGS ###
    QMAKE_CXXFLAGS += -Wpedantic -Wall -Wextra -Wformat=2 -Wmissing-include-dirs -Wswitch-default -Wswitch-enum -Wuninitialized -Wstrict-overflow=5
    QMAKE_CXXFLAGS += -Warray-bounds -Wtrampolines -Wfloat-equal -Wundef -Wshadow -Wunsafe-loop-optimizations -Wcast-qual -Wcast-align -Wconversion
    QMAKE_CXXFLAGS += -Wuseless-cast -Wlogical-op -Wno-aggressive-loop-optimizations -Wnormalized=nfkc -Wpacked -Wredundant-decls -Winline -Winvalid-pch
    QMAKE_CXXFLAGS += -Wvector-operation-performance -Wdisabled-optimization -Wnoexcept -Wstrict-null-sentinel -Wold-style-cast -Woverloaded-virtual -Wsign-promo
    QMAKE_CXXFLAGS += -Wzero-as-null-pointer-constant -Wmissing-declarations -Weffc++ -Wpadded -Waggregate-return
    # Not interesting
    #QMAKE_CXXFLAGS += -Wsuggest-attribute=pure -Wsuggest-attribute=const -Wsuggest-attribute=noreturn -Wsuggest-attribute=format
    QMAKE_CXXFLAGS += -Wno-missing-declarations -Wno-effc++ -Wno-padded -Wno-aggregate-return -Wno-unused-variable -Wno-unused-result -Wno-unused-parameter
    # QT problems
    QMAKE_CXXFLAGS += -Wno-zero-as-null-pointer-constant -Wno-long-long -Wno-useless-cast -Wno-conversion -Wno-float-equal -Wno-packed -Wno-switch-default -Wno-strict-overflow


    #CHECKS
    QMAKE_CXXFLAGS_DEBUG += -fcheck-new -ggdb3 -g3
    #QMAKE_CXXFLAGS_DEBUG += -fsanitize=address
    #QMAKE_LFLAGS_DEBUG   += -fsanitize=address
    #QMAKE_CXXFLAGS_DEBUG += -fsanitize=leak
    #QMAKE_LFLAGS_DEBUG   += -fsanitize=leak

    # OPTIMISATION
    QMAKE_CXXFLAGS_RELEASE += -Ofast -DNDEBUG -fPIC -fomit-frame-pointer
    #QMAKE_CXXFLAGS_RELEASE += -flto
    #QMAKE_LFLAGS_RELEASE   += -flto
}

#PROFILING
LIBS += -L/opt/local/gperftools/lib -lprofiler -ltcmalloc


