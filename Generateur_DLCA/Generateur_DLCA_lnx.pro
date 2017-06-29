#-------------------------------------------------
#
# Project created by QtCreator 2017-06-15T16:25:02
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = DLCA
TEMPLATE = app

SOURCES += main.cpp\
        mainwindow.cpp\
        Sphere.cpp

HEADERS += mainwindow.h\
        Sphere.h

FORMS += mainwindow.ui

RESOURCES += \
    mainwindow.qrc

#CONFIG += warn_on

#### INTEL ###
#QMAKE_CXXFLAGS_DEBUG += -fp-trap=common -ansi-alias-check -check=conversions,stack,uninit -fp-stack-check -fstack-security-check -ftrapuv -par-runtime-control=3 -vec-guard-write  -Wcheck
#QMAKE_CXXFLAGS_DEBUG += -check-pointers=rw -check-pointers-dangling=all -check-pointers-narrowing -check-pointers-undimensioned
#LIBS += -lchkp -lchkpwrap

QMAKE_CXXFLAGS_RELEASE = -Ofast -falign-functions=16 -ansi-alias -xHost -static -no-prec-div -DNDEBUG

### When not profiling
##QMAKE_CXXFLAGS_RELEASE += -lto

### When profiling
QMAKE_CXXFLAGS_RELEASE +=  -g -traceback #-fno-inline-functions -p -pg  -fno-omit-frame-pointer  -fno-optimize-sibling-calls   -fno-inline
QMAKE_LFLAGS_RELEASE += -g #-p -pg
LIBS += -L/opt/local/gperftools/lib -lprofiler -ltcmalloc





#### GNU ###
#QMAKE_CXXFLAGS_DEBUG += -fp-trap=common -ansi-alias-check -check=conversions,stack,uninit -fp-stack-check -fstack-security-check -ftrapuv -par-runtime-control=3 -vec-guard-write  -Wcheck
#QMAKE_CXXFLAGS_DEBUG += -check-pointers=rw -check-pointers-dangling=all -check-pointers-narrowing -check-pointers-undimensioned
#LIBS += -lchkp -lchkpwrap

#QMAKE_CXXFLAGS_RELEASE = -Ofast -march=native -falign-functions=16 -static -ffast-math -DNDEBUG -Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-unused-result
### When not profiling
##QMAKE_CXXFLAGS_RELEASE += -lto

### When profiling
#QMAKE_CXXFLAGS_RELEASE +=  -g -traceback #-fno-inline-functions -p -pg  -fno-omit-frame-pointer  -fno-optimize-sibling-calls   -fno-inline
#QMAKE_LFLAGS_RELEASE += -g #-p -pg
#LIBS += -L/opt/local/gperftools/lib -lprofiler -ltcmalloc
