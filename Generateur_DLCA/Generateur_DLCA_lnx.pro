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

#QMAKE_CXXFLAGS_RELEASE = -Ofast -falign-functions=16 -ansi-alias -fstrict-aliasing -xHost -static -no-prec-div
#QMAKE_CXXFLAGS_RELEASE += -lto
#QMAKE_CXXFLAGS_RELEASE +=  -g -traceback -fno-inline-functions -p -pg  -fno-omit-frame-pointer -DNDEBUG  -fno-inline-functions-called-once -fno-optimize-sibling-calls  -fno-default-inline    -fno-inline
#QMAKE_LFLAGS_RELEASE += -g -p -pg

