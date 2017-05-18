TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CFLAGS_DEBUG += -fopenmp

SOURCES += main.c



LIBS += -L/usr/local/lib -lgsl -lgslcblas -lm -fopenmp -Wall -ffast-math
