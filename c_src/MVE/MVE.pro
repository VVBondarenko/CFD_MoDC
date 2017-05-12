TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.c

LIBS += -L/usr/local/lib -lgsl -lgslcblas -lm -fopenmp
