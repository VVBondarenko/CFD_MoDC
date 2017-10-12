TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    viscousvortexdomainsolver.cpp \
    streamline.cpp

LIBS += -L/usr/local/lib -lgsl -lgslcblas -fopenmp
QMAKE_LFLAGS += -fopenmp
QMAKE_CXXFLAGS += -fopenmp

HEADERS += \
    viscousvortexdomainsolver.h \
    streamline.h
