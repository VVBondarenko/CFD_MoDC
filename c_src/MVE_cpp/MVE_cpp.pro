TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    viscousvortexdomainsolver.cpp \
    streamline.cpp

LIBS += -L/usr/local/lib -lgsl -lgslcblas -parallel -qopenmp -tbb -simd -O3#-fopenmp
#QMAKE_LFLAGS += -fopenmp
QMAKE_CXXFLAGS += -qopenmp -tbb -simd -O3
#-march=i686

HEADERS += \
    viscousvortexdomainsolver.h \
    streamline.h
