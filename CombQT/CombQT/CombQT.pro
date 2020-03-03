TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
# remove possible other optimization flags
QMAKE_CXXFLAGS_RELEASE -= -O
QMAKE_CXXFLAGS_RELEASE -= -O1
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE -= -O3
# add the desired -O3 if not present
QMAKE_CXXFLAGS_RELEASE += -O3

LIBS += -lfftw3

SOURCES += \
        analysisbank.cpp \
        generator.cpp \
        main.cpp

HEADERS += \
    analysisbank.h
