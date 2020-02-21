TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

#LIBS += -lfftw3
SOURCES += \
        ../CombQT/CombQT/main.cpp
       # main.cpp



#INCLUDEPATH += D:/Qt/fftw3/
#LIBS +=  D:/Qt/fftw3/libfftw3-3.a\
#         D:/Qt/fftw3/libfftw3f-3.a\
#        D:/Qt/fftw3/libfftw3l-3.a


win32: LIBS += -L$$PWD/fftw3/ -lfftw3-3

INCLUDEPATH += $$PWD/fftw3
DEPENDPATH += $$PWD/fftw3

win32: LIBS += -L$$PWD/fftw3/ -lfftw3f-3

INCLUDEPATH += $$PWD/fftw3
DEPENDPATH += $$PWD/fftw3

win32: LIBS += -L$$PWD/fftw3/ -lfftw3l-3

INCLUDEPATH += $$PWD/fftw3
DEPENDPATH += $$PWD/fftw3
