TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += gauss_elim.cpp

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib -llapack -larmadillo -lblas
