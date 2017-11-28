TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    euler.cpp \
    initialcondition.cpp \
    integrator.cpp \
    newtoniangravity.cpp \
    planet.cpp \
    potential.cpp \
    solarsystem.cpp \
    ten_body_problem.cpp \
    three_body_problem.cpp \
    two_body_problem.cpp \
    verlet.cpp \
    main.cpp

HEADERS += \
    euler.h \
    initialcondition.h \
    integrator.h \
    newtoniangravity.h \
    planet.h \
    potential.h \
    solarsystem.h \
    ten_body_problem.h \
    three_body_problem.h \
    two_body_problem.h \
    verlet.h \
    catch.hpp
