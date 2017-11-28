#pragma once
#include <vector>
#include <string>
#include "planet.h"

class InitialCondition {
public:
    InitialCondition() {}
    static void twoBodyProblem();
    static void twoBodyRange();
    static void threeBodyProblem();
    static void tenBodyProblem();
    virtual void setupParticles(class SolarSystem& solarSystem) = 0;
    virtual std::string getName();
};
