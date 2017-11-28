#pragma once
#include <vector>
#include <string>
#include "planet.h"

class InitialCondition {
public:
    InitialCondition() {}
    static void twoBodyProblem();
    //static void threeBodyProblem();
    //static void fourBodyProblem();
    virtual void setupParticles(class SolarSystem& solarSystem) = 0;
    virtual std::string getName();
};
