#pragma once
#include "initialcondition.h"
#include "planet.h"
#include <vector>
#include <string>


class two_body_problem : public InitialCondition {
public:
    two_body_problem() {}
    void setupParticles(class SolarSystem& solarSystem);
    std::string getName();
};

