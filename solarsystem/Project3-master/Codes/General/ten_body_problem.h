#pragma once
#include "initialcondition.h"
#include "planet.h"
#include <vector>
#include <string>


class ten_body_problem : public InitialCondition {
public:
    ten_body_problem() {}
    void setupParticles(class SolarSystem& solarSystem);
    std::string getName();
};
