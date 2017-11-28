#pragma once
#include "planet.h"
#include "potential.h"
#include <string>

class NewtonianGravity : public Potential {
private:
    double m_G;
public:
    NewtonianGravity(double G);
    void calculateForces(Planet* p1, Planet* p2, double& U);
    std::string getName();
};
