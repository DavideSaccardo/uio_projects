#pragma once
#include "planet.h"
#include <string>

class Potential {
protected:
    double m_potentialEnergy = 0;

public:
    Potential() {}
    virtual void calculateForces(Planet* p1, Planet* p2, double& U) = 0;
    virtual std::string getName();
    void   resetPotentialEnergy() { m_potentialEnergy = 0; }
    double getPotentialEnergy()   { return m_potentialEnergy; }
    void setPotentialEnergy(double potentialEnergy);
};
