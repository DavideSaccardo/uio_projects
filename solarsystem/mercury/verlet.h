#pragma once
#include "integrator.h"
#include "planet.h"
#include <string>

class Verlet : public Integrator {
private:
    bool m_firstStep = true;

public:
    Verlet(class SolarSystem* solarSystem);
    std::string getName();
    void integrateOneStep(std::vector<Planet*> planets);
};

