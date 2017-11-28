#pragma once
#include "integrator.h"
#include "planet.h"
#include <vector>

class Euler : public Integrator {
private:
    bool m_firstStep = true;
public:
    Euler(class SolarSystem* system);
    void integrateOneStep(std::vector<Planet*> planets);
    std::string getName();
};
