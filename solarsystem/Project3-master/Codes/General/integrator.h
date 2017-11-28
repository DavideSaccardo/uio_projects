#pragma once
#include <vector>
#include <string>
#include "planet.h"

class Integrator {
protected:
    double               m_dt             = 1e-3; //Here you canc hange the step size
    class SolarSystem*   m_solarSystem    = nullptr;

public:
    Integrator(class SolarSystem* solarSystem);
    void setDt(double dt);
    double getDt() { return m_dt; }
    virtual std::string getName();
    virtual void integrateOneStep(std::vector<Planet*> planets) = 0;
};
