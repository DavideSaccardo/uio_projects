#include "integrator.h"
#include "solarsystem.h"

Integrator::Integrator(SolarSystem* solarSystem) {
    m_solarSystem = solarSystem;
}

void Integrator::setDt(double dt) {
    m_dt = dt;
}

std::string Integrator::getName() {
    return "Unknown";
}

