#include "potential.h"

void Potential::setPotentialEnergy(double potentialEnergy)
{
    m_potentialEnergy = potentialEnergy;
}

std::string Potential::getName() {
    return "Unknown";
}

