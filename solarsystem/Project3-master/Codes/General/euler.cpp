#include "euler.h"
#include "planet.h"
#include "solarsystem.h"
#include <iostream>
#include <cmath>
using namespace std;
Euler::Euler(SolarSystem* solarSystem)
    : Integrator(solarSystem) {
}

void Euler::integrateOneStep(std::vector<Planet*> planets) {
    //Here we use the Euler method to compute positions and velocities and the planets
    m_solarSystem->calculateForces(); //before Euler we compute the forces acting on each planet
    for (int i=0; i<planets.size(); i++) {
        //computation of positions and velocities (next time we are going to use the class vec3 to do all this with vectors, sorry)
        planets.at(i)->setX(planets.at(i)->getX()+m_dt*planets.at(i)->getVx());
        planets.at(i)->setY(planets.at(i)->getY()+m_dt*planets.at(i)->getVy());
        planets.at(i)->setZ(planets.at(i)->getZ()+m_dt*planets.at(i)->getVz());
        planets.at(i)->setVx(planets.at(i)->getVx()+ m_dt*planets.at(i)->getFx()/planets.at(i)->getMass()-planets.at(0)->getVx());
        planets.at(i)->setVy(planets.at(i)->getVy()+ m_dt*planets.at(i)->getFy()/planets.at(i)->getMass()-planets.at(0)->getVy());
        planets.at(i)->setVz(planets.at(i)->getVz()+ m_dt*planets.at(i)->getFz()/planets.at(i)->getMass()-planets.at(0)->getVz());

        m_solarSystem->computeKineticEnergy(); //for each planet it computes the kinetic energy
    }
}

std::string Euler::getName() {
    return "Euler";
}
