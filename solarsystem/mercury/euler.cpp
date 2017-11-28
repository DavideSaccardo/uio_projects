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
    m_solarSystem->calculateForces();
    for (int i=0; i<planets.size(); i++) {

        planets.at(i)->setX(planets.at(i)->getX()+m_dt*planets.at(i)->getVx());
        planets.at(i)->setY(planets.at(i)->getY()+m_dt*planets.at(i)->getVy());
        planets.at(i)->setZ(planets.at(i)->getZ()+m_dt*planets.at(i)->getVz());
        planets.at(i)->setVx(planets.at(i)->getVx()+ m_dt*planets.at(i)->getFx()/planets.at(i)->getMass()-planets.at(0)->getVx());
        planets.at(i)->setVy(planets.at(i)->getVy()+ m_dt*planets.at(i)->getFy()/planets.at(i)->getMass()-planets.at(0)->getVy());
        planets.at(i)->setVz(planets.at(i)->getVz()+ m_dt*planets.at(i)->getFz()/planets.at(i)->getMass()-planets.at(0)->getVz());
                /*
         * This is where you need to update the positions and the velocities
         * and the velocities of each particle according to the Euler-Cromer
         * scheme.
         *
         * You can access the position of particle number i by
         *
         *      p->getPosition()
         *
         * Remember that the positions and velocities of all the particles are
         * vec3 vectors, which you can use the +=, -=, etc. operators on.
    */
           m_solarSystem->computeKineticEnergy();
        m_solarSystem->removeLinearMomentum();
    }
}

std::string Euler::getName() {
    return "Euler";
}
