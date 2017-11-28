#include "two_body_problem.h"
#include "solarsystem.h"
#include <cmath>

void two_body_problem::setupParticles(SolarSystem &solarSystem) {

    // all the mass should be given as ratio of the sun's mass
    //          mass   x  y  z  vx vy vz
    // velocity in AU/day

    double vfactor = 1.0; /* change earth speed */

    Planet* sun= new Planet (1,0,0,0,0,0,0); //here the Sun is fixed at the center of frame
    Planet* earth= new Planet (3e-6, 1, 0,0, 0, vfactor*2*M_PI/(365), 0 );

    solarSystem.add(sun); // has to be the first element
    solarSystem.add(earth);

}

std::string two_body_problem::getName() {
    return "Two-body";
}
