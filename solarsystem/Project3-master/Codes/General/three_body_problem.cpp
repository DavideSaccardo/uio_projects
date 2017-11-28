#include "three_body_problem.h"
#include "solarsystem.h"
#include <cmath>

void three_body_problem::setupParticles(SolarSystem &solarSystem) {


    // all the mass should be given as ratio of the sun's mass
    //          mass   x  y  z  vx vy vz
    // velocity in AU/day

    //Here we declare Sun, Earth, Jupiter.
    //If you want to fix the Sun you have to put its velocities equal to 0, i.e.
    //you have to delete vx, vy, vz of the sun and put 0,0,0, then you have to go to solarsystem->void SolarSystem::calculateForces()
    // and you have to uncomment " //m_planets.at(0)->resetForce(); " this keeps the Sun fixed at the center of our frame. Finally you have
    //to comment lines 118,119,120 and uncomment lines 108,109,110 in planet.cpp

    //As the system is configured with the tenbodyproblem, to have the sun not fixed (and CM fixed) here you need to uncomment lines 108,109,110 and
    //comment lines 118,119,120
    Planet* sun= new Planet (1,0,0,0,-0.00136455/365,+0.00206698/365,+2.20476e-005/365);
    Planet* earth= new Planet (3e-6,9.228E-01,3.885E-01,-1.511E-04, -6.888E-03,1.583E-02,-1.138E-06);
    Planet* jupiter= new Planet(9.5e-4,-4.589E+00,-2.915E+00,1.147E-01,3.957E-03,-6.011E-03,-6.358E-05);

    //We add the planets to the class
    solarSystem.add(sun);
    solarSystem.add(earth);
    solarSystem.add(jupiter);


    //solarSystem->add(saturn);
}


std::string three_body_problem::getName() {
    return "Three-body";
}
