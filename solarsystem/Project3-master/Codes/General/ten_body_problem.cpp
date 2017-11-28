#include "ten_body_problem.h"
#include "solarsystem.h"
#include <cmath>

void ten_body_problem::setupParticles(SolarSystem &solarSystem) {

    // all the mass should be given as ratio of the sun's mass
    //          mass   x  y  z  vx vy vz
    // velocity in AU/day

    //Here we declare all the bodies of the solar system.
    //If you want to fix the Sun you have to put its velocities equal to 0, i.e.
    //you have to delete vx, vy, vz of the sun and put 0,0,0, then you have to go to solarsystem->void SolarSystem::calculateForces()
    // and you have to uncomment " //m_planets.at(0)->resetForce(); " this keeps the Sun fixed at the center of our frame. Finally you have
    //to comment lines 113,114,115 and uncomment lines 108,109,110 in planet.cpp

    Planet* sun= new Planet (1,0,0,0,-0.00187358/365,+0.00199705/365, +4.36039e-005/365);
    Planet* earth= new Planet (3e-6,9.228E-01,3.885E-01,-1.511E-04, -6.888E-03,1.583E-02,-1.138E-06);
    Planet* jupiter= new Planet(9.5e-4,-4.589E+00,-2.915E+00,1.147E-01,3.957E-03,-6.011E-03,-6.358E-05);
    Planet* mars= new Planet(3.3e-7, -1.559E+00,5.845E-01,5.031E-02,-4.345E-03,-1.192E-02,-1.433E-04);
    Planet* saturn= new Planet(2.75e-4, -3.580E-01,-1.005E+01,1.890E-01, 5.268E-03,-2.170E-04,-2.057E-04);
    Planet* venus= new Planet(2.45e-6,-6.318E-01,3.418E-01,4.106E-02,-9.565E-03,-1.796E-02,3.053E-04);
    Planet* mercury= new Planet(1.65e-7, -3.413E-01,-2.718E-01,8.702E-03,1.193E-02,-2.062E-02,-2.780E-03);
    Planet* uranus= new Planet(4.4e-5,1.786E+01 ,8.806E+00,-1.987E-01,-1.768E-03,3.344E-03,3.540E-05);
    Planet* neptune= new Planet(5.15e-5,2.861E+01,-8.824E+00,-4.777E-01,9.047E-04,3.019E-03,-8.319E-05);
    Planet* pluto= new Planet(6.55e-9,1.054E+01,-3.171E+01,3.435E-01,3.040E-03,3.236E-04,-9.238E-04);

    //We add the planets to the class
    solarSystem.add(sun);
    solarSystem.add(mercury);
    solarSystem.add(venus);
    solarSystem.add(earth);
    solarSystem.add(mars);
    solarSystem.add(jupiter);
    solarSystem.add(saturn);
    solarSystem.add(uranus);
    solarSystem.add(neptune);
    solarSystem.add(pluto);

}


std::string ten_body_problem::getName() {
    return "Ten-body";
}
