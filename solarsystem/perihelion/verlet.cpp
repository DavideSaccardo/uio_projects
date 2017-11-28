#include "verlet.h"
#include "planet.h"
#include "solarsystem.h"
#include <iostream>
#include <cmath>

using namespace std;

Verlet::Verlet(){

}


void Verlet::integrate(class SolarSystem* system, double dt){

    // memory variable usful to remember the previous velocity
    double memory[system->planets.size()][3];

    // calculate all the forcces acting on each planet
    system->calculateForces();

    for(int i=0; i<system->planets.size(); i++ ){
        for(int j=0; j<3; j++){
            memory[i][j] = (dt*0.5*system->planets[i].force[j])/system->planets[i].mass;
        }
    }

    // new position of each object in the solarsystem
    for(int i=0; i<system->planets.size(); i++){
        for(int j=0; j<3; j++){
            system->planets[i].position[j] +=  dt*system->planets[i].velocity[j] + dt*dt*0.5*system->planets[i].force[j]/system->planets[i].mass;
        }
    }

    // calculate all the forcces in order to update velocity
    system->calculateForces();

    for(int i=0; i<system->planets.size(); i++ ){
        for(int j=0; j<3; j++){
            system->planets[i].velocity[j] += memory[i][j] + dt*0.5*system->planets[i].force[j]/system->planets[i].mass;
        }
    }


}

