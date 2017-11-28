#include "planet.h"
#include "planet.h"

// calss that contains all the properies of each planet
Planet::Planet(double m, double x, double y, double z, double vx, double vy, double vz){
    // mass in terms of the mass of the Sun
    mass = m;
    // position in AU
    position[0] = x;
    position[1] = y;
    position[2] = z;
    // velocity in AU/years
    velocity[0] = vx;
    velocity[1] = vy;
    velocity[2] = vz;

}

// calculate the magnitude of the total angular momentum over the mass of the planet
double Planet::angularmomentum(){
    double lx=0;
    double ly=0;
    double lz=0;
    double l=0;

    lx = position[1]*velocity[2] - position[2]*velocity[1] ;
    ly = position[2]*velocity[0] - position[0]*velocity[2];
    lz = position[0]*velocity[1] - position[1]*velocity[0];

    l = lx*lx + ly*ly + lz*lz;

    return l;

}

// set to zero the forces acting on the planet
void Planet::resetForce(){
    for(int i=0; i<3; i++)
        force[i] = 0;
}
