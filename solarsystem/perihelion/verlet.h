#ifndef VERLET_H
#define VERLET_H

class Verlet
{
public:
    Verlet();
    // changes the positions and the velocities of all the planets
    // in the solarsystem according to the evolution
    // and to Velocity Verlet algorithm
    void integrate(class SolarSystem *system, double dt);

};

#endif // VERLET_H
