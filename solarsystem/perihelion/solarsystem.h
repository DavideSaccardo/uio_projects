#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H
#include <vector>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <ios>
#include "verlet.h"
#include "planet.h"

using namespace std;

class SolarSystem
{
public:
    SolarSystem();
    void add(Planet &planet); // add a planet

    std::vector<class Planet> planets;
    void calculateForces();

    double target;
    double gamma;

    void write_to_file(double dt, std::string aa);
    void switch_correction(int n);
    void Perihelium(Planet &p1);
    void property();
    double Energy();
    double r_calculator(Planet p1, Planet p2);
private:
    void calculatePairForce(Planet &p1, Planet &p2);
    //void calculatePairForce(int i, int j);
    void zeroForces();
};

#endif // SOLARSYSTEM_H
