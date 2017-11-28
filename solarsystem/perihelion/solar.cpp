#include <iostream>
#include "planet.h"
#include "solarsystem.h"
#include "verlet.h"
#include <cmath>
using namespace std;

int main(){
    // study of the precession of Mercury's perihelion
    // we studied before the interval of the first perihelion
    // in order to easily find it in the interval (100.14, 100.16)
    // all the mass should be given as ratio of the sun's mass
    // velocity in AU/day vx vy vz

    //              (mass,   x , y , z , vx ,vy ,vz)
    Planet mercury(1.65e-7,0.3075,0,0,0,12.44,0);
    Planet sun(1, 0,0,0,0,0,0 );

    SolarSystem* solarSystem = new SolarSystem();
    solarSystem->add(sun); // has to be the first element
    solarSystem->add(mercury);

    Verlet* verlet = new Verlet();

    double h=1.0/2e7; // step-size
    double t=0; // initialization of the time evolution
    double years = 100.16; // how long you want the system evolve in years
    double start = 100.14; // when you want to start looking for the perihelion

    // initialization of the angle of the first perihelion after 1 century
    double angle = 0;
    // if switch_correction(1) you use the relativistic correction
    // if switch_correction(0) yuo use Newton's force
    solarSystem->switch_correction(1);

    double mem=1e5; // variable usful to find the the closest point to the sun

    // main cycle where you look for the first perihelion after one century
    while( t<years){

        t+=h;

        verlet->integrate(solarSystem, h);

        if(t>start ){
            //solarSystem->write_to_file(t, "einstein.txt");
            if( (solarSystem->r_calculator(solarSystem->planets[0], solarSystem->planets[1]) )<mem){
                mem = solarSystem->r_calculator(solarSystem->planets[0], solarSystem->planets[1]);
            }
        }
    }

    angle =  solarSystem->planets[1].position[1]/solarSystem->planets[1].position[0] ;

    // angle given in arcseconds
    cout << atan(angle)*180*3600/M_PI << "''" << endl;


    delete solarSystem;
    delete verlet;
    //delete euler;
}
