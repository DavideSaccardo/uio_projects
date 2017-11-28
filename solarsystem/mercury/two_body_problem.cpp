#include "two_body_problem.h"
#include "solarsystem.h"
#include <cmath>

void two_body_problem::setupParticles(SolarSystem &solarSystem) {

    // all the mass should be given as ratio of the sun's mass
    //          mass   x  y  z  vx vy vz
    // velocity in AU/day

    Planet* mercury = new Planet(1.65e-7,0.3075,0,0,0,12.44,0);
    Planet* sun= new Planet (1, 0,0,0,0,0,0 );//1.142200785006636E-03*365.15, -4.208281605957696E-05*365.15, -1.421958209178023E-04*365.15);
    //Planet jupiter(9.5e-4, -4.595050077098597E+00, -2.919444990679906E+00,  1.147977944379108E-01,  5.096929224228880E-03, -6.063762480109976E-03, -2.056733108221762E-04);
    //Planet saturn(2.75e-4,+6.164409988928,+6.366764207406,+2.364527399872,-4.426842485459E-03,+3.394095269623E-03,+1.592276937224E-03);


    solarSystem.add(sun); // has to be the first element
    solarSystem.add(mercury);
    //solarSystem.add(jupiter);
    //solarSystem->add(saturn);
}
//    Verlet* verlet = new Verlet();

//    int n=1000;
//    double h=1.0/n;
//    double t=0; // years


//    for(int i=1; i<=n; i++){
//        t=i*h;

//        verlet->integrate(solarSystem, t);
//        solarSystem->write_to_file(t);

//    }
//    cout << t << endl;
//    //for i in timesteps:
//    // verlet->integrate(solarSystem);
//    //

//    delete solarSystem;
//    delete verlet;
//    //delete euler;
//}

//}

std::string two_body_problem::getName() {
    return "Two-body";
}
