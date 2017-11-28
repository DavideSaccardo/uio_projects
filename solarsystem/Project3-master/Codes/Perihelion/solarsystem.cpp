#include "planet.h"
#include "solarsystem.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <ios>

using namespace std;
ofstream ofile;
SolarSystem::SolarSystem(){

}

// add planets to the vector planets
void SolarSystem::add(Planet& planet){
    planets.push_back(planet); // add an element at the end of the vector
}

// the following function sets alle the forces to zero
void SolarSystem::zeroForces(){
    for(int i = 0; i < planets.size(); i++){ // to compile : c++ -std=c++11 -o t.x *.cpp
        planets[i].resetForce();
    }
}

void SolarSystem::calculateForces(){
    // initializing to zero all the forces
    zeroForces();

    for(int i = 0; i < planets.size(); i++){
        for(int j = i+1; j < planets.size(); j++){
            /*
    for(Planet planets[0] : planets){ // to compile : c++ -std=c++11 -o t.x *.cpp
        for(Planet planets[0] : planets){ // it creates only a new planets
            */
            // no forces acting on itself

            calculatePairForce(planets[i], planets[j]);


            //cout << planets[j].mass << "mass loop" << endl;
            // calculate the force due to p2 acting on p1


        }
    }

}

// function to switch between relativistic correction and newtonian force
void SolarSystem::switch_correction(int n){
    if(n==1){
        target =1;
    }
    else{
        target =0;
    }
}

// calculator for the distance between two planets
double SolarSystem::r_calculator(Planet p1, Planet p2){
    // calcoli distanza come vettore
    double r;
    r=0;
    for(int i=0; i<3; i++){
        r += (p1.position[i]-p2.position[i])*(p1.position[i]-p2.position[i]);
    }
    r = sqrt(r);
    return r;
}

// function that calculates the forces
void SolarSystem::calculatePairForce(Planet &p1, Planet &p2){

    double r;
    double four_pi_2_massp1_massp2 = 4*M_PI*M_PI*p1.mass*p2.mass;
    double c = 63197.8;

    if(target == 1){
        gamma = (1 + 3*p2.angularmomentum()/(c*c*r*r));
    }
    else{
        gamma = 1;

    }

    // F = -four_pi_2 / r^3 (x_p1 - x_p2) M_p1 M_p2 /M_sun
    // every mass is in term of the mass of the sun so
    // M_p1 M_p2 /M_sun =M_p1 M_p2

    r = r_calculator(p1,p2);

    p1.force[0] = -four_pi_2_massp1_massp2*(p1.position[0]-p2.position[0])/(r*r*r);
    p1.force[1] = -four_pi_2_massp1_massp2*(p1.position[1]-p2.position[1])/(r*r*r);
    p1.force[2] = -four_pi_2_massp1_massp2*(p1.position[2]-p2.position[2])/(r*r*r);
    for(int i=0; i<3; i++){
        p2.force[i] = -p1.force[i]*gamma;
    }


}







// ADDITIONAL FUNCTIONS

// write to a file data

void SolarSystem::write_to_file(double dt, std::string aa){
    std::ofstream log(aa, std::ios_base::app | std::ios_base::out);

    log << setw(15) << setprecision(8) << dt << ", \t";

    for(int j=0; j < planets.size(); j++){
        for(int k=0; k < 3; k++){
            log << setw(15)  << planets[j].position[k] << ", \t";
        }
    }
    //log << setw(15) << setprecision(8) << r_calculator(planets[0], planets[1]) << ", \t";
    log << "\n";

}


// useful function if we want to control the conservation of energy
// or position and velocity

double SolarSystem::Energy(){
    double kinetic = 0;
    double potential = 0;

    for(int i = 0; i < planets.size(); i++){
        for(int j = 0; j < 3; j++){
            kinetic += planets[i].velocity[j]*planets[i].velocity[j]*planets[i].mass*0.5;
        }
    }

    for(int i = 0; i < planets.size(); i++){
        for(int j = i+1; j < planets.size(); j++){
            double four_pi_2_massp1_massp2 = 4*M_PI*M_PI*planets[i].mass*planets[j].mass;
            potential += -four_pi_2_massp1_massp2/r_calculator(planets[i], planets[j]);
        }
    }

    double E=kinetic+potential;
    return E;

}


void SolarSystem::property(){
    for(int i=0; i < planets.size(); i++){
        cout << planets[i].mass << "masses" << "\t";
    }
    cout <<"\n" ;
    for(int i=0; i < planets.size(); i++){
        cout << planets[i].velocity[1] << "velocity" << "\n";
        cout << planets[i].position[0] << "position" << "\n";
    }
    cout <<"\n" ;

    cout << "energy \t" << Energy() << endl;

}
