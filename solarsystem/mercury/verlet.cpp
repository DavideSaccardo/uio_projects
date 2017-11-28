#include "verlet.h"
#include "planet.h"
#include "solarsystem.h"
#include <iostream>
#include <cmath>

using namespace std;

Verlet::Verlet(SolarSystem* solarSystem)
    : Integrator(solarSystem) {
}

std::string Verlet::getName() {
    return "Velocity verlet";
}


void Verlet::integrateOneStep(std::vector<Planet *> planets){



    // mi serve
    // di chi è la prop.
    // devo salvarlo
    // mi potrebbe riservire?
    // posso dividere in sottoblocchi?


    // position

    // why a declaration of a pointer here changes variables
    // if i declare memory as a matrix pointer it will change the values in solarsystem
    // to try uncomment the last part and create the pointer
    // constant

    double ** memory = new double*[planets.size()];
    for(int i = 0; i <planets.size(); i++){
        memory[i] = new double[3];
    }

    //cout<<planets.size()<< endl;


    for(int i=0; i<planets.size(); i++ ){
            memory[i][0] = (m_dt*0.5*planets.at(i)->getFx())/planets.at(i)->getMass();
            memory[i][1] = (m_dt*0.5*planets.at(i)->getFy())/planets.at(i)->getMass();
            memory[i][2] = (m_dt*0.5*planets.at(i)->getFz())/planets.at(i)->getMass();

            //matrix which keeps in mind acceleration a_i
    }

    m_solarSystem->calculateForces(); //

    for(int i=0; i<planets.size(); i++){
        planets.at(i)->setX(planets.at(i)->getX()+m_dt*planets.at(i)->getVx() + m_dt*m_dt*0.5*planets.at(i)->getFx()/planets.at(i)->getMass());
        planets.at(i)->setY(planets.at(i)->getY()+m_dt*planets.at(i)->getVy() + m_dt*m_dt*0.5*planets.at(i)->getFy()/planets.at(i)->getMass());
        planets.at(i)->setZ(planets.at(i)->getZ()+m_dt*planets.at(i)->getVz() + m_dt*m_dt*0.5*planets.at(i)->getFz()/planets.at(i)->getMass());
    }

    m_solarSystem->calculateForces();

    for(int i=0; i<planets.size(); i++ ){
        planets.at(i)->setVx(planets.at(i)->getVx()+memory[i][0] + m_dt*0.5*planets.at(i)->getFx()/(planets.at(i)->getMass()));
        planets.at(i)->setVy(planets.at(i)->getVy()+memory[i][1] + m_dt*0.5*planets.at(i)->getFy()/(planets.at(i)->getMass()));
        planets.at(i)->setVz(planets.at(i)->getVz()+memory[i][2] + m_dt*0.5*planets.at(i)->getFz()/(planets.at(i)->getMass()));
            // v_(i+1)=memory+(h/2)*a_(i+1), memory=(h/2)*a_i-

    }

    //cout << "end" << endl;
    //system->property();

    m_solarSystem->computeKineticEnergy();
    // free space with problems

    //m_solarSystem->removeLinearMomentum();

    for (int i = 0; i <planets.size(); i++){
        delete[] memory[i];
    }

    delete[] memory;


    //cout<<"funziono"<<endl;


}

