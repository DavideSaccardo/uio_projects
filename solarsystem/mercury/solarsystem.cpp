#include "planet.h"
#include "solarsystem.h"
#include "initialcondition.h"
#include "integrator.h"
#include "newtoniangravity.h"
#include "potential.h"
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

void SolarSystem::add(Planet* planet){
    m_planets.push_back(planet); // add an element at the end of the vector
    m_number_of_planets+=1;
    // cout << m_number_of_planets << ": n_planets" << "\t";

}

void SolarSystem::property(){
    for(int i=0; i < m_number_of_planets; i++){
        cout << m_planets.at(i)->getMass() << "masses" << "\t";
    }
    cout <<"\n" ;
    for(int i=0; i < m_number_of_planets; i++){
        cout << m_planets.at(i)->getVx() << "v_x" << "\t";
        cout << m_planets.at(i)->getVy() << "v_y" << "\t";
        cout << m_planets.at(i)->getVz() << "v_z" << "\t";
    }
    cout <<"\n" ;
}

int SolarSystem::number_of_planets() const
{
    return m_number_of_planets;
}

void SolarSystem::zero_all_Forces(){
    for(int i = 0; i < m_number_of_planets; i++){ // to compile : c++ -std=c++11 -o t.x *.cpp
        m_planets.at(i)->resetForce();
    }
}


void SolarSystem::calculateForces(){
zero_all_Forces(); // voglio azzerare per non sommare alle vecchie forze
    m_potential->resetPotentialEnergy();
    for(int i = 0; i < m_number_of_planets; i++){

        for(int j = i+1; j < m_number_of_planets; j++){

            m_potential->calculateForces(m_planets.at(i), m_planets.at(j), U);

            //cout<<m_planets.at(2)->getFx()<<"Fx"<<i<<j<<endl;
//            if(i==0&&j==1){
//                cout<<"sole-terra"<<endl;
//            } else if(i==0&&j==2){
//                cout<<"sole-giove"<<endl;
//            } else {
//                cout<<"terra-giove"<<endl;
//            }
        m_potential->setPotentialEnergy(m_potential->getPotentialEnergy()+U);

            //void calculatePairForce(int i, int j);
            //cout << planets[j].mass << "mass loop" << endl;
            // calculate the force due to p2 acting on p1
            //}

        }
        //m_planets.at(i)->addForce(m_planets.at(i)->getFx(), m_planets.at(i)->getFy(),m_planets.at(i)->getFz());
        //cout<<m_planets.at(i)->getFx()<<"Fx"<<i<<endl;

    }
   // cout<<m_planets.at(1)->getFx()<<"Fx"<<endl;
}
/*
void SolarSystem::calculatePairForce(int i, int j){

}
*/

void SolarSystem::setPotential(Potential* potential) {
    m_potential = potential;
}

void SolarSystem::setIntegrator(Integrator* integrator) {
    m_integrator = integrator;
}

void SolarSystem::setInitialCondition(InitialCondition* initialCondition) {
    m_initialCondition = initialCondition;
    m_initialCondition->setupParticles(*this);
}

void SolarSystem::setDt(double dt) {
    m_integrator->setDt(dt);
}

void SolarSystem::integrate(long int numberOfSteps) {
    m_integrateSteps = numberOfSteps;
    time =0;
    printIntegrateInfo(0);
    for (int i=1; i<numberOfSteps+1; i++) {
        time += 0.0001;
        m_integrator->integrateOneStep(m_planets);
        printIntegrateInfo(i);
        writePositionsToFile();
    }
    closeOutFile();
}

double SolarSystem::computeKineticEnergy() {
    /*
     * Here, the kinetic energy of the entire system should be computed. Since
     * this is independent of the potential in use, we place this method
     * directly in the system class.
     *
     * Remember that you can access the mass and velocity of particle i by
     *
     *      m_particles.at(i)->getMass()
     *      m_particles.at(i)->getVelocity()
     *
     * Remember also that the Particle class has a built in method
     * Particle::velocitySquared which can be used here.
     */
    m_kineticEnergy=0;
    for(int i=0;i<m_number_of_planets;i++){
        m_kineticEnergy += 0.5*m_planets.at(i)->getMass()*m_planets.at(i)->velocitySquared();
    }
    return m_kineticEnergy;
}


void SolarSystem::printIntegrateInfo(int stepNumber) {
    if (stepNumber == 0) {
        cout << endl
             << " STARTING INTEGRATION "    << endl
             << "-------------------------" << endl
             << "  o Number of steps:     " << m_integrateSteps << endl
             << "  o Time step, dt:       " << m_integrator->getDt() << endl
             << "  o Initial condition:   " << m_initialCondition->getName() << endl
             << "  o Number of planets: " << m_planets.size() << endl
             << "  o Potential in use:    " << m_potential->getName() << endl
             << "  o Integrator in use:   " << m_integrator->getName() << endl
             << endl;
    }else if (stepNumber % 1000 ==0){
        /*
        m_potentialEnergy   = m_potential->getPotentialEnergy();
        m_totalEnergy       = m_kineticEnergy + m_potentialEnergy;
        printf("Step: %5d    E =%10.5f   Ek =%10.5f    Ep =%10.5f   \n",
               stepNumber,m_totalEnergy, m_kineticEnergy, m_potentialEnergy);
        fflush(stdout);*/

    }


    if( time >100){
        findperihelium();
    }


}


void SolarSystem::findperihelium(){

    double position = sqrt((m_planets.at(1)->getX()*m_planets.at(1)->getX()) + (m_planets.at(1)->getY()*m_planets.at(1)->getY())+(m_planets.at(1)->getZ()*m_planets.at(1)->getZ()));

//    if(position < sperihelium){
//        perihelium = position;
//        printf("ratio: %10.5f \n" , atan2(m_planets.at(1)->getY(),m_planets.at(1)->getX()) * 180 /M_PI);
//    }

}


//void SolarSystem::write_to_file( double dt){
//    std::ofstream log("data.txt", std::ios_base::app | std::ios_base::out);

//    log << setw(15) << setprecision(8) << dt << ", \t";

//    for(int j=0; j < m_number_of_planets; j++){
//            log << setw(15) << setprecision(8) << m_planets.at(j)->getX() << ", \t";
//             log << setw(15) << setprecision(8) << m_planets.at(j)->getY() << ", \t";
//              log << setw(15) << setprecision(8) << m_planets.at(j)->getZ() << ", \t";
//    }
//    log << "\n";

//}

void SolarSystem::removeLinearMomentum() {
    /*
     * Here you should remove the total momentum of the entire system, to
     * ensure the entire system does not drift away during long integration
     * times.
     *
     * Remember that you can access the mass and velocity of particle i by
     *
     *      m_particles.at(i)->getMass();
     *      m_particles.at(i)->getVelocity();
     *
     * Remember also that the vec3-vector class supports the += and -=
     * operators, so you can do
     *
     *      totalMomentum += p->getVelocity() * p->getMass();
     */
    double px=0;
    double py=0;
    double pz=0;
    for(int i=0;i<  m_number_of_planets;i++){
        px+=m_planets.at(i)->getMass()*(m_planets.at(i)->getVx());
        py+=m_planets.at(i)->getMass()*(m_planets.at(i)->getVy());
        pz+=m_planets.at(i)->getMass()*(m_planets.at(i)->getVz());
    }

    for(int i=0; i<m_number_of_planets;i++){
       m_planets.at(i)->setVx(m_planets.at(i)->getVx()-px/(1+1.65e-7));
       m_planets.at(i)->setVy(m_planets.at(i)->getVy()-py/(1+1.65e-7));
       m_planets.at(i)->setVz(m_planets.at(i)->getVz()-pz/(1+1.65e-7));
    }
    //    for(int i=0;i<  m_number_of_planets;i++){
    //        m_planets.at(i)->setVx(m_planets.at(i)->getVx()-)*(m_planets.at(i)->getVx());
    //        py+=m_planets.at(i)->getMass()*(m_planets.at(i)->getVy());
    //        pz+=m_planets.at(i)->getMass()*(m_planets.at(i)->getVz());
    //    }
    //        for(int i=0;i<  m_number_of_planets;i++){
    //            m_planets.at(i)->setVx(0);
    //            m_planets.at(i)->setVy(0);
    //            m_planets.at(i)->setVz(0);
    //     }


}

void SolarSystem::setFileWriting(bool writeToFile) {
    m_writeToFile = writeToFile;
}

void SolarSystem::writePositionsToFile() {
    if (m_outFileOpen == false) {
        m_outFile.open("positions_E.dat", std::ios::out);
        m_outFileOpen = true;
    }
    /*
     * This is where you should print the positions of each particle to file.
     * Note that the file, "positions.dat", is already open; it is opened in
     * the above if-test the first time this method is called in
     * System::Integrate.
     *
     * Which format you choose for the data file is up to you.
     */
    // std::ofstream log("test.txt", std::ios_base::app | std::ios_base::out);

    // m_outFile << setw(15) << setprecision(8) << m_integrator->getDt() << ", \t";

    for(int j=0; j < m_planets.size(); j++){

        m_outFile << setw(15) << setprecision(8) << m_planets.at(j)->getX() << " \t";
        m_outFile << setw(15) << setprecision(8) << m_planets.at(j)->getY() << " \t";
        m_outFile << setw(15) << setprecision(8) << m_planets.at(j)->getZ() << " \t";
    }
    m_outFile << "\n";
}

void SolarSystem::closeOutFile() {
    if (m_writeToFile == true) {
        m_outFile.close();
        m_outFileOpen = false;
    }
}

