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

//getter and setter
int SolarSystem::number_of_planets() const
{
    return m_number_of_planets;
}


double SolarSystem::getPx() const
{
    return m_px;
}

void SolarSystem::setPx(double value)
{
    m_px = value;
}

double SolarSystem::getPy() const
{
    return m_py;
}

void SolarSystem::setPy(double py)
{
    m_py = py;
}

double SolarSystem::getPz() const
{
    return m_pz;
}

void SolarSystem::setPz(double pz)
{
    m_pz = pz;
}

double SolarSystem::getTotalMass() const
{
    return m_totalMass;
}

void SolarSystem::setTotalMass(double totalMass)
{
    m_totalMass = totalMass;
}

int SolarSystem::getIntegrateSteps() const
{
    return m_integrateSteps;
}

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

int SolarSystem::getLabel() const
{
    return m_label;
}

void SolarSystem::setLabel(double label)
{
    m_label = label;
}

void SolarSystem::setChangeRange(bool changeRange){
    m_changeRange = changeRange;
}

void SolarSystem::add(Planet* planet){
    //function to add planet to the class Planets
    m_planets.push_back(planet); // add an element at the end of the vector
    m_number_of_planets+=1;

    // cout << m_number_of_planets << ": n_planets" << "\t";

}

void SolarSystem::property(){
    //function to get the mass, Vx,Vy,Vz of all the planets. It's useful to see if everything is going well.
    for(int i=0; i < m_number_of_planets; i++){
        cout << m_planets.at(i)->getMass() << "masses" << "\n";
    }
    cout <<"\n" ;
    for(int i=0; i < m_number_of_planets; i++){
        cout << m_planets.at(i)->getVx() <<"--" <<i << "v_x" << "\n";
        cout << m_planets.at(i)->getVy() <<"--"<< i << "v_y" << "\n";
        cout << m_planets.at(i)->getVz() <<"--"<<i << "v_z" << "\n";
    }
    cout <<"\n" ;
}

void SolarSystem::zero_all_Forces(){
    //function that delete all the forces acting on a planet
    for(int i = 0; i < m_number_of_planets; i++){
        m_planets.at(i)->resetForce();
    }
}

void SolarSystem::calculateForces(){
    //function to compute the forces
    zero_all_Forces(); // I delete all the old forces
    m_potential->resetPotentialEnergy();
    for(int i = 0; i < m_number_of_planets; i++){
        for(int j = i+1; j < m_number_of_planets; j++){

            m_potential->calculateForces(m_planets.at(i), m_planets.at(j), U); //compute the force between planet j and planet i
            //m_planets.at(0)->resetForce(); //this function delete the forces acting on the Sun. We use it to fix the Sun in the origin as CM.

            m_potential->setPotentialEnergy(m_potential->getPotentialEnergy()+U);
            //at each step it computes the potential energy. This function is here, because we have the perfect loop to calculate the potential energy without counting terms twice
        }
    }
}

void SolarSystem::integrate(int numberOfSteps) {
    //core of the program where are called all the important calculations
    m_integrateSteps = numberOfSteps;
    printIntegrateInfo(0); //print information on what we are doing at the beginning, just a brief summary
    for (int i=1; i<numberOfSteps+1; i++) {
        m_integrator->integrateOneStep(m_planets); //it uses Euler or Verlet depending on what you have declared at the beginning in "Initial conditions"
        printIntegrateInfo(i); //it prints info about energy, pot. energy, kin. energy, ang.mom at each step
        writePositionsToFile(i); //it opens a file and it prints positions of all the planetson it
    }
    closeOutFile();
}

double SolarSystem::computeKineticEnergy() {
    //function used to compute the kinetic energy of the system
    m_kineticEnergy=0;
    for(int i=0;i<m_number_of_planets;i++){
        m_kineticEnergy += 0.5*m_planets.at(i)->getMass()*m_planets.at(i)->velocitySquared();
    }
    return m_kineticEnergy;
}

double SolarSystem::computeTotalMass(){
    //function used to compute the total mass of the system
    m_totalMass=0;
    for(int i=0;i<m_number_of_planets;i++){
        m_totalMass += m_planets.at(i)->getMass();
    }
    return m_totalMass;
}

void SolarSystem::printIntegrateInfo(int stepNumber) {
    //function used to print useful information such as total energy, potential energy, kinetic energy, angular momentum,
    //coordinates of CM, linear momentum and so on

    if(stepNumber==0||stepNumber==300000){
        cout <<"Center of mass X:    "<<findCenterOfMassX()<<endl;
        cout <<"Center of mass Y:    "<<findCenterOfMassY()<<endl;
        cout <<"Center of mass Z:    "<<findCenterOfMassZ()<<endl;
        computeLinearMomentum();
    }

    if (stepNumber == 0) {
        m_totalEnergy=initialPotentialEnergy()+computeKineticEnergy();
        m_totalL=sqrt(AngularMomentumX()*AngularMomentumX()+AngularMomentumY()*AngularMomentumY()+AngularMomentumZ()*AngularMomentumZ());
        cout << endl
             << " STARTING INTEGRATION "    << endl
             << "-------------------------" << endl
             << "  o Number of steps:     " << m_integrateSteps << endl
             << "  o Time step, dt:       " << m_integrator->getDt() << endl
             << "  o Initial condition:   " << m_initialCondition->getName() << endl
             << "  o Number of planets: " << m_planets.size() << endl
             << "  o Potential in use:    " << m_potential->getName() << endl
             << "  o Integrator in use:   " << m_integrator->getName() << endl
             << "  o Initial angular momentum:   " << setprecision(9) << m_totalL << endl
             << "  o Initial kinetic energy:   " << computeKineticEnergy() << endl
             << "  o Initial potential energy:   " << initialPotentialEnergy() << endl
             << "  o Initial total energy:   " << setprecision(9) << m_totalEnergy << endl
             << endl;


        //        cout<<"Centro di massa x:   "<<m_planets.at(0)->getVx()<<endl;
        //        cout<<"Centro di massa y:   "<<m_planets.at(0)->getVy()<<endl;
        //        cout<<"Centro di massa z:   "<<m_planets.at(0)->getVz()<<endl;

    } else if (stepNumber % 1000 == 0) {
        m_potentialEnergy   = m_potential->getPotentialEnergy();
        m_totalEnergy   = m_kineticEnergy + m_potentialEnergy;
        stepNumber=stepNumber*m_integrator->getDt();
        printf("Year: %5d    E =%10.10f   Ek =%10.5f    Ep =%10.5f\n    L =%10.10f\n",
               stepNumber,m_totalEnergy, m_kineticEnergy, m_potentialEnergy, m_totalL);

        fflush(stdout);
    }

}

double SolarSystem::findCenterOfMassX(){
    //function used to find the X coordinate of the center of Mass
    double cx=0;
    for(int i=0;i<m_number_of_planets;i++){
        cx+=m_planets.at(i)->getX()*m_planets.at(i)->getMass();
    }
    computeTotalMass();
    return cx/(m_totalMass);
}

double SolarSystem::findCenterOfMassY(){
    //function used to find the Y coordinate of the center of Mass
    double cy=0;
    for(int i=0;i<m_number_of_planets;i++){
        cy+=m_planets.at(i)->getY()*m_planets.at(i)->getMass();
    }
    computeTotalMass();
    return cy/(m_totalMass);
}

double SolarSystem::findCenterOfMassZ(){
    //function used to find the Z coordinate of the center of Mass
    double cz=0;
    for(int i=0;i<m_number_of_planets;i++){
        cz+=m_planets.at(i)->getZ()*m_planets.at(i)->getMass();
    }
    computeTotalMass();
    return cz/(m_totalMass);
}

double SolarSystem::AngularMomentumX(){
    //fucntion used to findthe angular momentum along x
    double L_x=0;
    for(int i=0;i<m_number_of_planets;i++){
        L_x+= (m_planets.at(i)->getY()*m_planets.at(i)->getVz() - m_planets.at(i)->getZ()*m_planets.at(i)->getVy())*m_planets.at(i)->getMass();
    }
    return L_x;
}

double SolarSystem::AngularMomentumY(){
    //fucntion used to find the angular momentum along y
    double L_y=0;
    for(int i=0;i<m_number_of_planets;i++){
        L_y+= (m_planets.at(i)->getZ()*m_planets.at(i)->getVx() - m_planets.at(i)->getX()*m_planets.at(i)->getVz())*m_planets.at(i)->getMass();
    }
    return L_y;
}

double SolarSystem::AngularMomentumZ(){
    //fucntion used to find the angular momentum along z
    double L_z=0;
    for(int i=0;i<m_number_of_planets;i++){
        L_z+= (m_planets.at(i)->getX()*m_planets.at(i)->getVy() - m_planets.at(i)->getY()*m_planets.at(i)->getVx())*m_planets.at(i)->getMass();
    }
    return L_z;
}

double SolarSystem::initialPotentialEnergy(){
    //funtion used to find the initial potential energy
    zero_all_Forces(); // delete al forces and reset potential energy: I'm sure I'm calculating the initial
    m_potential->resetPotentialEnergy();
    for(int i = 0; i < m_number_of_planets; i++){
        for(int j = i+1; j < m_number_of_planets; j++){

            double r = (m_planets.at(i)->getX()-m_planets.at(j)->getX())*(m_planets.at(i)->getX()-m_planets.at(j)->getX())+(m_planets.at(i)->getY()-m_planets.at(j)->getY())*(m_planets.at(i)->getY()-m_planets.at(j)->getY())+(m_planets.at(i)->getZ()-m_planets.at(j)->getZ())*(m_planets.at(i)->getZ()-m_planets.at(j)->getZ());

            r = sqrt(r);
            double four_pi_2_massp1_massp2 =(4*M_PI*M_PI)*m_planets.at(i)->getMass()*m_planets.at(j)->getMass();
            U=-four_pi_2_massp1_massp2/r;
            m_potential->setPotentialEnergy(m_potential->getPotentialEnergy()+U);

        }


    }
    return m_potential->getPotentialEnergy();
}


void SolarSystem::computeLinearMomentum() {
    //function used to compute the total linear momentum of the system and we can print the x,y,z coordinates
    double px=0;
    double py=0;
    double pz=0;
    for(int i=0;i<  m_number_of_planets;i++){
        px+=m_planets.at(i)->getMass()*(m_planets.at(i)->getVx());
        py+=m_planets.at(i)->getMass()*(m_planets.at(i)->getVy());
        pz+=m_planets.at(i)->getMass()*(m_planets.at(i)->getVz());

    }

    //   cout<<"linear momentum along x:  "<<px<<endl;
    //   cout<<"linear momentum along y:  "<<py<<endl;
    //   cout<<"linear momentum along z:  "<<pz<<endl;

}

void SolarSystem::setFileWriting(bool writeToFile) {
    m_writeToFile = writeToFile;
}

void SolarSystem::writePositionsToFile(double i) {
    //function used to open a file and to write positions of all the planets in it
    if (m_outFileOpen == false) {
        m_outFile.open("positions.dat", std::ios::out);
        m_outFileOpen = true;
    }

    //    double initial_energy=-0.0044040677080726970544 ;
    //    double initial_L=156.865296259493618435953976586 ;
    //    m_outFile << setw(15) << setprecision(4) << i*m_integrator->getDt() << " \t";
    //    m_outFile << setw(15) << setprecision(30) << log10(fabs((totalAngularMomentum()-initial_L)/(initial_L))) << " \n";
    //above there is just what we used to compute the relative error of energy/ang momentum
    for(int j=0; j < m_planets.size(); j++){

        m_outFile << setw(15) << setprecision(8) << m_planets.at(j)->getX() << " \t";
        m_outFile << setw(15) << setprecision(8) << m_planets.at(j)->getY() << " \t";
        m_outFile << setw(15) << setprecision(8) << m_planets.at(j)->getZ() << " \t";
    }
    m_outFile << "\n";
}

void SolarSystem::closeOutFile() {
    //function to close the file, which was opened with the fucntion above
    if (m_writeToFile == true) {
        m_outFile.close();
        m_outFileOpen = false;
    }
}

