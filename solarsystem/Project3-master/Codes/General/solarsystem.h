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
//Here are defined all the variables and all the functions in solarsystem
class SolarSystem
{
private:
    double                      m_label             = 0.0;
    bool                        m_changeRange       = false;
    bool                        m_writeToFile       = false;
    bool                        m_outFileOpen       = false;
    int                         m_integrateSteps    = 0;
    int                         m_number_of_planets = 0;
    double                      m_kineticEnergy     = 0;
    double                      m_totalEnergy       = 0;
    double                      m_potentialEnergy   = 0;
    double                      m_totalMass         = 0;
    double                      m_totalL            = 0;
    class Integrator*           m_integrator        = nullptr;
    class Potential*            m_potential         = nullptr;
    class InitialCondition*     m_initialCondition  = nullptr;
    std::ofstream               m_outFile;
    std::vector<Planet*>        m_planets;
    double                      m_px                = 0;
    double                      m_py                = 0;
    double                      m_pz                = 0;

public:
    SolarSystem() {}
    double                      U                   = 0;
    void setChangeRange         (bool changeRange);
    void calculateForces        ();
    void zero_all_Forces        ();
    void setPotential           (class Potential* potential);
    void setIntegrator          (class Integrator* integrator);
    void setInitialCondition    (class InitialCondition* initialCondition);
    void setDt                  (double dt);
    void integrate              (int numberOfSteps);
    void add                    (Planet* p);
    void printIntegrateInfo     (int stepNumber);
    void computeLinearMomentum  ();
    void setFileWriting         (bool writeToFile);
    void writePositionsToFile   (double i);
    void closeOutFile           ();
    double computeKineticEnergy ();
    double computeTotalMass     ();
    void property               ();
    double findCenterOfMassX    ();
    double findCenterOfMassY    ();
    double findCenterOfMassZ    ();
    double totalAngularMomentum ();
    double initialPotentialEnergy();
    double AngularMomentumX     ();
    double AngularMomentumY     ();
    double AngularMomentumZ     ();

    //getters and setters
    int number_of_planets() const;
    double getPx() const;
    void setPx(double value);
    double getPy() const;
    void setPy(double py);
    double getPz() const;
    void setPz(double pz);
    double getTotalMass() const;
    void setTotalMass(double totalMass);
    int getIntegrateSteps() const;
    int getLabel() const;
    void setLabel(double label);
};

#endif // SOLARSYSTEM_H
