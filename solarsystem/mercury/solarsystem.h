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
    private:
        bool                        m_writeToFile       = false;
        bool                        m_outFileOpen       = false;
        long int                         m_integrateSteps    = 0;
        int                         m_number_of_planets = 0;
        double                      m_kineticEnergy     = 0;
        double                      m_angularmomentum   = 0;
        double                      m_totalEnergy       = 0;
        double                      m_potentialEnergy   = 0;
        class Integrator*           m_integrator        = nullptr;
        class Potential*            m_potential         = nullptr;
        class InitialCondition*     m_initialCondition  = nullptr;
        std::ofstream               m_outFile;
        std::vector<Planet*>        m_planets;

    public:
        SolarSystem() {}
        double                      U                   = 0;

        double perihelium =100;
        double time;
        void calculateForces        ();
        void zero_all_Forces        ();
        void setPotential           (class Potential* potential);
        void setIntegrator          (class Integrator* integrator);
        void setInitialCondition    (class InitialCondition* initialCondition);
        void setDt                  (double dt);
        void integrate              (long numberOfSteps);
        void add                    (Planet* p);
        void printIntegrateInfo     (int stepNumber);
        void removeLinearMomentum   ();
        void setFileWriting         (bool writeToFile);
        void writePositionsToFile   ();
        void closeOutFile           ();
        double computeKineticEnergy ();
        void property               ();


        //getter
        int number_of_planets() const;
        void findperihelium();
};




// void write_to_file(double dt);




#endif // SOLARSYSTEM_H
