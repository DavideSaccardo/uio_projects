#include "initialcondition.h"
#include "solarsystem.h"
#include "planet.h"
#include "euler.h"
#include "verlet.h"
#include "newtoniangravity.h"
#include "two_body_problem.h"
#include "three_body_problem.h"
#include "ten_body_problem.h"
#include "time.h"
#include <iostream>
#include <cmath>

std::string InitialCondition::getName() {
    return "Unkown";
}

//Here are defined all the possible configurations implemented in our code
//G is just a constant that we use to accomplish the task about changing the dependence on r of the force from r^2 to r^3, so it is set as zero everywhere except for the configuration prepared to actually do all the calculations about it
//For each configuration you can choose between the solver Verlet or Euler and if you define another class with another potential you can just replace "Newtonian gravity" without affecting all the system
//Moreover you can defide if you want to print or not the positions
//setChangeRange changes the labels of the file where we write to the positions. It's used for the task about beta (dependence on r)
//IN integrate you can choose the number of steps. Note: Number of year=Number of steps* stepsize. The stepsize (dt) is in integrator.h
void InitialCondition::twoBodyProblem() {
    double G = 0.0;

    SolarSystem* twoBodySystem = new SolarSystem();
    twoBodySystem->setIntegrator        (new Verlet(twoBodySystem));
    twoBodySystem->setPotential         (new NewtonianGravity(G));
    twoBodySystem->setInitialCondition  (new two_body_problem());
    twoBodySystem->setFileWriting       (true);
    twoBodySystem->setChangeRange       (false);
    //For this configurations, it computes also the time and it prints it
    clock_t t = clock();
    twoBodySystem->integrate            (1e5);
    double timing = (clock()-t)/CLOCKS_PER_SEC;
    cout <<"Time to run:  " <<setprecision(8)<<timing <<" seconds" <<endl;
}

void InitialCondition::twoBodyRange() {
    //Here we accomplish the task about increasing the exponent of r in the expression of force
    int N = 20;
    double stepsize = 1.0/(double)N;

    for(double kk = 0; kk <= 1.0 + stepsize; kk+=stepsize) {

        double beta = kk;
        double label = N*kk;
        cout<<label<<endl;
        SolarSystem* twoBodySystem = new SolarSystem();
        twoBodySystem->setIntegrator        (new Verlet(twoBodySystem));
        twoBodySystem->setPotential         (new NewtonianGravity(beta));
        twoBodySystem->setInitialCondition  (new two_body_problem());
        twoBodySystem->setFileWriting       (true);
        twoBodySystem->setChangeRange       (true);
        twoBodySystem->setLabel             (label); //it increases the exponent of r for each loop
        twoBodySystem->integrate            (1e5);
    }
}

void InitialCondition::threeBodyProblem() {
    //Here we consider Sun, Earth and Jupiter. Inside three_body_problem you can set the positions to have the Sun fixed or not. By default Sun is not fixed.
    double G = 0.0;

    SolarSystem* threeBodySystem = new SolarSystem();
    threeBodySystem->setIntegrator        (new Verlet(threeBodySystem));
    threeBodySystem->setPotential         (new NewtonianGravity(G));
    threeBodySystem->setInitialCondition  (new three_body_problem());
    threeBodySystem->setFileWriting       (true);
    threeBodySystem->setChangeRange       (false);
    threeBodySystem->integrate            (1e5);
}

void InitialCondition::tenBodyProblem() {
    //Here we consider all the planets (plus Pluto). Inside ten_body_problem you can set the positions to have the Sun fixed or not. By default Sun is not fixed.
    double G = 0.0;

    SolarSystem* tenBodySystem = new SolarSystem();
    tenBodySystem->setIntegrator        (new Verlet(tenBodySystem));
    tenBodySystem->setPotential         (new NewtonianGravity(G));
    tenBodySystem->setInitialCondition  (new ten_body_problem());
    tenBodySystem->setFileWriting       (true);
    tenBodySystem->setChangeRange       (false);
    tenBodySystem->integrate            (300000); //by default these are 300 years, enough for Pluto to finish one revolution

}


