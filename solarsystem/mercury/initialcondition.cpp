#include "initialcondition.h"
#include "solarsystem.h"
#include "planet.h"
#include "euler.h"
#include "verlet.h"
#include "newtoniangravity.h"
#include "two_body_problem.h"

#include <iostream>
#include <cmath>

std::string InitialCondition::getName() {
    return "Unkown";
}
// f(, verlet, false, 100)
void InitialCondition::twoBodyProblem() {
    double G = 1.0;

    SolarSystem* twoBodySystem = new SolarSystem();
    twoBodySystem->setIntegrator        (new Verlet(twoBodySystem));
    twoBodySystem->setPotential         (new NewtonianGravity(G));
    twoBodySystem->setInitialCondition  (new two_body_problem());
    twoBodySystem->setFileWriting       (true);
    //twoBodySystem->removeLinearMomentum ();
    twoBodySystem->integrate            (1.0016E+10);
}
