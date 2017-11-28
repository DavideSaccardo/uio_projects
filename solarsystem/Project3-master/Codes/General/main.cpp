#include <iostream>
#include "planet.h"
#include "solarsystem.h"
#include "verlet.h"
#include "integrator.h"
#include "initialcondition.h"

using namespace std;

int main(){

    double label = 0.0;
    //Here you can decide which configuration you want to consider, just replace tenBodyProblem by:
    // twoBodyProblem() or twoBodyRange() or threeBodyProblem()

    //When you decide which configuration you want, you have to go inside the equivalent file .cpp and you have to see what it is wrriten inside
    InitialCondition::tenBodyProblem();

    //To choose stepsize go to integral.h
    //To choose number of steps go to initialconditions.cpp, and here you can choose the number of steps for each configuration
    //Remember: number of years=number of steps X step size

    return 0;

}
