#include "newtoniangravity.h"
#include <iostream>
#include <cmath>

using namespace std;

NewtonianGravity::NewtonianGravity(double G) : m_G(G) {

}

void NewtonianGravity::calculateForces(Planet* p1, Planet* p2, double& U) {
    //Function used to compute the force between planet p1 and planet p2

    double r; //Here it computes the distances
    r = (p1->getX()-p2->getX())*(p1->getX()-p2->getX())+(p1->getY()-p2->getY())*(p1->getY()-p2->getY())+(p1->getZ()-p2->getZ())*(p1->getZ()-p2->getZ());
    r = sqrt(r);

    // calculate force
    // F = -four_pi_2 / r^3 (x_p1 - x_p2) M_p1 M_p2 /M_sun
    // every mass is in term of the mass of the sun so
    // M_p1 M_p2 /M_sun =M_p1 M_p2

    double rfactor = 1/pow(r,m_G); //useful parameter for the task about increasing the exponential of r in the force from 2 to three, in all the other cases it's just rfactor=1
    double four_pi_2_massp1_massp2 =rfactor*(4*M_PI*M_PI)*p1->getMass()*p2->getMass();

    //it uses add forces to add to a Planet all the forces acting on it
    p1->addForce(-four_pi_2_massp1_massp2*(p1->getX()-p2->getX())/(r*r*r),-four_pi_2_massp1_massp2*(p1->getY()-p2->getY())/(r*r*r),-four_pi_2_massp1_massp2*(p1->getZ()-p2->getZ())/(r*r*r));
    p2->addForce(+four_pi_2_massp1_massp2*(p1->getX()-p2->getX())/(r*r*r),+four_pi_2_massp1_massp2*(p1->getY()-p2->getY())/(r*r*r),+four_pi_2_massp1_massp2*(p1->getZ()-p2->getZ())/(r*r*r));

    U=-four_pi_2_massp1_massp2/r*rfactor; //It computes the potential energy
}

std::string NewtonianGravity::getName() {
    return "Newtonian gravity";
}
