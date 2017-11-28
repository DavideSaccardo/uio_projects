#include "newtoniangravity.h"
#include <iostream>
#include <cmath>

using namespace std;

NewtonianGravity::NewtonianGravity(double G) : m_G(G) {

}

void NewtonianGravity::calculateForces(Planet* p1, Planet* p2, double& U) {
    /*
     * This is where the ordinary Newtoninan gravity forces and potential
     * energies should be calculated. This method is called by the System
     * class in System::computeForces, for all particle pairs a and b.
     *
     * Note that you may access the mass and the position of the particles a
     * and b by
     *
     *      a.getMass();       b.getMass();
     *      a.getPosition();   b.getPosition();
     *
     * In order to apply the forces to each particle, it is easiest to use the
     * Particle::addForce method.
     *
     * Since calculating the forces between a and b almost inevitably involves
     * calculating the potential energy, V(r), it is assumed by the Potential
     * class that this quantity is calculated here and added to the
     * m_potentialEnergy variable. Note: You may skip this until you have a
     * working two-body problem, since the calculation of the potential energy
     * is only neccessary for verification purposes later.
     */

    // ...
    //m_potentialEnergy += V;
    //a.addForce(dFx, dFy, dFz);
    //b.addForce(...);



    // calcoli distanza come vettore
    double r;

    r = (p1->getX()-p2->getX())*(p1->getX()-p2->getX())+(p1->getY()-p2->getY())*(p1->getY()-p2->getY())+(p1->getZ()-p2->getZ())*(p1->getZ()-p2->getZ());

    r = sqrt(r);

    // calclolo modulo distanza
    // calculate forza
    // F = -four_pi_2 / r^3 (x_p1 - x_p2) M_p1 M_p2 /M_sun
    // every mass is in term of the mass of the sun so
    // M_p1 M_p2 /M_sun =M_p1 M_p2
    double four_pi_2_massp1_massp2 =(4*M_PI*M_PI)*p1->getMass()*p2->getMass();

    //    cout << p1->getMass() <<"mass pair1"<< endl;
    //    cout << p2->getMass() <<"mass pair2"<< endl;

    //    cout << four_pi_2_massp1_massp2 <<"cost" << endl;

    double c = 63197.8; // speed of light in AU/y

    double gamma = (1 + 3*p2->angularmomentum()*p2->angularmomentum()/(r*r*c*c));

    //cout << 3*p2->angularmomentum()*p2->angularmomentum()/(r*r*c*c) << endl;

    p1->addForce(-four_pi_2_massp1_massp2*(p1->getX()-p2->getX())/(r*r*r),-four_pi_2_massp1_massp2*(p1->getY()-p2->getY())/(r*r*r),-four_pi_2_massp1_massp2*(p1->getZ()-p2->getZ())/(r*r*r));

    p2->addForce(+four_pi_2_massp1_massp2*(p1->getX()-p2->getX())*gamma/(r*r*r),+four_pi_2_massp1_massp2*(p1->getY()-p2->getY())*gamma/(r*r*r),+four_pi_2_massp1_massp2*(p1->getZ()-p2->getZ())*gamma/(r*r*r));

   U=-four_pi_2_massp1_massp2/r;

}

std::string NewtonianGravity::getName() {
    return "Newtonian gravity";
}
