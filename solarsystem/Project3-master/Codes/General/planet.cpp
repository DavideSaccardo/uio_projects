#include "planet.h"
#include <iostream>
#include <cmath>

using namespace std;

//getter and setter
double Planet::getX() const
{
    return m_x;
}

double Planet::getY() const
{
    return m_y;
}

double Planet::getZ() const
{
    return m_z;
}

double Planet::getVx() const
{
    return m_vx;
}

double Planet::getVy() const
{
    return m_vy;
}

double Planet::getVz() const
{
    return m_vz;
}

double Planet::getFx() const
{
    return m_fx;
}

double Planet::getFy() const
{
    return m_fy;
}

double Planet::getFz() const
{
    return m_fz;
}


double Planet::getMass() const
{
    return m_mass;
}

void Planet::setFx(double fx)
{
    m_fx = fx;
}

void Planet::setFy(double fy)
{
    m_fy = fy;
}

void Planet::setFz(double fz)
{
    m_fz = fz;
}

void Planet::setX(double x)
{
    m_x = x;
}

void Planet::setY(double y)
{
    m_y = y;
}

void Planet::setZ(double z)
{
    m_z = z;
}

void Planet::setVx(double vx)
{
    m_vx = vx;
}

void Planet::setVy(double vy)
{
    m_vy = vy;
}

void Planet::setVz(double vz)
{
    m_vz = vz;
}

Planet::Planet(double m, double x, double y, double z, double vx, double vy, double vz){
    //function to give to the planet we are considering the mass, positions, velocities in the class Planet.h
    m_mass=m;
    //use this when you want the sun fixed at the origin of our system
//    m_x=x;
//    m_y=y;
//    m_z=z;

    //TEN BODIES: use this with all the bodies when you want the CM fixed and Sun free to move
        m_x=+0.00219511;
        m_y=+0.00559061;
        m_z=-0.000127547;

    //THREE BODIES:: use this with three bodies when you want the CM fixed and Sun free to move
    //    m_x=+0.00435263;
    //    m_y=+0.00276545;
    //    m_z=-0.000108861

    m_vx=vx*365;
    m_vy=vy*365;
    m_vz=vz*365;
}

void Planet::resetForce(){
    //function to reset forces acting on a planet
    m_fx=0;
    m_fy=0;
    m_fz=0;
}

void Planet::addForce(double dFx, double dFy, double dFz) {
    //function to add forces on a planet
    m_fx += dFx;
    m_fy += dFy;
    m_fz += dFz;
}

double Planet::velocitySquared() {
    return m_vx*m_vx +
            m_vy*m_vy +
            m_vz*m_vz;
}


