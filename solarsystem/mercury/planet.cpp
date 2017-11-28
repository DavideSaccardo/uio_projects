#include "planet.h"
#include <iostream>
#include <cmath>
using namespace std;

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
    m_mass=m;
    m_x=x;
    m_y=y;
    m_z=z;
    m_vx=vx;
    m_vy=vy;
    m_vz=vz;
}


void Planet::resetForce(){
    m_fx=0;
    m_fy=0;
    m_fz=0;
}

void Planet::addForce(double dFx, double dFy, double dFz) {
    m_fx += dFx;
    m_fy += dFy;
    m_fz += dFz;
}

double Planet::velocitySquared() {
    return m_vx*m_vx +
           m_vy*m_vy +
           m_vz*m_vz;
}

double Planet::angularmomentum(){
    double lx=0;
    double ly=0;
    double lz=0;
    double l=0;

    lx = m_y*m_vz - m_z*m_vy;
    ly = m_z*m_vx - m_x*m_vz;
    lz = m_x*m_vy - m_y*m_vx;

    l = sqrt(lx*lx + ly*ly + lz*lz);

    return l;

}
