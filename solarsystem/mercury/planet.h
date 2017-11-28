#ifndef PLANET_H
#define PLANET_H


class Planet{

private:
    double m_mass=0;
    double m_x=0;
    double m_y=0;
    double m_z=0;
    double m_vx=0;
    double m_vy=0;
    double m_vz=0;
    double m_fx=0;
    double m_fy=0;
    double m_fz=0;

public:
    /*Planet(double m, double* initialPosition, double* initialVelocity);*/
    Planet(double m, double x, double y, double z, double vx, double vy, double vz);


    void resetForce();
    void addForce(double dFx, double dFy, double dFz);
    double velocitySquared();
    //getters
    double getMass() const;
    double getX() const;
    double getY() const;
    double getZ() const;
    double getVx() const;
    double getVy() const;
    double getVz() const;
    double getFx() const;
    double getFy() const;
    double getFz() const;

    //setters
    void setFx(double fx);
    void setFy(double fy);
    void setFz(double fz);
    void setX(double x);
    void setY(double y);
    void setZ(double z);
    void setVx(double vx);
    void setVy(double vy);
    void setVz(double vz);
    double angularmomentum();
};

#endif // PLANET_H
