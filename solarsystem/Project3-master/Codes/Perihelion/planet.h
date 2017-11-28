#ifndef PLANET_H
#define PLANET_H


class Planet{
public:
    /*Planet(double m, double* initialPosition, double* initialVelocity);*/
    Planet(double m, double x, double y, double z, double vx, double vy, double vz);
    double mass;
    double position[3];
    double velocity[3];
    double force[3];
    void resetForce();
    double angularmomentum();
private:
};

#endif // PLANET_H
