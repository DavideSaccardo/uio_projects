#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include "write_to_file.h"

using namespace std;


double angularmomentum(double m_x, double m_y,double m_z,double m_vx, double m_vy, double m_vz){
    double lx=0;
    double ly=0;
    double lz=0;
    double l=0;

    lx = m_y*m_vz - m_z*m_vy;
    ly = m_z*m_vx - m_x*m_vz;
    lz = m_x*m_vy - m_y*m_vx;

    return l = m_mass*sqrt(lx*lx + ly*ly + lz*lz);
}



int main(){

    // number of points
    int n=1000;
    // step size
    double h = 1.0/n;
    // time variable
    double* t = new double[n+1];
    for(int i=0; i<=n; i++){
        t[i] = i*h;
    }
    double x_0 = 0.3075;
    double y_0 = 0;
    double vx_0 = 0;
    double vy_0 = 12.44;

    double* x = new double[n+1];
    double* y = new double[n+1];

    double* vx = new double[n+1];
    double* vy = new double[n+1];

    double* ax = new double[n+1];
    double* ay = new double[n+1];
    // initial condition
    x[0] = x_0;
    y[0] = y_0;
    vx[0] = vx_0;
    vy[0] = vy_0;

    // define the problem
    // F = - k x
    double four_pi_2 = 4*3.1415926*3.1415926;
    double* r_3 = new double[n+1];

    for(int i=0; i<n; i++){
        r_3[i] = (x[i]*x[i] + y[i]*y[i])*sqrt(x[i]*x[i] + y[i]*y[i]);
        ax[i] = - four_pi_2*x[i]/r_3[i];
        ay[i] = - four_pi_2*y[i]/r_3[i];

        x[i+1] = x[i] + t[i+1]*vx[i] + t[i+1]*t[i+1]*0.5*ax[i];
        y[i+1] = y[i] + t[i+1]*vy[i] + t[i+1]*t[i+1]*0.5*ay[i];
        r_3[i+1] = (x[i+1]*x[i+1] + y[i+1]*y[i+1])*sqrt(x[i+1]*x[i+1] + y[i+1]*y[i+1]);

        ax[i+1] = - four_pi_2*x[i+1]/r_3[i+1];
        ay[i+1] = - four_pi_2*y[i+1]/r_3[i+1];

        vx[i+1] = vx[i] + t[i+1]*0.5*(ax[i+1] + ax[i]);
        vy[i+1] = vy[i] + t[i+1]*0.5*(ay[i+1] + ay[i]);
    }


    // G M_sun = four_pi_2

    // earth mass =0.01
    double energy = 0.5*0.01*(vx[0]*vx[0] + vy[0]*vy[0]) - four_pi_2*0.01/sqrt(x[0]*x[0] + y[0]*y[0]); // 1/2 omega^2 m x_0
    cout << energy << endl;
    energy = 0.5*0.01*(vx[99]*vx[99] + vy[99]*vy[99]) - four_pi_2*0.01/sqrt(x[99]*x[99] + y[99]*y[99]);
    cout << energy << endl;


    // matrix of data
    double** data = new double*[n+1];

    for(int i=0; i<n+1; i++){
        data[i] = new double[n+1];
    }

    for(int i=0; i<n+1; i++){
        for(int j=0; j<n+1; j++){
            if(j==0){
                data[i][j] = t[i];
            }
            else if(j==1){
                data[i][j] = x[i];
            }
            else if(j==2){
                data[i][j] = y[i];
            }
            else{
                data[i][j] = 0;
            }

        }
    }

    string pi = "data.txt";
    write_to_file(n+1, data, pi);

    // free space
    for(int i=0; i<=n; i++){
        delete[] data[i];
    }
    delete[] data;



    delete[] r_3;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
    delete[] ax;
    delete[] ay;

    return 0;
}
