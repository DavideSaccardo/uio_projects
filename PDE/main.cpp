#include <iostream>
#include <math.h>
#include "data_analysis.h"
using namespace std;


void tridiagonal_solver(double A, double B, double C, int n, double* S, double* Q, double cost){

    // setting the value of our problem
    double *d = new double[n];
    d[0] = B; // to avoid division by zero in the next loop for
    double *q_new = new double[n + 1];
    q_new[0] = 0;

    // making a new diagonal elements
    for(int i=1; i<n-1; i++){
        d[i] = B - A*C / d[i-1];
        q_new[i] = Q[i] - A*q_new[i-1] / d[i-1];
    }

    d[n-1]=B-(A)*C/d[n-2];
    q_new[n-1] = Q[n-1] - (A)*q_new[n-2] / d[n-2] + cost;
    // substitution for finding S
    S[n-1] = q_new[n-1] / d[n-1];
    for(int i=n-2; i>=0; i--){
        S[i] = ( q_new[i] - C*S[i+1] )/ d[i];
    }

    delete[] d;
    delete[] q_new;

}

int main()
{

    int n=1000;
    double dt, dx;
    int max_index_t =30000;
    dx = 1.0/(n);
    dt = 1e-5;
    double alfa = dt/(dx*dx);
    double* x = new double[n+1];
    for(int i=0; i<=n; i++){
        x[i] = i*dx;
    }
    double* t = new double[max_index_t];
    for(int i=0; i<max_index_t; i++){
        t[i] = i*dt;
    }

    double* v = new double[n+1];
    double* u = new double[n+1];
    //    for(int i=0; i<=n+1; i++){
    //        u[i] = new double[max_index_t];
    //    }

    // initial condition at t=0;
    for(int i=0; i<=n; i++){
        u[i] =0; //exp(-x[i]*x[i]/(0.01))/sqrt(0.01);
    }

    // initial condition at x=0 and x=L
    //for(int i=0; i<max_index_t; i++){
    u[0] = 0;
    u[n] = 1;
    //}

    data_analysis* dat = new data_analysis;
    dat->open_file("pde.txt");

    double r_1 = 2-2*alfa;
    double r_2 = 2+2*alfa;
    for(int j=0; j<max_index_t; j++){
        for(int i=1; i<n; i++){
            v[i] = alfa*u[i-1] + r_1*u[i] + alfa*u[i+1];
        }
        v[0] = 0;//10*exp(-1.0/(j+0.01))/sqrt(j+0.01);
        v[n] = 1;//exp(-1.0/(j+0.01))/sqrt(j+0.01);

        tridiagonal_solver(-alfa, 2+2*alfa, -alfa, n, u, v, alfa);
        u[0] = 0;//10*exp(-1.0/(j+0.01))/sqrt(j+0.01);;
        u[n] = 1;//exp(-1.0/(j+0.01))/sqrt(j+0.01);;
        if(j%100){
            for(int i=0; i<=n; i++){
                dat->write(u[i]);
                //dat->printscreen(u[i], "aa");
                dat->add_column();
            }
            dat->new_row();
        }

    }
    dat->close_file();


    delete dat;
    delete[] x;
    delete[] t;
    //for(int i; i<=n+1; i++) delete[] u[i];
    delete[] u;
    delete[] v;

}

/*
    // test
    double* u = new double[n];
    double* q = new double[n];
    q[0] = 4;
    q[1] = 5;
    // double A, double B, double C, int n, double* S, double* Q, double cost
    tridiagonal_solver(2.0, 1.0, 2.0, 2, u, q, 0.0);
    dat->printscreen(u[0], "a");
    dat->printscreen(u[1], "a");
*/
