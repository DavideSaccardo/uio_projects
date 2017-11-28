#include "pde_diffusion.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath> //libreria migliore
#include <string>
#include <array>
#include "data_analysis.h"

using namespace std;
using std::array;

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


pde_diffusion::pde_diffusion(int n, double a)
{
    N = n;
    alpha = a;
    r_cn = 2-2*alpha;
    r_i = 1 + 2*alpha;
    r_e = 1 -2*alpha;
    r_2d = 1+4*alpha;
}

void pde_diffusion::two_dimension(double** u, double tolerance, int cutoff){
    double v[N+1][N+1], temp[N+1][N+1], delta;
    for(int i=0; i<N+1;i++){
        for(int j=0; j<N+1;j++){
            v[i][j]=u[i][j];
        }
    }

    int iterations = 0;
    double difference=tolerance+1;

    while((iterations<cutoff)&&(difference>tolerance)){
        difference=0;


        for(int i=0; i<N+1;i++){
            for(int j=0; j<N+1;j++){
                temp[i][j]=u[i][j];
            }
        }

        //cout<<"ehi"<<endl;
        for(int i=1;i<N;i++){
            for(int j=1;j<N;j++){
                delta=alpha*(temp[i+1][j]+temp[i-1][j]+temp[i][j-1]+temp[i][j+1]);
                u[i][j]=(delta+v[i][j])/(r_2d);
                difference=fabs(u[i][j]-temp[i][j]);
            }
        }
        iterations++;
        if(iterations>cutoff&&difference>tolerance) cout<<"Failure"<<endl;
        difference/=pow(N-1,2.0);
    }


}



void pde_diffusion::Crank_nicolson(double* u)
{
    double v[N+1];

    for(int i=1; i<N; i++){
        v[i] = alpha*u[i-1] + r_cn*u[i] + alpha*u[i+1];
    }
    v[0] = 0;//10*exp(-1.0/(j+0.01))/sqrt(j+0.01);
    v[N] = 1;//exp(-1.0/(j+0.01))/sqrt(j+0.01);

    tridiagonal_solver(-alpha, 2+2*alpha, -alpha, N, u, v, alpha);
    u[0] = 0;//10*exp(-1.0/(j+0.01))/sqrt(j+0.01);;
    u[N] = 1;//exp(-1.0/(j+0.01))/sqrt(j+0.01);;

}

void pde_diffusion::Implicit(double* u){

    double v[N+1];
    for(int i=1; i<N+1;i++){
        v[i]=u[i];
    }

    tridiagonal_solver(-alpha,r_i,-alpha,N+1,u,v,alpha);
    //cout<<time<<endl;
    u[0]=0;
    u[N]=1;

}

void pde_diffusion::Explicit(double* u){
    double v[N+1];
    for(int i=1; i<N+1;i++){
        v[i]=u[i];
    }
    for(int i=1; i<N; i++){
        u[i] = alpha*v[i-1] + r_e*v[i] + alpha*v[i+1];
    }

}


