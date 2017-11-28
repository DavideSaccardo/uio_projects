#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>

void tridiagonal_solver(double, double, double, int, double*, double*, double);

using namespace std;
ofstream ofile;

int main()
{
    //set variables
    int n=1e2;
    int time_resol=3e3;
    double initial_x=0;
    double final_x=1;
    double step_x=(final_x-initial_x)/(n);
    double initial_time=0;
    double final_time=3;
    double step_time=(final_time-initial_time)/(time_resol);
    double alpha=step_time/(step_x*step_x);
    //cout<<alpha<<endl;
    double* x=new double[n+1];
    double* t=new double[n+1];
    for(int i=0; i<n+1;i++){
        x[i]=i*step_x;
        t[i]=i*step_time;
    }


    double* u=new double[n+1];
    double* v=new double[n+1];

    //set initial conditions
    for(int i=0;i<n;i++){
        u[i]=v[i]=0;
    }
    u[n]=v[n]=1;


    double a,b, c;
    a=c=-alpha;
    b=1+2*alpha;

    string fileout = "backward2.txt";
    ofile.open(fileout);

    for(double time=initial_time;time<=final_time;time+=step_time){

        tridiagonal_solver(a,b,c,n+1,u,v,alpha);
            //cout<<time<<endl;
        u[0]=0;
        u[n]=1;
        for(int i=1; i<n+1;i++){
            v[i]=u[i];
        }
            ofile << setiosflags(ios::showpoint | ios::uppercase);
           // ofile << " time:              x:           u(x): " << endl;
            for (int i=0;i<n+1;i++) {
            ofile << setw(15) << setprecision(8) << x[i];
            ofile << setw(15) << setprecision(8) << u[i] << endl;
            }


    }
    ofile.close();

    delete[] u;
    delete[] v;
    delete[] x;
    delete[] t;

    return 0;
}

void tridiagonal_solver(double A, double B, double C, int n, double* S, double* Q, double cost){

    // setting the value of our problem
    double *d = new double[n];
    d[0] = B; // to avoid division by zero in the next loop for
    double *q_new = new double[n + 1];
    q_new[0] = 0;

    // making a new diagonal elements
    for(int i=1; i<n-1; i++){
        d[i] = B - A*C / d[i-1];
        //cout<<Q[i]<<endl;
        q_new[i] = Q[i] - A*q_new[i-1] / d[i-1];
    }
    d[n-1]=B-(A)*C/d[n-2];
    q_new[n-1] = Q[n-1] - (A)*q_new[n-2] / d[n-2]+cost;
    // substitution for finding S
    S[n-1] = q_new[n-1] / d[n-1];
    for(int i=n-2; i>=0; i--){
        S[i] = ( q_new[i] - C*S[i+1] )/ d[i];
    }

    delete[] d;
    delete[] q_new;

}

