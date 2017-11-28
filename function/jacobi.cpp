#include "jacobi.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include "time.h"

using namespace std;

// function which orders the eigenvalues inside eig_v and eigenvectors inside R
void order_eigpairs(int n, double* eig_v, double** R){
    int it=0;
    double a, b;
    double a_c[n-1], b_c[n-1];
    while(it<n-1){
        for(int i=1; i<n-1; i++){
            a = eig_v[i-1];
            b = eig_v[i];
            for(int r=0; r<n-1; r++){
                a_c[r] = R[r][i-1];
                b_c[r] = R[r][i];
            }
            if(b < a){
                eig_v[i] = a;
                eig_v[i-1] = b;
                for(int g=0; g<n-1; g++){
                    R[g][i-1] = b_c[g];
                    R[g][i] = a_c[g];
                }
            }
        }
        it++;
    }
}


//function which gets the maximum off-diagonal element of A
double get_max(double ** A, int& k, int& l, int n){
    double max=0;
    for(int i=0; i<n-1; i++){
        for(int j=0; j<n-1; j++){
            if(fabs(A[i][j]) > max && i!=j){
                max = fabs(A[i][j]);
                k = i;
                l = j;
            }
        }
    }
    return max;
}


//define the rotation
// A is the matrix you want to trasform
// R is the matrix of eigenvector
void rotate(double** A, double**R, int k, int l, int n){
    double t=0;
    double tau=0;
    double c=0;
    double s=0;
    tau = (A[l][l]-A[k][k])/(2*A[k][l]);

    if(A[k][l] != 0){
        if(tau >= 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        }
        else {
            t = -1.0/(-tau +sqrt(1.0 + tau*tau));
        }

        c = 1/sqrt(1+t*t);
        s = t*c;
    }
    else{
        c=1;
        s=0;
    }

    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A[k][k];
    a_ll = A[l][l];
    A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
    A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
    A[k][l] = 0.0; // hard-coding non-diagonal elements by hand
    A[l][k] = 0.0; // same here
    for ( int i = 0; i < n; i++ ) {
        if ( i != k && i != l ) {
            a_ik = A[i][k];
            a_il = A[i][l];
            A[i][k] = c*a_ik - s*a_il;
            A[k][i] = A[i][k];
            A[i][l] = c*a_il + s*a_ik;
            A[l][i] = A[i][l];
        }

        //  And finally the new eigenvectors
        r_ik = R[i][k];
        r_il = R[i][l];
        R[i][k] = c*r_ik - s*r_il;
        R[i][l] = c*r_il + s*r_ik;
    }

}


void jacobi(int n, double** A, double** R, double* eig_v, double eps, int method, int& iterations,double& timeused){

    int p,q;
    double max_off=get_max(A,p,q,n-1);

// start time
    clock_t start, finish;
    start = clock();

//function jacobi performs brute-force Jacobi algorithm if method is declared to be 0
    if(method==0){
        while(max_off > eps && iterations < 1e6){ //brute force
            int p, q;
            max_off = get_max(A, p, q, n-1);
            rotate(A, R, p, q, n-1);
            iterations++;
        }
    }

//function jacobi performs brute-force Jacobi algorithm if method is declared to be 0
    if(method==1){
        while(max_off> eps && iterations < 1e6){ //cycling
            for(int i=1;i<n-1;i++) {
                for(int j=i+1;j<n-2;j++){
                    rotate(A, R, i, j, n-1);
                    iterations++;
                }
            }
        }
    }

    finish = clock();
    timeused = (double) (finish - start)/(CLOCKS_PER_SEC );

//assign the diagonal of A (eigenvalues) to a vector

//cout <<"Eigenvalues: "<<;
    for (int i = 0; i < n-1; i++){
        eig_v[i] = A[i][i];
      //cout << eig_v[i] << endl; //check results
    }

//function to order eig_v (eigenvalues) and R (eigenvectors)
    order_eigpairs(n, eig_v, R);
}
