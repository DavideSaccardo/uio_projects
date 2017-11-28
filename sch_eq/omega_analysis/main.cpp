#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include "time.h"
#include <armadillo>
#include <string>
#include <sstream>


using namespace std;
using namespace arma;
ofstream ofile;

// order
void order(double* v, int n){
    double max=v[0];
    int it=0;
    double a, b;
    while(it<n){
        for(int i=1; i<n; i++){
            a = v[i-1];
            b = v[i];
            if(b < a){
                v[i] = a;
                v[i-1] = b;
            }
        }
        it++;
    }
}

// define the off_norm
double off_norm(double ** A, int& k, int& l, int n){
    double max=0;
    double eps = 1e5;
    for(int i=0; i<n; i++){
        for(int j=i+1; j<n; j++){
            if(fabs(A[i][j]) > max){
                max = fabs(A[i][j]); // always positive
                k = i;
                l = j;
            }
        }
    }
    return max;
}

// define the rotation
// A is the matrix you want to trasform
// R will be the matrix of eigenvector
void rotate(double** A, double** R, int k, int l, int n){
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


// print matrix
void print_matrix(double ** A, int n){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout << A[i][j] << "\t";
        }
        cout << "\n";
    }
}

// frobenius norm
double frob_norm(double ** A, int n){
    double sum;
    sum=0;
    for(int i=0; i<n; i++){
        for(int q=0; q<n; q++){
            sum += A[i][q]*A[i][q];
        }
    }
    cout << "frobenius norm" << "\t" << sqrt(sum) << endl;
    return sqrt(sum);
}

// order eigenvalues and eigenvectors
// dimension of R and eig_v is n-1
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







// main function
// This program solve two time-independent Schrodinger's equations (already scaled)
// One characterized by Harmonic Potential
// Another characterized by Harmonic and Coulomb Potentials with different value of omega

int main(int argc, char *argv[])
{
    int n; // number of points of the problem
    // the matrix for solving the eigenvalue problem
    // will be (n-1)X(n-1)

    if (argc<=1){
        cout << "Bad usage! You need the number of points n" << endl;
    }
    else{
        n =atoi(argv[1]);
    }

    // set omega_r of the problem
    double omega[4] = {0.01, 0.5, 1, 5};
    // define the max value of the variable rho
    double rho_n[5] = {50, 7, 4.8, 2.3, 5};

    // define a loop where we will write to file different solutions
    // of the schrodinger's equations
    // four interacting cases V = rho[i]*rho[i]*omega[ii]*omega[ii] + 1/rho[i];
    // and one for the harmonic oscillator
    // V = rho[i]*rho[i]
    for(int ii = 0; ii < 5; ii++){

        // initialize matrix of eigenvectors
        // every columns of R will be an eigenvector
        double eig_v[n-1];
        double ** R = new double*[n-1];
        for(int i = 0; i < n-1; i++){
            R[i] = new double[n-1];
        }

        for(int i=0; i<n-1; i++){
            for(int j=0; j<n-1; j++){
                if(i==j){
                    R[i][j] = 1;
                }
                else{
                    R[i][j] = 0;
                }
            }
        }


        //radial schrodinger's equation with potential V
        // conditions on indipendent variable rho
        double rho_0 = 0.0;
        // aproximation rho max = infinity
        double h = (rho_n[ii] - rho_0)/n;
        double *rho = new double[n+1];
        // we set the first and the last value of rho
        // but they will not be used in the jacobi algorithm
        rho[0] = rho_0;
        rho[n] = rho_n[ii];
        for(int i=1; i<n; i++){
            rho[i] = i*h; // rho[i] = rho_0 + i*h & rho[n] = rho_n
        }

        // define the needed value to build our matrix (n-1)X(n-1)
        // define the potential V
        // we use n-1 elements
        double *V = new double[n];
        V[0] = 0;
        for(int i = 1; i<n; i++){
            if(ii < 4){
                V[i] = rho[i]*rho[i]*omega[ii]*omega[ii] + 1/rho[i];
            }
            else{
                V[i] = rho[i]*rho[i];

            }
        }



        // diagonal elements
        double *d = new double[n];
        d[0] = 0;
        for(int i =1; i<n; i++){
            d[i] = V[i] + 2.0/(h*h) ; //
        }
        // upper and lower diagonal elements
        double e = -1.0/(h*h) ;

        // initialize matrix
        double ** A = new double*[n-1];
        for(int i = 0; i < n-1; i++){
            A[i] = new double[n-1];
        }

        // define the matrix
        for(int i=0; i<n-1; i++){
            for(int j=0; j<n-1; j++){
                if(i==j){
                    A[i][j] = d[i+1];
                }
                else if (fabs(i-j)==1) {
                    A[i][j]  = e;
                }
                else{
                    A[i][j]  = 0;
                }
            }
        }

        int k, l;
        double max_off;
        max_off = off_norm(A, k, l, n-1);
        int iterations = 0;
        double eps = 1.0E-10; // define tolerance

        while(max_off > eps && iterations < 1e6){
            int p, q;
            max_off = off_norm(A, p, q, n-1);
            rotate(A, R, p, q, n-1);
            iterations++;
        }

        // define the eigenvalues
        for (int i = 0; i < n-1; i++){
            eig_v[i] = A[i][i];
            //cout << eig_v[i] << endl;
        }

        // order eigenvalues and eigenvectors
        order_eigpairs(n, eig_v, R);

        // write to a file
        // first column -> variable rho
        // second column -> eigenvalue
        // third column -> first eigenvector
        // fourth column -> second eigenvector
        // and so on
        std::string pi;
        if(ii == 4){    // non interacting
            std::string pi = "rho_val_vec_noninter.txt";
            ofile.open(pi);
            ofile << setiosflags(ios::showpoint | ios::uppercase);
            ofile << "rho:   eigenvalue:    eivec" << endl;
            for (int i=0; i < n-1; i++) {
                ofile << setw(15) << setprecision(8) << rho[i] << ", \t";
                ofile << setw(15) << setprecision(8) << eig_v[i] << ", \t";
                for(int j=0; j < n-1; j++){
                    ofile << setw(15) << setprecision(8) << R[i][j] << ", \t";
                }
                ofile << "\n";
            }

            ofile.close();
        }
        else{   // interacting
            std::string pi = "rho_val_vec" + std::to_string(omega[ii]) + ".txt" ;
            ofile.open(pi);
            ofile << setiosflags(ios::showpoint | ios::uppercase);
            ofile << "rho:   eigenvalue:    eivec" << endl;
            for (int i=0; i < n-1; i++) {
                ofile << setw(15) << setprecision(8) << rho[i] << ", \t";
                ofile << setw(15) << setprecision(8) << eig_v[i] << ", \t";
                for(int j=0; j < n-1; j++){
                    ofile << setw(15) << setprecision(8) << R[i][j] << ", \t";
                }
                ofile << "\n";
            }

            ofile.close();
        }


        // free space
        for (int i = 0; i < n-1; i++){
            delete[] A[i];
            delete[] R[i];
        }

        delete[] A;
        delete[] R;
        delete[] rho;
        delete[] V;
        delete[] d;

    }

    return 0;

}
