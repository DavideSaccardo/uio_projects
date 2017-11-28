#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <sstream>
#include "write_to_file.h";

using namespace std;
ofstream ofile;

// print matrix
void print_matrix(double ** A, int n){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout << A[i][j] << "\t";
        }
        cout << "\n";
    }
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


void write_to_file(int n, double ** A, char * name){


    // interacting

    ofile.open(name);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i=0; i < n; i++) {
        for(int j=0; j < n; j++){
            ofile << setw(15) << setprecision(8) << A[i][j] << ", \t";
        }
        ofile << "\n";
    }

    ofile.close();

}




