#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath>
#include "time.h"
#include <armadillo>


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
//            if(eps < A[i][j]){
//                k = i;
//                l = j;
//                return A[i][j];
//            }
            if(fabs(A[i][j]) > max){
                max = fabs(A[i][j]); // always positive
                k = i;
                l = j;
            }
        }
    }
    return max;
}

//define the rotation
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
// frob norm
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


// main function

int main()//int argc, char *argv[])
{
    int n;
    // 30 no possible
    int loop_r=2; // number of times the loop will aplly the algorithm
    /*
    if (argc<=1){
        cout << "Bad usage! you need the number of points n" << endl;
    }
    else{
        n = atoi(argv[1]);
    }
    */

    int p_p; // point for fixed precision

    double* repeat = new double[loop_r];
    double* point_numb = new double[loop_r];
    double* t_elapsed = new double[loop_r];

    for(int ii = 1; ii < loop_r; ii++){

        n=400*ii;
        double eig_v[n-1];
        double ** R = new double*[n-1]; // initialize matrix of eigenvectors
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
        // diui + ei−1ui−1 + ei+1ui+1 = λui,
        // conditions on indipendent variable rho
        double rho_0 = 0.0;
        double rho_n = 11; // aproximation rho max= infinity
        double h = (rho_n - rho_0)/n;
        double *rho = new double[n+1];
        rho[0] = rho_0;
        rho[n] = rho_n;
        for(int i=1; i<n; i++){
            rho[i] = i*h; // rho[i] = rho_0 + i*h & rho[n] = rho_n
        }

        // define the matrix of our problem
        // define the potential V
        double *V = new double[n];
        double omega = 0.25;
        V[0] = rho[0];
        for(int i = 1; i<n; i++){
            V[i] = rho[i]*rho[i]*omega*omega + 1/rho[i];   //rho[i]*rho[i];//
        }
        // diagonal elements
        double *d = new double[n+1];
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
        /*
        double f_i; // initial frobenius norm
        f_i = frob_norm(A, n-1);
        */
        // Start timing
        clock_t start, finish;
        start = clock();

        // loop that diagonalizes A
        while(max_off > eps && iterations < 1e6){
            int p, q;
            max_off = off_norm(A, p, q, n-1);
            rotate(A, R, p, q, n-1);
            iterations++;
        }

        finish = clock();
        double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
        // print time execution
        // cout << setprecision(10) << setw(20) << "Time used =" << timeused  << endl;
        t_elapsed[ii] = timeused;

        //
        cout <<  n << endl;


        // diagonal elements
        // cout << "diagonal elements" << endl;
        // eigenvalue
        //order(eig_v, n-1);

        for (int i = 0; i < n-1; i++){
            eig_v[i] = A[i][i];
           //cout << eig_v[i] << endl;
        }
        // order(eig_v, n-1);

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


        for (int i = 0; i < n-1; i++){

          // cout << eig_v[i] << endl;
        }

        // print eig_v
        if(ii==loop_r-1){
            ofile.open("rho_eig_v.txt");
            ofile << setiosflags(ios::showpoint | ios::uppercase);
            ofile << "rho:   eigenvalue" << endl;
            for (int i=0; i < n-1; i++) {
                ofile << setw(15) << setprecision(8) << rho[i] << ", \t";
                ofile << setw(15) << setprecision(8) << eig_v[i] << "," << endl;
            }

            ofile.close();

            ofile.open("eigenkets.txt");
            ofile << setiosflags(ios::showpoint | ios::uppercase);
            ofile << "eigenkets in coloumn" << endl;
            for (int i=0; i < n-1; i++) {
                for(int j=0; j < n-1; j++){
                    ofile << setw(15) << setprecision(8) << R[i][j] << ", \t";
                }
                ofile << "\n";
            }

            ofile.close();
        }

        // print_matrix(A, n);
        // print_matrix(R, n);

        /*
        double f_f; // final frobenius norm
        f_f = frob_norm(A, n-1);
        cout << "the difference between the initial and final frobenius norm is " << fabs(f_f - f_i) << endl;
        */
        // free space


        if(ii==loop_r-1){
            for (int i = 0; i < n; i++){
                delete[] A[i];
                delete[] R[i];
            }

            delete[] A;
            delete[] R;
            delete[] rho;
            delete[] V;
            delete[] d;
        }

        // variables of iterations

        repeat[ii] = iterations;
        point_numb[ii] = n;

        // how many points n we need to have
        // in order to get the lowest three
        // eigenvalues with approximately four leading digits
        // after the decimal point
        /*
        double first;
        first = floor(A[0][0]*1E4)/1E4;
        double second;
        second = floor(A[1][1]*1E4)/1E4;
        double third;
        third = floor(A[2][2]*1E4)/1E4;

        if(first == 3){
            if(second == 11){
                if(third == 19){
                    p_p = n;
                }
            }
        }
        */

    }


    // write to file txt
    ofile.open("prova.txt");
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << " numb:       repeat:      time:" << endl;
    for (int i=0; i < loop_r; i++) {
        ofile << setw(15) << setprecision(8) << point_numb[i] << ", \t";
        ofile << setw(15) << setprecision(8) << repeat[i] << ", \t";
        ofile << setw(15) << setprecision(8) << t_elapsed[i] << "," << endl;
    }

    ofile.close();

    // cout << p_p << endl;


    // armadillo version to
    /*
    mat H(n-1,n-1);
    for(int i=0; i<n-1; i++){
        for(int j=0; j<n-1; j++){
            if(i==j){
                H(i,j) = d[i+1];;
            }
            else if (fabs(i-j)==1) {
                H(i,j) = e;
            }
            else{
                H(i,j) = 0;
            }
        }
    }
    vec Eigval(n-1);
    eig_sym(Eigval, H);
    Eigval.print();
    */

    delete[] t_elapsed;
    delete[] repeat;
    delete[] point_numb;

    return 0;

}
