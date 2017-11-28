#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath> //libreria migliore
#include "time.h"
#include <armadillo>

using namespace std;
using namespace arma;
ofstream ofile;

double f(double x) {return 100*exp(-10*x);}
double solution(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

int main(int argc, char* argv[])
{

    // Declaration of initial variables:
    char *outfilename;
    int n;

    if (argc<=2){
        cout << "bad usage! you need the number of points n and the outputfile" << endl;
    }
    else{
        n = atoi(argv[1]);
        outfilename = argv[2]; // if you want to import data in matlab use the name
    }                          // A for n=10, B for n=100, C for n=1000

    // variable of the problem
    // A v = q
    // x is the variable of q(x[i])
    // the solution start with the point x=h and finish  when x = 1.0/(1.0 + n) * (n+1)
    double h = 1.0/(1.0 + n);
    double *x = new double[n+2];
    double *q = new double[n+1];
    q[0] = 0;
    double *s = new double[n+1];
    // points of integration
    for(int i=0; i<=n+1; i++){
        x[i] = h*i;
    }

    // A = b[1] c[1] 0 ...
    //     a[1] b[2] c[2]  0
    //     0    a[2] b[3] c[3]

    double *b = new double[n + 1];
    double *c = new double[n + 1];
    double *a = new double[n + 1];

    // setting the last component of non-diagonal elements
    b[0] = 0;
    a[0] = 0;
    a[n] = 0;
    c[0] = 0;
    c[n] = 0;

    // asking the value of a, b,c vectors
    // for(int i = 1; i<n; i++){
    //
    // printf("what is the value of a[%d] \n",i);
    // scanf("%f", &a[i]);
    // printf("what is the value of c[%d] \n",i);
    // scanf("%f", &c[i]);
    // } the same for b but ..

    // setting the value of our problem
    for(int i=1; i<=n; i++){
        b[i] = 2;
        q[i] = h*h*f(x[i]);
        s[i] = solution(x[i]);
    }

    for(int j=1; j<n; j++){
        a[j] = -1;
        c[j] = -1;
    }

    double *d = new double[n + 1];
    d[0] = 1; // to avoid division by zero in the next loop for
    double *q_new = new double[n + 1];
    q_new[0] = 0;

    // Start timing
    clock_t start, finish;
    start = clock();

    // making a new diagonal elements
    for(int i=1; i<=n; i++){
        d[i] = b[i] - a[i-1]*c[i-1] / d[i-1];
        q_new[i] = q[i] - a[i-1]*q_new[i-1] / d[i-1];
    }
    // substitution for finding v
    double *v = new double[n+1];
    v[0] = 0;
    v[n] = q_new[n] / d[n];
    for(int i=n-1; i>=1; i--){
        v[i] = ( q_new[i] - c[i] * v[i+1] )/ d[i];
    }

    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
    // print time execution
    cout << setprecision(10) << setw(20) << "Time used  for norm computation=" << timeused  << endl;

    // print the result
    /*
    for (int i=1; i<=n; i++){
        cout << v[i] << endl;
        cout << s[i] << endl;
        printf("\n");
    }
    */


    // Open file and write results to file:
    ofile.open(outfilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << " x:       s(x):       v(x): " << endl;
    for (int i=1;i<=n;i++) {
        ofile << setw(15) << setprecision(8) << x[i] << ",";
        ofile << setw(15) << setprecision(8) << s[i] << ",";
        ofile << setw(15) << setprecision(8) << v[i] << "," << endl;
        }

    ofile.close();

    delete [] x;
    delete [] q;
    delete [] q_new;
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] v;
    delete [] s;


return 0;
}
