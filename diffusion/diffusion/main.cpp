#include <iostream>
#include <array>
#include "pde_diffusion.h"
#include "data_analysis.h"

using std::array;
using namespace std;

int main()
{
    // define general parameters valid for both
    // one and two dimensions
    // number of points in the spatial coordinate(s)
    //  x: 0 -> 1       y: 0 -> 1
    int N = 100;
    // step size of the spatial coordinates
    double h = 1.0/N;
    // time step size
    double dt = 1e-3;
    // final time
    double t_final = 5;
    // constant of the partial differential equation
    double alpha = dt/(h*h);

    // initial conditions for one dimensional
    double* u = new double[N+1];

    for(int i=0;i<N;i++){
       u[i]=0;
    }
    u[N]=1;


    double** u_2d = new double*[N+1];
    for(int i=0;i<N+1;i++){
       u_2d[i]=new double[N+1];
    }

    for(int i=0; i<N+1;i++){
        for(int j=0; j<N+1;j++){
            u_2d[i][j] = 0;
        }
        u_2d[i][N] = 1;
    }

    data_analysis* dat = new data_analysis;
    dat->open_file("pde.txt");

    dat->printscreen(alpha, "alpha ");
    pde_diffusion* solver = new pde_diffusion(N, alpha);

    double tolerance=1e-6;
    int cutoff=1e4;
    for(int j=0; j<=(int) t_final/dt ; j++){
        solver->two_dimension(u_2d, tolerance, cutoff);

        if(j%(3000)){
            for(int i=0; i<=N; i++){
                for(int j=0; j<=N; j++){
                    dat->write(u_2d[i][j]);
                    dat->add_column();
                }
                dat->new_row();
            }
            dat->new_row();
        }
    }

    dat->close_file();
    delete dat;
    delete solver;
    delete[] u;
    for(int i=0; i<=N+1; i++) delete u[i];
    delete u;







    return 0;
}
