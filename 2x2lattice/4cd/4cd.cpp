#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <string>
#include <iostream>
#include "data_analysis.h"

using namespace std;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}


int main(int argc, char* argv[])
{

    data_analysis* dat = new data_analysis();
    int lattice_dim=20;

    // Temperature initial condition
    double T_initial = 1.0;
    double T_final=3.8;
    double T_step=0.7;
    int cut_off=50000;
    int MCcycles=1e6;

    double E = 0;
    double M = 0;

    double** spin_matrix;
    spin_matrix = new double*[lattice_dim];

    for(int i=0; i<lattice_dim; i++){
        spin_matrix[i] = new double[lattice_dim];
    }

    // random generator
    std::random_device rd;
    std::mt19937_64 gen(rd());
    // Then set up the uniform distribution for x \in [[0, 1]
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    // random initial conditions for the lattice
    for(int x=0; x<lattice_dim; x++){
        for(int y=0; y<lattice_dim; y++){

            //if((distribution(gen)<0.5)){
             //   spin_matrix[x][y] =+1;
            //}else
                spin_matrix[x][y] =-1;


            M += (double) spin_matrix[x][y];
        }
    }

    for(int x=0; x<lattice_dim; x++){
        for(int y=0; y<lattice_dim; y++){
            E -= (double) spin_matrix[x][y]*
                    (spin_matrix[periodic(x,lattice_dim,-1)][y] +
                    spin_matrix[x][periodic(y,lattice_dim,-1)]);

        }
    }


    double EnergyDifference[5];
    // variable for the different expectation values
    double local_expectation_values[5];
    double total_expectation_values[5];
    // numeber of accepted configurations
    int accepted;
    // histogram array
    int histArray[lattice_dim*lattice_dim+1];
    // possible energy differences
    // k_B = 1


    // summerize
    dat->printscreen(E, "initial enegry");
    dat->printscreen((double) lattice_dim,"Lattice dimension LxL with L = ");
    dat->printscreen((double) MCcycles,"MC cycles ");
    dat->printscreen(T_initial, "Initial temperature");
    dat->printscreen(T_final, "Final temperature");
    dat->printscreen(T_step, "Temperature step size");




    for(double T=T_initial; T<=T_final; T += T_step ){

        dat->open_file(std::to_string(MCcycles-cut_off)+"_"
                       +std::to_string(T)+"_"
                       +std::to_string(lattice_dim)+".txt");

        for( int de =-8; de <= 8; de+=4) EnergyDifference[de/4 + 2] = exp(-de/T);


        accepted=0;
        for(int i = 0; i < lattice_dim*lattice_dim+1; i++)
        {
            histArray[i] = 0;
        }

        for(int i=0; i<5; i++){
            local_expectation_values[i]=0;
            total_expectation_values[i]=0;
        }

        for (int cycles = 1; cycles <= MCcycles; cycles++){
            double norm=1.0/(cycles);


            for(int x=0; x<lattice_dim; x++){
                for(int y=0; y<lattice_dim; y++){
                    // random position in the lattice
                    int xr = (int) (distribution(gen) * (double)lattice_dim);
                    int yr = (int) (distribution(gen) * (double)lattice_dim);

                    // energy difference from the previous configurations
                    int deltaE;
                    deltaE =(int) 2*spin_matrix[xr][yr]*(
                                spin_matrix[xr][periodic(yr,lattice_dim,-1)] +
                            spin_matrix[periodic(xr,lattice_dim,-1)][yr] +
                            spin_matrix[xr][periodic(yr,lattice_dim,1)] +
                            spin_matrix[periodic(xr,lattice_dim,1)][yr] );

                    // metropolis algorithm with A(i->j) = min(1, exp(-deltaE/beta))
                    if(distribution(gen) <= EnergyDifference[deltaE/4 + 2]){
                        spin_matrix[xr][yr] *= -1;
                        // updating values
                        E += (double) deltaE;
                        M += (double) 2*spin_matrix[xr][yr];
                        accepted++;
                    }

                }
            }

            local_expectation_values[0] += E;
            local_expectation_values[1] += E*E;
            local_expectation_values[2] += fabs(M);
            local_expectation_values[3] += M*M;
            local_expectation_values[4] += M;

            dat->write((double) cycles);
            dat->add_column();
            dat->write(local_expectation_values[0]*norm/lattice_dim*lattice_dim);
            dat->add_column();
            dat->write((local_expectation_values[1]-
                       local_expectation_values[0]*norm*local_expectation_values[0])*norm
                    /(lattice_dim*lattice_dim*T*T));
            dat->add_column();
            dat->write(local_expectation_values[2]*norm/lattice_dim*lattice_dim);
            dat->add_column();
            dat->write( (local_expectation_values[3]-
                        local_expectation_values[2]*norm*local_expectation_values[2])*norm
                    /(lattice_dim*lattice_dim*T));
            dat->add_column();
            dat->write(local_expectation_values[4]*norm/lattice_dim*lattice_dim);
            dat->add_column();
            dat->write(accepted);
            dat->new_row();

            if(cycles>cut_off){
                histArray[(int) (E+2*lattice_dim*lattice_dim)/4]++;
            }
        } // MC cycle loop
        dat->close_file();

        dat->open_file("hist_"+std::to_string(T)+".txt");
        for(int i=0; i<=lattice_dim*lattice_dim; i++){
            dat->write((double) 4*i-2*lattice_dim*lattice_dim);
            dat->add_column();
            dat->write((double) histArray[i]/((double) MCcycles-cut_off));
            dat->new_row();
        }
        dat->close_file();

    } // Temperature loop





    for(int i=0; i<lattice_dim; i++) delete[] spin_matrix[i];
    delete[] spin_matrix;
    delete dat;
    return 0;
}



/*
double z;
double exp_meanE=0;
double exp_meanE_2=0;
double exp_meanM =0;
double exp_meanM_2 =0;

z = 12 + 2*exp(8/T_final) - 2*exp(-8/T_final);

exp_meanE = 2*8*exp(-8/T_final) + 2*(-8)*exp(8/T_final);
exp_meanE_2 = 2*8*8*exp(-8/T_final) + 2*(-8)*(-8)*exp(8/T_final);
// mean of absulute value of M
exp_meanM =2*4*exp(+8/T_final)  + 2*4*2 ;//-2*4
//four deg
exp_meanM_2 = 4*4*exp(+8/T_final) + (-4)*(-4)*exp(+8/T_final)  + 4*(2)*(2) + 4*(-2)*(-2);

dat->printscreen(meanM,"calculated <|M|> ");
dat->printscreen( exp_meanM/z ,"theoric <|M|> " );
printf("\n");
dat->printscreen( meanE ,"calculated <E> " );
dat->printscreen( exp_meanE/z ,"theoric <E> " );
printf("\n");
dat->printscreen( meanM_2/2 ,"calculated " );
dat->printscreen( exp_meanM_2/z ,"theoric " );
printf("\n");
dat->printscreen( meanE_2 ,"calculated " );
dat->printscreen( exp_meanE_2/z ,"theoric " );
*/
