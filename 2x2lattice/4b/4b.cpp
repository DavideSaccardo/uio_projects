#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <string>
#include <iostream>
#include "data_analysis.h"

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}


int main(int argc, char* argv[])
{

    int lattice_dim=2;
    double E = 0;
    double M = 0;

    data_analysis* dat = new data_analysis();
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

           // if((distribution(gen)<0.5)){
                spin_matrix[x][y] =+1;
           // }else spin_matrix[x][y] =-1;


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

    // Temperature initial condition
    double T_initial = 1;
    double T_final=1;
    double T_step=1.4;

    double EnergyDifference[5];
    // possible energy differences
    // k_B = 1

    // variable for the different expectation values
    double local_expectation_values[5];
    // numeber of accepted configurations
    int accepted;
    // histogram array
    int histArray[lattice_dim*lattice_dim+1];


    int MCcycles=10000000;
    // summerize
    dat->printscreen(E, "initial enegry");
    dat->printscreen((double) lattice_dim,"Lattice dimension LxL with L = ");
    dat->printscreen((double) MCcycles,"MC cycles ");
    dat->printscreen(T_initial, "Initial temperature");
    dat->printscreen(T_final, "Final temperature");
    dat->printscreen(T_step, "Temperature step size");

    dat->open_file("dT1_L2_up.txt");
    for(double T=T_initial; T<=T_final; T += T_step ){

        for( int de =-8; de <= 8; de+=4) EnergyDifference[de/4 + 2] = exp(-de/T);

        accepted=0;

        for(int i = 0; i < lattice_dim*lattice_dim+1; i++)
        {
            histArray[i] = 0;
        }

        for(int i=0; i<5; i++){
            local_expectation_values[i]=0;
        }

        for (int cycles = 1; cycles <= MCcycles; cycles++){

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
                        //                        accepted++;
                        //                        if(cycles>100000){
                        //                            histArray[(int) (E+2*lattice_dim*lattice_dim)/4]++;
                        //                        }

                    }
                }
            }

            local_expectation_values[0] += E;
            local_expectation_values[1] += E*E;
            local_expectation_values[2] += fabs(M);
            local_expectation_values[3] += M*M;
            local_expectation_values[4] += M;

            if(cycles%100==0){
            dat->write(cycles);
            dat->add_column();

            dat->write(local_expectation_values[0]/cycles);
            dat->add_column();
            dat->write((local_expectation_values[1]-
                       local_expectation_values[0]*local_expectation_values[0]/cycles
                    )/(T*T*cycles));
            dat->add_column();
            dat->write(local_expectation_values[2]/cycles);
            dat->add_column();
            dat->write((local_expectation_values[3]-
                       local_expectation_values[4]*local_expectation_values[4]/cycles
                    )/(T*cycles));
            dat->add_column();
            dat->write(local_expectation_values[4]/cycles);
            dat->add_column();
            dat->write((local_expectation_values[3]-
                       local_expectation_values[2]*local_expectation_values[2]/cycles
                    )/(T*cycles));
            dat->add_column();
            dat->new_row();
            }

        } // MC cycle loop

    } // Temperature loop
    dat->close_file();


    //    dat->open_file("hist3_5.txt");
    //    for(int i=0; i<=lattice_dim*lattice_dim; i++){
    //        dat->write((double) 4*i-2*lattice_dim*lattice_dim);
    //        dat->add_column();
    //        dat->write(histArray[i]/(MCcycles-100000));
    //        dat->new_row();
    //    }
    //    dat->close_file();


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

dat->printscreen(local_expectation_values[2]/MCcycles,"calculated <|M|> ");
dat->printscreen( exp_meanM/z ,"theoric <|M|> " );
printf("\n");
dat->printscreen( local_expectation_values[0]/MCcycles ,"calculated <E> " );
dat->printscreen( exp_meanE/z ,"theoric <E> " );
printf("\n");
dat->printscreen( local_expectation_values[3]/MCcycles ,"calculated <M*M>" );
dat->printscreen( exp_meanM_2/z ,"theoric <M*M> " );
printf("\n");
dat->printscreen( local_expectation_values[1]/MCcycles ,"calculated <E*E>" );
dat->printscreen( exp_meanE_2/z ,"theoric <M*M>" );
*/
