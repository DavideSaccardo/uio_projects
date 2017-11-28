#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <string>
#include <iostream>
#include "data_analysis.h"
#include <mpi.h>

using namespace std;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
    return (i+limit+add) % (limit);
}

// This is the main function used to find:
// mean energy
// mean absulute magnetic moment
// susceptibility
// heat capacity
// mean magnetic moment
// as a function of the temperature
// Boltzmann constant k_B = 1
int main(int argc, char* argv[])
{

    // data_analysis class is usful to write to file
    // and to print to screen
    data_analysis* dat = new data_analysis();

    // MPI initialization
    // total number of processors
    int numprocs;
    // processor's rank
    int my_rank;
    // when we run the program this gets
    // total number of processors we want to use
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    // this loop is used to generate data for
    // different lattice dimensions L
    for(int lattice_dim=40; lattice_dim<=100; lattice_dim +=20){

        // we can measure the running time per each lattice dimension
        // here we memorize the time start
        double timestart = MPI_Wtime();

        // define the temperature range
        double T_initial = 2.2;
        double T_final=2.3;
        // temperature step size
        double T_step=0.005;

        // the next constant can be used to start
        // writing to file data after a
        // number cut_off of Monte Marlo cycles
        int cut_off=0;

        // this is the total number of Monte Carlo
        // cycles we want to perform
        int MCcycles=1e3;

        // we calculate the normalization factor
        // considering the number of processors,
        // Monte Carlo cycles and cut_off
        double norm=1.0/((MCcycles-cut_off)*numprocs);

        // we write to file opening the file with
        // the processor of rank=0 writes to file
        if(my_rank==0){
            // the name of the file gaves the main information about data
            dat->open_file(std::to_string(MCcycles-cut_off)+"_"
                           +std::to_string(T_initial)+"to"+std::to_string(T_final)+"_"
                           +std::to_string(lattice_dim)+".txt");
        }


        // variables that are used to keep track of
        // the actual configuration of the system
        // energy configuration
        double E = 0;
        // magnetization configuration
        double M = 0;
        // lattice configuration
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
                if((distribution(gen)<0.5)){
                    spin_matrix[x][y] =+1;
                }else spin_matrix[x][y] =-1;
                M += (double) spin_matrix[x][y];
            }
        }

        // initial energy configuration
        for(int x=0; x<lattice_dim; x++){
            for(int y=0; y<lattice_dim; y++){
                E -= (double) spin_matrix[x][y]*
                        (spin_matrix[periodic(x,lattice_dim,-1)][y] +
                        spin_matrix[x][periodic(y,lattice_dim,-1)]);

            }
        }

        // vector to precalculate the exp(- delta energy/ (k_b T))
        double EnergyDifference[5];
        // variable for the different expectation values
        // local is used by each processor
        double local_expectation_values[5];
        // total is used to gather together all the results
        // of different processors
        double total_expectation_values[5];

        // more variables
        // number of accepted configurations
        int accepted;
        // histogram array that counts how many times
        // an energy appears
        int histArray[lattice_dim*lattice_dim+1];

        // summerize what we are doing
        if(my_rank==0){
            dat->printscreen(E, "initial enegry");
            dat->printscreen((double) lattice_dim,"Lattice dimension LxL with L = ");
            dat->printscreen((double) MCcycles,"MC cycles ");
            dat->printscreen(T_initial, "Initial temperature");
            dat->printscreen(T_final, "Final temperature");
            dat->printscreen(T_step, "Temperature step size");
        }

        // loop over the temperature range
        for(double T=T_initial; T<=T_final; T += T_step ){

            // precalculate the exponential
            for( int de =-8; de <= 8; de+=4) EnergyDifference[de/4 + 2] = exp(-de/T);

            /*
             * intialization of variables
            accepted=0;
            for(int i = 0; i < lattice_dim*lattice_dim+1; i++)
            {
                histArray[i] = 0;
            }
            */

            //initialize per each temperature
            for(int i=0; i<5; i++){
                local_expectation_values[i]=0;
                total_expectation_values[i]=0;
            }

            // Monte Carlo experiments
            for (int cycles = 1; cycles <= MCcycles; cycles++){

                // loop over all spins
                for(int x=0; x<lattice_dim; x++){
                    for(int y=0; y<lattice_dim; y++){

                        //take a random position in the lattice
                        int xr = (int) (distribution(gen) * (double)lattice_dim);
                        int yr = (int) (distribution(gen) * (double)lattice_dim);

                        // energy difference from the previous configurations
                        int deltaE;
                        deltaE =(int) 2*spin_matrix[xr][yr]*(
                                    spin_matrix[xr][periodic(yr,lattice_dim,-1)] +
                                spin_matrix[periodic(xr,lattice_dim,-1)][yr] +
                                spin_matrix[xr][periodic(yr,lattice_dim,1)] +
                                spin_matrix[periodic(xr,lattice_dim,1)][yr] );

                        // Metropolis algorithm to choose if we accpet or not the flip
                        if(distribution(gen) <= EnergyDifference[deltaE/4 + 2]){
                            spin_matrix[xr][yr] *= -1;

                            // updating actual configurations
                            E += (double) deltaE;
                            M += (double) 2*spin_matrix[xr][yr];
                            //accepted++;
                        }
                        //histArray[(int) (E+2*lattice_dim*lattice_dim)/4]++;

                    } // end y loop
                } // end x loop

                // updating interested variables at the end of each
                // monte carlo experiments
                if(cycles>cut_off){
                    local_expectation_values[0] += E;
                    local_expectation_values[1] += E*E;
                    local_expectation_values[2] += fabs(M);
                    local_expectation_values[3] += M*M;
                    local_expectation_values[4] += M;
                }

            } // MC cycle loop

            // the following loop gathers the results from different
            // processors into the total_expectation_values vector
            for(int i=0; i<5; i++){
                MPI_Reduce(&local_expectation_values[i], &total_expectation_values[i]
                           , 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            }

            // write to file all the results
            if(my_rank==0){
                dat->write((double) T);
                dat->add_column();
                // expectation value of energy
                dat->write(total_expectation_values[0]*norm/(lattice_dim*lattice_dim));
                dat->add_column();
                // heat capacity
                dat->write((total_expectation_values[1]-
                           total_expectation_values[0]*norm*total_expectation_values[0])*norm
                        /(lattice_dim*lattice_dim*T*T));
                dat->add_column();
                // expectation value of absulute value of magnetization
                dat->write(total_expectation_values[2]*norm/(lattice_dim*lattice_dim));
                dat->add_column();
                // susceptibility
                dat->write( (total_expectation_values[3]-
                            total_expectation_values[2]*norm*total_expectation_values[2])*norm
                        /(lattice_dim*lattice_dim*T));
                dat->add_column();
                // expectation value of magnetization
                dat->write(total_expectation_values[4]*norm/(lattice_dim*lattice_dim));
                dat->new_row();
            }

        } // end Temperature loop

        if(my_rank==0) dat->close_file();

        // time used
        double timeend = MPI_Wtime();
        double timeused = (double) (timeend - timestart);
        dat->printscreen((double) (timeend - timestart), "running time");

        for(int i=0; i<lattice_dim; i++) delete[] spin_matrix[i];
        delete[] spin_matrix;

    } // loop lattice

    // free space

    delete dat;


    MPI_Finalize ();

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
