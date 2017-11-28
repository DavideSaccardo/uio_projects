#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <array>
#include <cmath>
#define N 100
using namespace std;
using std::array;
ofstream ofile;

int main()
{
    //set variables
   //int n=1e2;
    int time_resol=1e4;
    double initial_x=0;
    double final_x=1;
    double step_x=(final_x-initial_x)/(N);
    double step_y=step_x;
    double initial_time=0;
    double final_time=10;
    double step_time=(final_time-initial_time)/(time_resol);
    double alpha=step_time/(step_x*step_x);
    //cout<<alpha<<endl;
    array<double, N+1> x;
    array<double, N+1> t;
    for(int i=0; i<N+1;i++){
        x[i]=i*step_x;
        t[i]=i*step_time;
    }

    array<array<double, N+1>, N+1> u{}; //create matrix filled with zeros, we are guessing the values inside u
    array<array<double, N+1>, N+1> v{};
    array<array<double, N+1>, N+1> temp;
    //set initial conditions
    for(int i=0;i<N+1;i++){
        v[i][N]=u[i][N]=1;

    }

//    for (int i=0;i<N+1;i++) {
//        for(int j=0;j<N+1;j++){
//            cout <<setw(15)<< v[i][j];
//            if(j==N) cout<<endl;
//        }
//    }


    string fileout = "testlol.txt";
    ofile.open(fileout);
//    for (int i=0;i<n+1;i++) {
//        ofile << setw(15) << setprecision(8) << x[i]<<endl;
//        }
//    ofile.close();

    double delta;
    int cutoff=1e4;
    double tolerance=1e-6;
   // double difference=tolerance+1;
    int iterations =0;
    double cost=1+4*alpha;
    for(double time=initial_time;time<=final_time;time+=step_time){
        double difference=tolerance+1;
        v=u;
        //cout<<"ehi"<<endl;
        while((iterations<cutoff)&&(difference>tolerance)){
            difference=0;
            temp=u;
            double delta = 0;
            //cout<<"ehi"<<endl;
            for(int i=1;i<N;i++){
                for(int j=1;j<N;j++){
                    delta=alpha*(temp[i+1][j]+temp[i-1][j]+temp[i][j-1]+temp[i][j+1]);
                    u[i][j]=(delta+v[i][j])/(cost);
                    //cout<<u[i][j]<<endl;
                    difference=fabs(u[i][j]-temp[i][j]);
                }
            }
            iterations++;
            if(iterations>cutoff&&difference>tolerance) cout<<"Failure"<<endl;
            difference/=pow(N-1,2.0);
        }
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        // ofile << " time:              x:           u(x): " << endl;
        for (int i=0;i<N+1;i++) {
            for(int j=0;j<N+1;j++){
                ofile << setw(15) << setprecision(8) << u[i][j];
            }
            ofile <<setw(15) << setprecision(8) <<endl;
        }
    }
    //cout<<iterations<<endl;
    ofile.close();

    return 0;
}


