#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath> //libreria migliore
#include <string>

using namespace std;
ofstream ofile;

#include "write_to_file.h"



void write_to_file(int n, double ** A, string name){


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




