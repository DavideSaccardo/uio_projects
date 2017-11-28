#include "data_analysis.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <cmath> //libreria migliore
#include <string>

using namespace std;
ofstream ofile;


data_analysis::data_analysis()
{

}

void data_analysis::open_file(std::string name_file){
    name = name_file;
    file.open(name, std::ios::out);

}
void data_analysis::close_file(){
    file.close();
}

void data_analysis::printscreen(double aa, std::string name_file){
    cout <<name_file << aa << endl;
}

void data_analysis::write(double v){
            file << setw(15)  << v << ",";
}

void data_analysis::add_column(){
    file << "\t";
}
void data_analysis::new_row(){
    file << "\n";
}
