#ifndef DATA_ANALYSIS_H
#define DATA_ANALYSIS_H
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <ios>

using namespace std;

class data_analysis
{
public:
    data_analysis();

    std::string name;
    std::ofstream file;

    void open_file(string name_file);
    void close_file();
    void write(double v);
    void add_column();
    void new_row();
    void printscreen(double aa, std::string name_file);
};

#endif // DATA_ANALYSIS_H
