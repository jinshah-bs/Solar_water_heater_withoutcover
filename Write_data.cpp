//
// Created by thinkcomputational on 01/12/21.
//

#include "Write_data.h"
#include <iostream>
#include <fstream>
#include <sstream>
void write(std::vector<double> data, int time, const std::string filename) {
    std::string fNAME, folder;
    std::ofstream files;
    std::stringstream min;
    unsigned int N = data.size();
    min << time*5;
    folder = "mkdir - p Results/" + min.str();
    system(folder.c_str()); // Creates the directory
    fNAME = "Results" + min.str() + filename + ".csv";
    files.open(fNAME.c_str(), std::ios::out);
    for (auto item: data) {
        files << item << std::endl;
    }
    files.close();
}
