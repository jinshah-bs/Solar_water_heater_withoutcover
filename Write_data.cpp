//
// Created by thinkcomputational on 01/12/21.
//

#include "Write_data.h"
#include <iostream>
#include <fstream>
#include <sstream>
void write(std::vector<double> &data, int time, const std::string filename) {
    std::string fNAME, folder;
    std::ofstream files;
    std::stringstream min;
    unsigned int N = data.size();
    min << time;
    folder = "mkdir -p Results/" + min.str();
    system(folder.c_str()); // Creates the directory
    fNAME = "Results/" + min.str() +  "/" + filename + ".txt";
    files.open(fNAME.c_str(), std::ios::out);
    for (auto item: data) {
        files << item << std::endl;
    }
    files.close();
}

void write(std::vector<std::vector<double>> &data, std::vector<double> &time, const std::string filename) {
    std::string folder, Fname;
    std::ofstream files;
    folder = "mkdir -p Results";
    system(folder.c_str());
    Fname = "Results/" + filename + ".csv";
    unsigned int N = time.size();
    files.open(Fname.c_str(), std::ios::out);
    for (auto t:time) {
        files << "t= " << t << " , ";
    }
    files << std::endl;
    for (int i = 0; i < data[0].size(); ++i) {
        for (int j = 0; j < N; ++j) {
            files << data[j][i] << " , ";
        }
        files <<std::endl;
    }
    files.close();
}

void write(std::vector<std::vector<double>> &data, const std::string filename) {
    std::string folder, Fname;
    std::ofstream files;
    folder = "mkdir -p Results";
    system(folder.c_str());
    Fname = "Results/" + filename + ".csv";
    files.open(Fname.c_str(), std::ios::out);
    files << "Time" << " , " <<"q_dni" << " , " << "Wind_vel" << " , " << "T_ambient" << " , "
          << "Q_loss" <<  " , " << "DeltaH" << " , " << "Eff"
          << " , " << "T_out" << " , " << "T_in" << std::endl;

    for (auto row: data) {
        for (auto col : row) {
            files << col << " , ";
        }
        files <<std::endl;
    }

}
