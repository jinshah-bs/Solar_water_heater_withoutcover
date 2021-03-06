//
// Created by thinkcomputational on 01/12/21.
//

#ifndef SOLAR_ANALYSIS_WRITE_DATA_H
#define SOLAR_ANALYSIS_WRITE_DATA_H
#include <vector>
#include <string>
void write(std::vector<double> &data, int time, const std::string filename);
void write(std::vector<std::vector<double>> &data,
           std::vector<double>&time,
           const std::string filename);
void write(std::vector<std::vector<double>> &data, const std::string filename);
#endif //SOLAR_ANALYSIS_WRITE_DATA_H
