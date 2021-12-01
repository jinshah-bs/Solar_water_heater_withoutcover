//
// Created by thinkcomputational on 25/11/21.
//

#ifndef SOLAR_ANALYSIS_READ_DATA_H
#define SOLAR_ANALYSIS_READ_DATA_H
#include <vector>
extern std::vector <std::vector<double>> MetData;
extern double dt;

std::vector<std::vector<double>> parse2DCsvFile(std::string inputFileName);
double InterpolateWhetherData(int loc, int prop, int iter);
#endif //SOLAR_ANALYSIS_READ_DATA_H
