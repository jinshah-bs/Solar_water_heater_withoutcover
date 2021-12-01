//
// Created by thinkcomputational on 25/11/21.
//

#ifndef SOLAR_ANALYSIS_MAT_PROP_H
#define SOLAR_ANALYSIS_MAT_PROP_H
#include <vector>

std::vector<double>EnthalpyToTemp(std::vector<double>h);
double calcTemp(double T);

double TempToEnthalpy(double T);

std::vector<double>dencity(std::vector<double>T);
double calcDencity(double T);

std::vector<double>spHeat(std::vector<double>T);
double calcSpHeat(double T);

std::vector<double>theConductivity(std::vector<double>T);
double calcConductivity(double T);

void updateThermalProp(std::vector<double> &rho,
                       std::vector<double> &K,
                       std::vector<double> &Cp,
                       std::vector<double> T);

double mu(double T);

#endif //SOLAR_ANALYSIS_MAT_PROP_H
