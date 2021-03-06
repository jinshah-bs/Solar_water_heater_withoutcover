//
// Created by thinkcomputational on 25/11/21.
//

#include "mat_prop.h"
#include <cmath>
#include <vector>
double calcTemp(double h) {
    double T = 6.8539e-17 * pow(h,3.0) - 2.0833e-10 * pow(h,2.0) + 5.436e-4 * h - 16.6521;
    return T;
}

std::vector<double> EnthalpyToTemp(std::vector<double> h) {
    std::vector<double>T;
    for (auto i:h){
        T.push_back(calcTemp(i));
    }
    return T;
}

std::vector<double> dencity(std::vector<double> T) {
    std::vector<double>rho;
    for (auto i:T){
        rho.push_back(calcDencity(i));
    }
    return rho;
}

double calcDencity(double T) {
    return -8.2737e-7 * pow(T, 3.0) + 1.405e-4 * pow(T,2.0) - 0.67915 * T + 885.1784;
}

std::vector<double> spHeat(std::vector<double> T) {
    std::vector<double>spheat;
    for (auto i:T){
        spheat.push_back(calcSpHeat(i));
    }
    return spheat;
}

double calcSpHeat(double T) {
    return 2.00854e-6 * pow(T, 3.0) - 1.01954e-3 * pow(T,2.0) + 3.67468 * T + 1832.244;
}

std::vector<double> theConductivity(std::vector<double> T) {
    std::vector<double>cond;
    for (auto i:T){
        cond.push_back(calcConductivity(i));
    }
    return cond;
}

double calcConductivity(double T) {
    return 4.31257e-12 * pow(T, 3.0) - 4.98256e-9 * pow(T,2.0) - 1.15789e-4 * T + 0.130732;
}

void updateThermalProp(std::vector<double> &rho,
                       std::vector<double> &K,
                       std::vector<double> &Cp,
                       std::vector<double> T) {
    rho = dencity(T);
    K = theConductivity(T);
    Cp = spHeat(T);
}

double mu(double T) {
    double sum = 119.1*exp(-0.07702*T) + 23.67* exp(-0.02033*T);
    return sum*0.001;
}

double TempToEnthalpy(double T) {
    double h = 1.8436e-4 *pow(T, 3.0) + 1.645175*pow(T, 2.0) + 1857.6064*T + 31438.1971;
    return h;
}
