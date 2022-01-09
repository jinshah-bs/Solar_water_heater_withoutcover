//
// Created by thinkcomputational on 25/11/21.
//

#include "mat_prop.h"
#include <cmath>
#include <vector>
double calcTemp(double h) {
    return 2204 * pow(h, 1.01) - 8628;
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
    return -8.2737e-7 * pow(T, 3.0) + 1.405e-4 * pow(T,2.0) - 0.67915 * T + 917.2;
}

std::vector<double> spHeat(std::vector<double> T) {
    std::vector<double>spheat;
    for (auto i:T){
        spheat.push_back(calcSpHeat(i));
    }
    return spheat;
}

double calcSpHeat(double T) {
    return (-5.35 * pow(T,-0.7875) + 2.36)*1000;
}

std::vector<double> theConductivity(std::vector<double> T) {
    std::vector<double>cond;
    for (auto i:T){
        cond.push_back(calcConductivity(i));
    }
    return cond;
}

double calcConductivity(double T) {
    return 7.288e-8 * pow(T, 3.0) + 7.1e-6 * pow(T,2.0) - 6.663e-4 * T + 0.14064;
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
    return (120.22 * exp(-0.03056*T))/100;
}

double TempToEnthalpy(double T) {
    return calcSpHeat(T)*T;
}
