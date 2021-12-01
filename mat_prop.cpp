//
// Created by thinkcomputational on 25/11/21.
//

#include "mat_prop.h"
#include <cmath>
#include <vector>
double calcTemp(double h) {
    return 7.81408e-8 * pow(h,3.0) - 2.2223e-4 * pow(h,2.0) + 0.55 * h - 17.0745;
}

std::vector<double> EnthalpyToTemp(std::vector<double> h) {
    std::vector<double>T;
    for (auto i:h){
        T.push_back(calcTemp(i/1000.0));
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
    std::vector<double>C = {5.99e-16, -8.5072e-13,
                            5.04056e-10, -1.61879e-7,
                            3.06161e-5, -3.47386e-3,
                            0.231286, -8.45598, 141.571};
    double sum = 0;

    for (int i=0; i<C.size(); ++i){
        sum = sum + C[8-i]* pow(T, static_cast<double>(i));
    }
    return sum*0.001;
}

double TempToEnthalpy(double T) {
    double h = -5.3268e-8 * pow(h, 3.0) + 1.7913e-3 * pow(h, 2.0) + 1.83415 * h + 32.06383;
    return h*1000;
}
