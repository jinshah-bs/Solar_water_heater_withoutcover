//
// Created by thinkcomputational on 26/11/21.
//

#ifndef SOLAR_ANALYSIS_HTC_H
#define SOLAR_ANALYSIS_HTC_H
#include <bits/stdc++.h>
#include "mat_prop.h"
#include <vector>
#define PI 3.14156
#define g 9.81
extern double eps_ab, StefanConst;
extern double l, Di, Do, dx, mdot;
double NuG(double Re, double Pr1, double Pr2);
double NuH(double Re, double Pr, double DbyL);
double skyT(double Tamb);
double h_rad(double T, double Tamb);
double ConvNoWind(double Ra, double Pr);
double ConvWind(double Re, double Pr1, double Pr2);
void updateHTCF(std::vector<double> Tf,
               std::vector<double> Ta,
               std::vector<double> &htcf);
void updateHTCA(std::vector<double> Ta,
                double w, double Tamb,
                std::vector<double> &htca);
void updateHTCR(std::vector<double> Ta, double Tamb, std::vector<double> &htcr);
double calculateHloss(double dT, double L, double &mdot);
double calcQloss(std::vector<double> &T, double Tr, std::vector<double> &htc);

#endif //SOLAR_ANALYSIS_HTC_H
