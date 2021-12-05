//
// Created by thinkcomputational on 04/12/21.
//

#ifndef SOLAR_ANALYSIS_SOLAR_DATA_H
#define SOLAR_ANALYSIS_SOLAR_DATA_H
#include <vector>
double E0(int day_no);
std::vector<double>HourAngle(double t_start, double dt, unsigned int N);
double declination(int day_no);
double cosTheta(double delta, double phi, double omega);
void GHItoDNI(int D, double lat, double t_start, double dt, std::vector<double> &q);
#endif //SOLAR_ANALYSIS_SOLAR_DATA_H
