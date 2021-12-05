//
// Created by thinkcomputational on 04/12/21.
//

#include "Solar_data.h"
#include <math.h>
#define Isc = 1367;
std::vector<double> HourAngle(double t_start, double dt,
                              unsigned int N) {
    std::vector<double>hr_angle;
    double start = (12.0-t_start)*15.0;
    double inter = 15.0*dt/60;
    for (int i = 0; i < N; ++i) {
        hr_angle.emplace_back((i+1)*inter);
    }
    return hr_angle;
}

double declination(int day_no) {
    double delta = sin((360.0/365.0)*(static_cast<double>(day_no) + 284.0)*3.1415/180);
    return 23.45*delta;
}

double E0(double day_no) {
    double g = (2*3.1416*(day_no-1)/365)*3.14156/180;
    double Eo = 1.000110 + 0.034221*cos(g) + 0.001280*sin(2.0*g) + 0.000719*cos(2*g) + 0.000077*sin(2.0*g);
    return Eo;
}

double cosTheta(double delta, double phi, double omega) {
    double d, p, o;
    d = delta * 3.1415/180;
    p = phi * 3.1415/180;
    o = omega * 3.1415/180;
    return sin(d)* sin(p) + cos(d)* cos(p)* cos(o);
}

void GHItoDNI(int D, double lat, double t_start, double dt,
              std::vector<double> &q) {
    unsigned int n = q.size();
    std::vector<double> GHI = q,
                        I0(n), omega(n);
    omega = HourAngle(t_start, dt, n);
    double delta = declination(D);

}


