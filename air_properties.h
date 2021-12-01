//
// Created by thinkcomputational on 29/11/21.
//

#ifndef SOLAR_ANALYSIS_AIR_PROPERTIES_H
#define SOLAR_ANALYSIS_AIR_PROPERTIES_H
namespace airProp {
    double mu(double &Ti) {
        double T = Ti + 273.15;
        double mu = 2.5914 * pow(10.0, -15.0) * pow(T, 3.0) - 1.4346 * pow(10.0, -11.0) * pow(T, 2.0) + 5.0523 * pow(10.0, -8.0) * T +
                    4.1130 * pow(10.0, -6.0);
        return (mu);
    }

    double k(double &Ti) {
        double T = Ti + 273.15;
        double k =
                1.5797 * pow(10.0, -17.0) * pow(T, 5.0) + 9.46 * pow(10.0, -14.0) * pow(T, 4.0) + 2.2012 * pow(10.0, -10.0) * pow(T, 3.0) -
                2.3758 * pow(10.0, -7.0) * pow(T, 2.0) + 1.7082 * pow(10.0, -4.0) * pow(T, 1.0) - 7.488 * pow(10.0, -3.0);
        return (k);
    }

    double cp(double &Ti) {
        double T = Ti + 273.15;
        double cp = 1.3864 * pow(10.0, -13.0) * pow(T, 4.0) - 6.4747 * pow(10.0, -10.0) * pow(T, 3.0) +
                    1.0234 * pow(10.0, -6.0) * pow(T, 2.0) - 4.3282 * pow(10.0, -4.0) * pow(T, 1.0) + 1.0613;
        return (cp * 1000.0);
    }

    double rho(double &Ti) {
        double T = Ti + 273.15;
        double rho = 345.57 / (T - 2.6884);
        return (rho);
    }
}
#endif //SOLAR_ANALYSIS_AIR_PROPERTIES_H
