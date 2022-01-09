//
// Created by thinkcomputational on 26/11/21.
//
#include "HTC.h"
#include "air_properties.h"
#include <cmath>
double NuG(double Re, double Pr1, double Pr2) {
    double fby8 = std::pow(1.82*std::log10(Re) - 1.64, -2.0)/8.0;
    double num = fby8*(Re - 1000.0)*Pr1;
    double den = 1 + 12.7*std::pow(fby8, 0.5)*(std::pow(Pr1,2.0/3.0) - 1.0);
    return (num/den)* std::pow(Pr1/Pr2, 0.11);
}
double NuH(double Re, double Pr, double DbyL) {
    double num = 0.0668*DbyL*Re*Pr;
    double den = 1.0 + 0.04*std::pow(DbyL*Re*Pr, 0.67);
    return 3.66 + num/den;
}

double skyT(double Tamb) { //tamb in degree celcius
    double eps, phi, Tdp;
    phi = 17.27*Tamb/(Tamb + 273.3);
    Tdp = 273.3*(log(phi) + phi)/(17.27 - log(phi) - phi);
    eps = 0.742 + 0.0062*Tdp;
    return Tamb*pow(eps,0.25);
}

double h_rad(double T, double Tamb) {
    double Tsky;
    Tsky = skyT(Tamb) + 273.15;
    double Ta = T + 273.15;
    return eps_ab*StefanConst*(Ta*Ta + Tsky*Tsky)*(Ta + Tsky);
}

double ConvNoWind(double Ra, double Pr) {
    double num = 0.387* pow(Ra, 1.0/6.0);
    double den = pow(1 + pow(0.559/Pr, 9.0/16.0), 8.0/27.0);
    return pow(0.6 + num/den, 2.0);
}

double ConvWind(double Re, double Pr1, double Pr2){
    double c, m, n;
    if (Re<40){
        c = 0.75;
        m = 0.4;
    } else if (Re<1000){
        c = 0.51;
        m = 0.5;
    }else if (Re < 2e5){
        c = 0.26;
        m = 0.6;
    } else{
        c = 0.076;
        m = 0.7;
    }
    return c * pow(Re, m) * pow(Pr1, 0.365) * pow(Pr1/Pr2, 1.0/4.0);
}

void updateHTCF(std::vector<double> Tf, std::vector<double> Ta,
                std::vector<double> &htcf) {
    double k = 2.0;
    std::vector<double> Tave(Tf.size());
    std::transform(Tf.begin(),Tf.end(),Ta.begin(),Tave.begin(), std::plus<double>());
    std::transform(Tave.begin(), Tave.end(), Tave.begin(),[k](double &c){return c/k;});
    double A = PI*Di*Di/4.0;
    double Re, Prf, Prs, vis;
    std::vector<double>K,Cp;
    K = theConductivity(Tave);
    Cp = spHeat(Tave);
    for (int i = 0; i < Tave.size(); ++i) {
        vis = mu(Tave[i]);
        Re = mdot*Di/A/vis;
        Prf = vis*Cp[i]/K[i];
        if (Re<=2400){
//            double DbyL = Di/(dx*(1+i)*dx);
//            if (1/DbyL/Re/Prf < 0.05){
//                htcf[i] = NuH(Re,Prf, DbyL)* calcConductivity(Tf[i])/Di;
//            } else
                htcf[i] = 3.66* calcConductivity(Tf[i])/Di;
        } else{
            Prs = mu(Ta[i])* calcSpHeat(Ta[i])/ calcConductivity(Ta[i]);
            htcf[i] = NuG(Re, Prf, Prs)*calcConductivity(Tf[i])/Di;
        }
    }
}

void updateHTCA(std::vector<double> Ta,
                double w, double Tamb,
                std::vector<double> &htca){
    double k = 2.0;
    double A = PI*Di*Di/4.0;
    double Re, Prf, Prs, vis, Ra, alpha;
    unsigned int N = Ta.size();
    std::vector<double> Tave(Ta.size()), beta(Ta.size());
    std::transform(Ta.begin(), Ta.end(), Tave.begin(), std::bind2nd(std::plus<double>(), Tamb));
    std::transform(Tave.begin(), Tave.end(), Tave.begin(),[k](double &c){return c/k;});
    std::transform(Tave.begin(), Tave.end(), beta.begin(), [](double &c){return 1.0/(c+273.15);});
    double K, Cp, rho;
    if (w==0){
        for (int i = 0; i < N; ++i) {
            vis = airProp::mu(Tave[i]);
            Cp = airProp::cp(Tave[i]);
            rho = airProp::rho(Tave[i]);
            K = airProp::k(Tave[i]);
            Prf = vis*Cp/K;
            alpha = K/Cp/rho;
            Ra = g*beta[i]*std::fabs(Tave[i]-Tamb)*pow(Do, 3.0)*rho/alpha/vis;
            htca[i] = ConvNoWind(Ra, Prf) * K/Do;
        }
    } else{
        for (int i = 0; i < N; ++i) {
            vis = mu(Tave[i]);
            Re = mdot*Di/A/vis;
            Prf = vis*Cp/K;
            Prs = airProp::mu(Ta[i])*airProp::cp(Ta[i])/airProp::k(Ta[i]);
            htca[i] = ConvWind(Re, Prf, Prs) * K/Do;
        }
    }
}

void updateHTCR(std::vector<double> Ta, double Tamb,
                std::vector<double> &htcr) {
    unsigned int N = Ta.size();
    for (int i = 0; i < N; ++i) {
        htcr[i] = h_rad(Ta[i], Tamb);
    }
}

double calculateHloss(double dT, double L, double &mdot) {
    //L in meter, dT in temperature difference
    return (0.156*dT)*L*mdot;
}

double calcQloss(std::vector<double> &T, double Tr, std::vector<double> &htc) {
    double sum = 0.0;
    for (int i = 0; i < T.size(); ++i) {
        sum += htc[i]*(T[i] - Tr);
    }
    return sum;
}

