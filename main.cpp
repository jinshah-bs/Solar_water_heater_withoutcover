#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include "read_data.h"
#include "TDMA.h"
#include "HTC.h"
#include "Write_data.h"
#include "Solar_data.h"
double eps_ab = 0.9,
       alpha_ab = 0.9,
       StefanConst = 5.67e-8;
double l = 2.5,
        Di = 0.02355,
        Do = 0.02675,
        CR = 33;
int n = 50;
double mdot, dx, lpm;
std::vector <std::vector<double>> MetData;
double t = 0, dt = 1.0;
double Ac, Aa, Pi, Po, Vi;
//Details of tank
double Dt = 0.3,
       Lt = 0.707;
int nt = 10;
std::vector<std::vector<double>> outabs, outFluid, outtank, outdata;
std::vector<double>outtime;
int main()
{
    //input whether data
    MetData = parse2DCsvFile("Whether_data.csv");
    GHItoDNI(274, 8.8932, 6.0, 5.0, MetData);
    //[i][0]= Tamb, [i][1]= wind, [i][2]=Q
    // Details of grid
    lpm = 10; //liter per minute
    dx = l/n;
    mdot = lpm * calcDencity(MetData[0][0])/60/1000;
    Ac = PI*Di*Di/4.0;
    Aa = PI*(Do*Do - Di*Di)/4.0;
    Pi = PI*Di;
    Po = PI*Do;
    Vi = Ac*dx;
    double dxT = Lt/nt,
           At = PI*Dt*Dt/4.0,
           Vt = At*dxT;
    std::vector<double>Tf(n,MetData[0][0]),
                       Ta(n,MetData[0][0]),
                       Tt(nt, MetData[0][0]);
    std::vector<double>h(n, TempToEnthalpy(MetData[0][0])),
                       ht(nt, h[0]); //enthalpy initial value 30 degree
    std::vector<double>rho = dencity(Tf),
                       rhoT = dencity(Tt),
                       K = theConductivity(Tf),
                       KT = theConductivity(Tt),
                       Cp = spHeat(Tf),
                       CpT = spHeat(Tt);
    // thermophysical properties of copper
    double rhoAB = 8960,
           kAB = 386,
           cpAB = 3850;

    /* here the equation is such that c*x(i-1) + a*x(i) + b*x(i+1) = d */
    std::vector<double> a(n,2.0*kAB*dt/dx/dx + rhoAB*cpAB), b(n,-kAB*dt/dx/dx),
                        c(n, -kAB*dt/dx/dx), d(n, 0.0);
    c[0] = 0.0; b[n-1] = 0.0, a[0] = a[n-1] = a[0] - kAB*dt/dx/dx; //bundary conditoin
    std::vector<double> htcf(n, 0.0), htca(n, 0.0), htcr(n,0.0);
    std::vector<double> Tfo = Tf,
                        Tfi = Tf,
                        Tao = Ta,
                        Tai(n, 100.0),
                        ho = h,
                        hto = ht;
    // defining time variable
    int dwrite = 300/static_cast<int>(dt);
    int nouter = 1,
        ninner = 1,
        ntheta = 1,
        nloop = 0;
    double q, w, Tamb, Ki, theta, h_mul, q_loss_tank;
    std::vector<double> DT(n,0.0);
    outFluid.emplace_back(Tf);
    outabs.emplace_back(Ta);
    outtank.emplace_back(Tt);
    outtime.emplace_back(0.0);
    outdata.emplace_back(std::vector<double>{10, 12, 13, 14});
    while (nouter<140){
        Ki = -2.2307e-4 *theta - 1.1e-4 * pow(theta, 2.0)  + 3.18596e-6 * pow(theta, 3.0) - 4.85509e-8 * pow(theta, 4.0);
        Ki = std::max(0.0, std::min(1.0, 1.0 + Ki));
        q = InterpolateWhetherData(nouter, 2, ninner)*0.82340*Ki*alpha_ab*CR;
        w = InterpolateWhetherData(nouter, 1, ninner);
        Tamb = InterpolateWhetherData(nouter, 0, ninner);
        theta = 15.0*PI*dt*ntheta/180.0/3600.0;
        ho = h; Tfo = Tf; Tao = Ta;
        updateThermalProp(rho, K, Cp, Tf);
        double error_T = 10.0;
        std::vector<double> error(n,0.0);
        // Inner loop
        h[0] = ht[nt-1] - calculateHloss(std::fabs(Tamb-Tt[nt-1]), 7.0, mdot); // Change this value
        // Update the diffusion term
        for (int i = 1; i < n-1; ++i) {
            DT[i] = (K[i+1]*Tf[i+1] - (K[i+1] + K[i])*Tf[i] + K[i-1]*Tf[i-1])/dx/dx;
        }
        DT[0] = (K[1]*Tf[1] - (K[1] + K[0])*Tf[0] + K[0]*Tf[0])/dx/dx; // Change Tf[0] with outlet temperature
        DT[n-1] = (-K[n-1]*Tf[n-1] + K[n-2]*Tf[n-2])/dx/dx; // Change Tf[0] with outlet temperature
        // loop Starts here
        updateHTCF(Tf, Ta, htcf);
        updateHTCA(Ta, w, Tamb, htca);
        updateHTCR(Ta, Tamb, htcr);
        nloop = 1;
        while (error_T > 1.0){
            // Solve for absorber
            for (int i = 0; i < n; ++i) {
                d[i] = Po/Aa *(q - htca[i]*(Tao[i] - Tamb) - htcr[i]*(Tao[i]- skyT(Tamb)))*dt - Pi*htcf[i]*(Tao[i]-Tfo[i])*dt/Aa + rhoAB*cpAB*Tao[i];
            }
            // Boundary condition
            solve(c, a, b, d, Tai);
            std::transform(Ta.begin(), Ta.end(), Tai.begin(), error.begin(), [&](double l, double r)
            {
                return std::abs(l - r);
            });
            error_T = *std::max_element(error.begin(), error.end());
            Ta = Tai;
            // Solve fluid domain
            for (int i = 1; i < n; ++i) {
                h_mul = 1.0/(rho[i] + mdot*dt/Vi);
                h[i] = (rho[i]*ho[i] + (mdot*dt/Vi)*h[i-1] + DT[i]*dt + htcf[i]*(Tao[i] - Tfo[i])*Pi*dt/Ac)*h_mul;
            }
            Tfi = EnthalpyToTemp(h);
            std::transform(Tf.begin(), Tf.end(), Tfi.begin(), error.begin(), [&](double l, double r)
            {
                return std::abs(l - r);
            });
            error_T = std::max(*std::max_element(error.begin(), error.end()), error_T);
            std::cout << "Error in T is " << error_T << "at iteration" << nloop << "and n_inner" << ninner <<std::endl;
            Tf = Tfi;
            nloop++;
        }

        //calculate for Tank.
        hto = ht;
        ht[0] = h[n-1] - calculateHloss(std::fabs(Ta[n-1]-Tamb), 3.0, mdot);
        updateThermalProp(rhoT, KT, CpT, Tt);
        for (int i = 1; i < nt-1; ++i) {
            DT[i] = (KT[i+1]*Tt[i+1] - (KT[i+1] + KT[i])*Tt[i] + KT[i-1]*Tt[i-1])/dxT/dxT;
        }
        DT[0] = (KT[1]*Tt[1] - (KT[1] + KT[0])*Tt[0] + KT[0]*Tt[0])/dxT/dxT; // Change Tf[0] with outlet temperature
        DT[nt-1] = (-KT[nt-1]*Tt[nt-1] + KT[nt-2]*Tt[nt-2])/dxT/dxT; // Change Tf[0] with outlet temperature
        for (int i = 1; i < nt; ++i) {
            h_mul = 1.0/(rhoT[i] + mdot*dt/Vt);
            q_loss_tank = 11.3048* pow(Tt[i] - Tamb, 0.2823)*mdot*dt/At;
            ht[i] = (rhoT[i] * hto[i] + (mdot*dt/Vt)*ht[i-1] + DT[i]*dt)*h_mul;
        }
        Tt = EnthalpyToTemp(ht);
        ninner += 1;
        ntheta += 1;
        if (ninner > 300.0/dt){
//            write(Ta, nouter*5, "T_absorber");
//            write(Tf, nouter*5, "T_fluid");
//            write(Tt, nouter*5, "T_tank");
            outFluid.emplace_back(Tf);
            outabs.emplace_back(Ta);
            outtank.emplace_back(Tt);
            outtime.emplace_back(nouter*5.0);
            ninner = 0;
            nouter += 1;
        }
        if (ntheta > static_cast<int>(3600.0/dt)){
            ntheta = 1;
        }
        std::cout << "max abs T " << *std::max_element(Tf.begin(), Tf.end()) << std::endl;
        std::cout << "outer iterator is " << nouter << std::endl;
        std::cout << "Efficiency of the system is " << mdot*(h[n-1] - h[0])/(q*Pi*Do*l/0.82340/Ki/alpha_ab)
                  << std::endl;
        std::cout << "-------------------------------------------------" << std::endl;
    }
    write(outFluid, outtime, "Fluid_temperature");
    write(outabs, outtime, "Absorber_temperature");
    write(outtank, outtime, "Tank_temperature");
    write(outdata, outtime, "Data");

}