//
// Created by thinkcomputational on 25/11/21.
//

#ifndef SOLAR_ANALYSIS_TDMA_H
#define SOLAR_ANALYSIS_TDMA_H
#include <vector>
/* here the equation is such that c*x(i-1) + a*x(i) + b*x(i+1) = d */
// Create the temporary vectors
// Note that this is inefficient as it is possible to call
// this function many times. A better implementation would
// pass these temporary matrices by non-const reference to
// save excess allocation and deallocation
void solve(std::vector<double> c,
           std::vector<double> a,
           std::vector<double> b,
           std::vector<double> d,
           std::vector<double>& f);
#endif //SOLAR_ANALYSIS_TDMA_H
