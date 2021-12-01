//
// Created by thinkcomputational on 25/11/21.
//

#include "TDMA.h"

void solve(std::vector<double> c,
           std::vector<double> a,
           std::vector<double> b,
           std::vector<double> d,
           std::vector<double>& f) {

    int N = d.size();

    std::vector<double> P(N, 0.0);
    std::vector<double> Q(N, 0.0);

    // This updates the coefficients in the first row
    // Note that we should be checking for division by zero here
    P[0] = -b[0]/a[0];
    Q[0] = d[0]/a[0];
    // Create the c_star and d_star coefficients in the forward sweep
    for (int i=1; i < N; i++)
    {
        double m = 1.0 / (a[i] + c[i] * P[i-1]);
        P[i] = -b[i] * m;
        Q[i] = (d[i] - c[i] * Q[i-1]) * m;
    }
    f[N-1] = Q[N-1];
    // This is the reverse sweep, used to update the solution vector f
    for (int j=N-1; j >= 0; j-- )
    {
        f[j] = P[j] * f[j+1] + Q[j];
    }
}