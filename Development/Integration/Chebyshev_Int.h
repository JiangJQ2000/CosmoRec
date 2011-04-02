//-----------------------------
// Author Jens Chluba Jan 2010
//-----------------------------
#ifndef CHEBYSHEV_INT_H
#define CHEBYSHEV_INT_H

#include <iostream>
#include <string>

using namespace std;

void compute_nested_Chebyshev_table(int pmax);

double Chebyshev_test_fcn(double x);
double Chebyshev_test_fcn_log(double lgx);

//====================================================================================
// Driver that loops over the different Chebyshev rules
//====================================================================================
int Integrate_using_Chebyshev(double a, double x0, double b, double epsrel, double epsabs, double (*fptr)(double), int *neval, double *r, string mess="");

#endif
