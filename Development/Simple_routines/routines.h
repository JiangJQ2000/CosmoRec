//======================================================================================
// Author: Jens Chluba 
// last modified: Oct 2010
//
// Purpose: collection of several simple routines
//======================================================================================

#ifndef ROUTINES_H
#define ROUTINES_H

#include <iostream>
#include <string>

using namespace std;

//======================================================================================
// Special functions
//======================================================================================
double Dawson_Int(const double &x); // Dawson integral; Based on Antia
double erf_JC(const double &x);     // Error function for a real argument
double erfc_JC(const double &x);    // complementary Error function for a real argument

double max(double a, double b);
double min(double a, double b);

double factorial(int n);
double log10factorial(int n);
double log10factorial_full(int n);
double factorial_corrfac(int n);

//======================================================================================
// checking for nan
//======================================================================================
bool isnan_JC(double a);

//======================================================================================
// routines for interpolation; based in GSL
//======================================================================================
int calc_spline_coeffies_JC(int nxi, const double *za, const double *ya, string variable="");
void update_spline_coeffies_JC(int memindex, int nxi, const double *za, const double *ya, string variable="");
double calc_spline_JC(double x, int memindex, string mess="");

void free_spline_JC(int memindex, string mess="");
void free_all_splines_JC();
void show_spline_memory();

//======================================================================================
// npol-1 is degree of the interpolating polynomial
//======================================================================================
void polint_JC(const double *xa, const double *ya, int na, const double x, int npol, double *y, double *dy);
void polint_JC(const double *xa, const double *ya, int na, const double x, int &istart, int npol, double *y, double *dy);

//======================================================================================
// simple grid setup
//======================================================================================
void init_xarr(double x0, double xm, double *xarr, int npts, int method_flag, int mess_flg);
void init_xarr(double x0, double xm, double *xarr, int npts, int method_flag);

//======================================================================================
void wait_f_r();
void wait_f_r(string mess);
void wait_f_r(int num);
void wait_f_r(double num);

//======================================================================================
void locate_JC(const double xx[], unsigned long n, double x, unsigned long *j);
void hunt(const double xx[], unsigned long n, double x, unsigned long *jlo);

//======================================================================================
// i/o modules
//======================================================================================
string int_to_string(int i);
string int_to_string(int i, int ni);

//======================================================================================
// root-finding methods
//======================================================================================
double find_root(double (* func)(double *), double x1, double x2, double xacc);
double find_root_brent(double (* func)(double *), double x1, double x2, double xacc);
double find_root_brent(double (* func)(double *, void *p), void *p, double x1, double x2, double xacc);

#endif
