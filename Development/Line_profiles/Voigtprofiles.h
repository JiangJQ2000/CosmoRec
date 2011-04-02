//==========================================================================================
// Author: Jens Chluba 
// last modified: Oct 2010
//
// Purpose: computation of the Voigt-profile and integrals over it
//==========================================================================================
// note: x  = (nu-nu0)/ DnuT
//       a == Voigt parameter
//==========================================================================================

#ifndef VOIGTPROFILE_NEW_H
#define VOIGTPROFILE_NEW_H

#include "physical_consts.h"

class Voigtprofile_Base
{
private:
    
    // --------------- Atomic data ------------
    double nu21;                                   // Hz
    double lambda21;                               // cm
    double A21;                                    // 1/s
    double f12;                                    // 1/s  
    double Gamma;                                  // total line width Gamma=Sum A_ij
    double AM;                                     // atomic mass in units of mp
    double k_mHc2;
    
public:
    
    //======================================================================================
    // Konstructors and Destructors
    //======================================================================================
    Voigtprofile_Base();
    ~Voigtprofile_Base();
    
    void Set_atomic_data(double nu21, double lambda21, double A21, double f12, double AM=1);
    void Set_atomic_data_Gamma(double nu21, double lambda21, double A21, double f12, double Gamma, double AM=1);
    double Get_nu21() const { return nu21; } 
    double Get_A21() const { return A21; } 
    double Get_f12() const { return f12; } 
    double Get_lambda21() const { return lambda21; } 
    double Get_Gamma() const { return Gamma; } 
    double Get_AM() const { return AM; } 
    
    void display_atomic_data() const;
    
    double DnuT(double Tm){ return nu21*sqrt( 2.0*k_mHc2*Tm/AM );} // Doppler width
    double DnuT_nu12(double Tm){ return sqrt( 2.0*k_mHc2*Tm/AM );} // Doppler width 
    double aVoigt(double Tm){ return Gamma/FOURPI/DnuT(Tm); }
    double a2Dn_n0(double a){ return Gamma/FOURPI/a/nu21; }
    double nu2x(double nu, double Tm){ return (nu-nu21)/DnuT(Tm);}
    double x2nu(double xD, double Tm){ return xD*DnuT(Tm)+nu21;}
    
    //======================================================================================
    // simple approximations of Voigt-profile
    //======================================================================================
    double phi_DC(double x);                        // normal Doppler core
    double phi_DC(double x, double a);              // renormalized Doppler core
    double phi_wings(double x, double a);           // profile in the wings

    //======================================================================================
    // simple approximations for the 1. derivative of Voigt-profile
    //======================================================================================
    double dphi_dx_DC(double x);                    // normal Doppler core
    double dphi_dx_DC(double x, double a);          // renormalized Doppler core
    double dphi_dx_wings(double x, double a);       // profile in the wings
    
    double xi_Int_wings(double x, double a);        // Integral over wing region Int_-infinity^x phi_wings dx
};

class Voigtprofile_Dawson : public Voigtprofile_Base
{
private:
    
    double F, H[9], dH[9], IH[9];
    double xi_Int_Dawson(double xmin, double xmax, double a);   
    
public:
    
    //======================================================================================
    // Konstructors and Destructors
    //======================================================================================
    Voigtprofile_Dawson();
    ~Voigtprofile_Dawson();
    
    double Dawson(double x);             // Dawson-integral exp(-x^2) * int_0^x  exp(y^2) dy
    double H0(double x){ return exp(-x*x); }
    double H1(double x){ return 2.0/SQRT_PI*(2.0*x*Dawson(x)-1.0); }
    double H2(double x){ return exp(-x*x)*(1.0-2.0*x*x); }
    double H3(double x){ return 4.0/3.0/SQRT_PI*( x*x-1.0 + Dawson(x)*x*(3.0-2.0*x*x) ); }
    double H4(double x){ return exp(-x*x)*(3.0+4.0*x*x*(x*x-3.0) )/6.0; }
    
    void eval_Hn(double x);
    void display_Hn();
    void eval_dHn(double x);
    void display_dHn();
    void eval_IHn(double x1, double x2);
    void display_IHn();
    
    //======================================================================================
    // different versions for the Voigt-profile
    //======================================================================================
    double phi(double x, double a);                          // Voigt-profile based on expansion in a (Mihalas) 
    double dphi_dx(double x, double a);                      // 1. derivative of Voigt-profile  
    
    //======================================================================================
    // integrals over the Voigt-function
    //======================================================================================
    double xi_Int(double xmin, double xmax, double a);       // int_{-xmin}^xmax phi(y) dy
    double xi_Int(double x, double a);                       // int_{-infinity}^x phi(y) dy
    double xi_Int_xc(double x, double a);                    // int_{-xc}^x phi(y) dy with xc=-nu0/DnuD
    double chi_Int(double x, double a);                      // int_x^{infinity} phi(y) dy
};

#endif
