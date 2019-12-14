//===========================================================================================
// Author: Jens Chluba
// First implementation: August 2018
//===========================================================================================

//===========================================================================================
// This simple module to add reionization
//===========================================================================================

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "Recfast++.h"
#include "cosmology.Recfast.h"
#include "evalode.Recfast.h"
#include "DM_annihilation_decay.Recfast.h"
#include "Reionization.Recfast.h"

using namespace std;
using namespace RECFAST_physical_constants; // defined in constants.Recfast.h
using namespace Recfastpp_Cosmology;        // defined in cosmology.Recfast.h

//===========================================================================================
// evaluate reionization terms following Poulin et al.
//===========================================================================================
double evaluate_dEdz_JC_old(const RecPars &Rec, double z, double Hz)
{
    double p=0.4;         // dE/dt ~ (1+z)^p
    double z_reion=6.0;    // redshift of reionization
    double f_reion=0.6;    // effective eV per HI and Hubble at z_reion

    double Hz_fac=H_z(z_reion)/Hz;
    double Fcal=pow((1.0+z)/(1.0+z_reion), p) * Hz_fac
                *exp( -0.12*(pow((1.0+z)/(1.0+z_reion), 4.0)-1.0));
    
    return -f_reion*Fcal;                             // '-' from trans. dt --> dz
}

//===========================================================================================
double evaluate_dEdz_JC(const RecPars &Rec, double z, double Hz)
{
    double p=1.5;         // dE/dt ~ (1+z)^p
    double z_reion=6.0;    // redshift of reionization
    double f_reion=15.0;    // effective eV per HI and Hubble at z_reion
    
    double Hz_fac=H_z(z_reion)/Hz;
    double Fcal=pow((1.0+z)/(1.0+z_reion), p) * Hz_fac
//                *exp( -0.12*pow((1.0+z)/(1.0+z_reion), 4.0) );
                 *exp( -pow((1.0+z)/(1.0+z_reion), 6.0) );

    return -f_reion*Fcal;                             // '-' from trans. dt --> dz
}

//===========================================================================================
double evaluate_dEdz_Poulin(const RecPars &Rec, double z, double Hz)
{
    double As=1.0, fesc=0.2;
    double xi_ion=1.38e+53*13.6; // eV sec^-1 (Msun/yr)^-1
    
    double a=0.01376;            // sec^-1 (Msun/yr) Mpc^-3
    double b=3.26, c=2.59, d=5.68;
    double rho_ion=a*pow(1.0+z, b)/(1.0+pow((1.0+z)/c, d));
    
    double dEdtdV=As*fesc*xi_ion*rho_ion*pow(1.0+z, 3);

    // Mpc^-3 --> m^-3
    dEdtdV/=pow(RF_Mpc, 3);
    // -Hz*(1+z) from dt --> dz and NH since everything is relative to Hydrogen number dens.
    return dEdtdV/( -Hz*(1.0+z)*NH(z) ) * ( 1.0+tanh(25.0-z) )/2.0;
}

//===========================================================================================
double evaluate_dEdz_Poulin_Xray(const RecPars &Rec, double z, double Hz)
{
    double fX=0.2;
    double A_X=3.4e+40 * 6.242e+11;  // eV sec^-1 (Msun/yr)^-1
    
    double a=0.01376;                   // sec^-1 (Msun/yr) Mpc^-3
    double b=3.26, c=2.59, d=5.68;
    double rho_ion=a*pow(1.0+z, b)/(1.0+pow((1.0+z)/c, d));
    
    double dEdtdV=A_X*fX*rho_ion*pow(1.0+z, 3);
    
    // Mpc^-3 --> m^-3
    dEdtdV/=pow(RF_Mpc, 3);
    // -Hz*(1+z) from dt --> dz and NH since everything is relative to Hydrogen number dens.
    return dEdtdV/( -Hz*(1.0+z)*NH(z) );
}

//===========================================================================================
// evaluation of ODE terms from reionization
//===========================================================================================
void evaluate_reionization_terms(double z, double Hz,
                                 double xp, double xHep,
                                 double CHe, double CH,
                                 double fvec[])
{
    RecPars Rec=Get_evalode_Vars();

    //=======================================================================================
    // using similar methods as for DM decay and annihilation
    //=======================================================================================
    double dE_dz=0.0;

    if(Rec.reion_model==1) dE_dz=evaluate_dEdz_Poulin(Rec, z, Hz);
    else if(Rec.reion_model==2) dE_dz=evaluate_dEdz_JC(Rec, z, Hz);

    //=======================================================================================
    // add corresponding terms
    //=======================================================================================
    double fvectemp[]={0.0, 0.0, 0.0};
    evaluate_ionization_heating_excitation_terms(z, Hz, xp, xHep, CHe, CH, dE_dz, fvectemp);
    
    if(Rec.reion_model==1)
    {
        if(z<200.0) fvec[0]=0.0;
        fvec[1]+=fvectemp[1];

        // evaluate X-ray heating
        dE_dz=evaluate_dEdz_Poulin_Xray(Rec, z, Hz);
        
        fvectemp[0]=fvectemp[1]=fvectemp[2]=0.0;
        evaluate_ionization_heating_excitation_terms(z, Hz, xp, xHep, CHe, CH, dE_dz, fvectemp);
        
        fvec[1]+=fvectemp[1];
        fvec[2]+=fvectemp[2];
    }
    else
    {
        fvec[0]+=fvectemp[0];
        fvec[1]+=fvectemp[1];
        fvec[2]+=fvectemp[2];
    }
    
    return;
}

//===========================================================================================
//===========================================================================================

