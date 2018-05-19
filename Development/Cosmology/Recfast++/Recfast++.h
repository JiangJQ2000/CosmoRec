//====================================================================================================================
// Author: Jens Chluba & Luke Hart
// First implementation: Oct 2010
// CITA, University of Toronto
// All rights reserved.
//====================================================================================================================
// 23.04.2018: added structures to controll all the parameters and data passed on
// 21.02.2018: added cosmic time setup for internal use and decaying particles
// 10.01.2018: compute better high-z Te solution (JC)
// 03.04.2017: added switches using Saha and changed redshift grid (JC)
// 04.01.2017: added new routines for separate H and HeI correction functions (LH)
// 05.12.2016: added dependence of Saha_HeIII on fundamental constants

//====================================================================================================================
//
// Main Recfast++ header with all required function to access recombination problem as standalone code or from
// CosmoRec. the solution to the RECFAST ODE-system is computed. Different driver routines are provided.
//
//--------------------------------------------------------------------------------------------------------------------
// important note: all free particle fractions are relative to the number oh hydrogen
// nuclei. I.e. xe=ne/nH, xp=np/nH etc.
//====================================================================================================================

#ifndef RECFAST_PP_H
#define RECFAST_PP_H

#define RECFASTPP_VERSION "Recfast++ v1.2"

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>

#include "constants.Recfast.h"

using namespace std;

//====================================================================================================================
struct CosmoPars
{
    double Yp;
    double T0;               // Temperature of CMB at z=0
    double Omega_m;          // Omega matter
    double Omega_b;          // Omega Baryons
    double Omega_k;          // Omega Curvature
    double Omega_L;          // Omega Lambda
    double h100;
    double Neff;             // effective number of neutrinos
    
    //========================================================================================
    // set default parameters
    //========================================================================================
    CosmoPars()
    {
        h100=0.71; T0=2.726; Yp=0.24, Neff=3.046;
        Omega_m=0.26, Omega_b=0.044, Omega_L=0.0; Omega_k=0.0;
    }
    
    //========================================================================================
    // show variables
    //========================================================================================
    void show()
    {
        cout << " Yp= " << Yp << " T0= " << T0 << endl;
        cout << " Omega_m= " << Omega_m   << " Omega_b= " << Omega_b
             << " Omega_L= " << Omega_L   << " Omega_k= " << Omega_k
             << " Neff= " << Neff << endl;
        cout << " Hubble constant in units of 100 km s^-1 Mpc^-1: " << h100 << endl;
    }
};

//====================================================================================================================
struct RecPars
{
    int verbosity;
    int npts;
    double zstart, zend;
    
    double A2s1s;              // Hydrogen 2gamma rate
    double F;                  // Recfast fudge factor
    bool include_CF;           // correction function switch
    double f_ann;              // annihilation cases
    double f_dec, Gamma_dec;   // decay cases
    
    double B0, nB;             // primordial magnetic fields
    bool turb, ambi;
    
    double aS, mS, pS;         // alp/alp_r, me/me_r, power for (1+z)^p
    int RM;                    // mode
    
    //========================================================================================
    // parameters that are set internally (depending on runmode)
    //========================================================================================
    bool do_separate_H_He;     // boolean parameter for correction function settings
                               // 0 for the f_h+f_he [default]
                               // 1 for the new individual correction functions
                               // for variation of fundamental constants should set ==1

    bool eval_ion_Tr;          // boolean parameter for temperature evaluation of beta_phot
                               // 0: Photo-ionization with T=Te (Recfast default, incorrect)
                               // 1: T=Tg (phys.correct)
                               // for large heating (e.g., PMF and Decaying particles) set 1

    //========================================================================================
    // set default parameters
    //========================================================================================
    RecPars()
    {
        verbosity=0;
        npts=10000; zstart=2.5e+4; zend=0.001;
        A2s1s=8.22458;         // Recfast default value
        F=1.14;
        include_CF=1;
        f_ann=0.0;
        f_dec=0.0; Gamma_dec=0.0;
        B0=nB=0.0; turb=ambi=0;
        aS=mS=1.0; pS=0.0; RM=0;
        //-------------------------------------//
        do_separate_H_He=eval_ion_Tr=0;
    }
    
    //========================================================================================
    // show variables
    //========================================================================================
    void show()
    {
        cout << " Recombination verbosity = " << verbosity << endl;
        cout << " zstart = " << zstart << " zend = " << zend << " npts= " << npts << endl;
        cout << " Hydrogen 2gamma rate = " << A2s1s << endl;
        cout << " Recfast fudge F = " << F << endl;
        cout << " Include correction function: " << include_CF << endl;
        cout << " f_ann = " << f_ann << endl;
        cout << " f_dec = " << f_dec << " Gamma_dec = " << Gamma_dec << endl;
        cout << " B0 = " << B0 << " nB = " << nB << endl;
        cout << " turb = " << turb << " ambi = " << ambi << endl;
        cout << " alp/alp_r = " << aS << " alp/alp_r = " << mS << " power for (1+z)^p = " << pS << endl;
        cout << " mode for variation of constants: " << RM << endl;
    }
};

//====================================================================================================================
// compute ionization history like in Recfast-module
//====================================================================================================================
int Xe_frac(const CosmoPars &Cosmo, const RecPars &Rec, vector<double> &zarr,
            vector<double> &Xe_H, vector<double> &Xe_He, vector<double> &Xe,
            vector<double> &TM);

int Xe_frac(const CosmoPars &Cosmo, const RecPars &Rec, vector<double> &zarr,
            vector<double> &Xe_H, vector<double> &Xe_He, vector<double> &Xe,
            vector<double> &dXe, vector<double> &TM);

int Xe_frac(const CosmoPars &Cosmo, const RecPars &Rec, vector<double> &zarr,
            vector<double> &Xe_H, vector<double> &Xe_He, vector<double> &Xe,
            vector<double> &dXe, vector<double> &dX_H, vector<double> &TM);

//====================================================================================================================
// to start computation with a given initial solution at low z
// zarr contains the points at which the solution should be obtained
// rescaling of f_Xe from initial condition is used; For this dXei is needed
//====================================================================================================================
int Xe_frac_rescaled(const CosmoPars &Cosmo, const RecPars &Rec, const vector<double> &zarr,
                     vector<double> &Xe_H, vector<double> &Xe_He, vector<double> &Xe,
                     vector<double> &TM,
                     double Xe_Hi, double Xe_Hei, double Xei, double TMi, double dXei);

#endif
//====================================================================================================================
//====================================================================================================================
