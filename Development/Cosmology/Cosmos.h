//========================================================================================
// Author: Jens Chluba (May 2003)
// CITA, University of Toronto
// All rights reserved.
//========================================================================================
// 19.05.2018: updated Recfast++ part with better parameter control [JC]
// 09.02.2018: initialization parameters now passed on as structure [JC]
// 10.01.2018: added DTe-function for high redshifts
// 04.12.2017: restructured Saha solution parts
// 16.05.2017: added sigT(z) to explicitly account for variations of fundamental constants
// 06.05.2017: added several useful cosmology functions [JC & Luke Hart]

//========================================================================================
// 22.07.2015
//========================================================================================
// fixed Saha populations for helium (not used before...)
//========================================================================================

//========================================================================================
// 20.10.2014
//========================================================================================
// y-parameter of photon and electron temperature
//========================================================================================

//========================================================================================
// 18.06.2012
//========================================================================================
// added possibility to set Hubble factor using external function pointer
//========================================================================================

//========================================================================================
// 08.06.2012
//========================================================================================
// added possibility to load Hubble factor from some external table
//========================================================================================

//========================================================================================
// 15.10.2011
//========================================================================================
// added simple routines that are needed to describe acoustic damping of modes
//========================================================================================

//========================================================================================
// 28.05.2008
//========================================================================================
// Set the variable fac_mHemH in "physical_consts.h" to take into
// account the fact that the helium mass is not 4*mH (Wong et al 2008)
// However, we only changed those variables that are important for the
// recombination computations.
//========================================================================================

#ifndef COSMOS_H
#define COSMOS_H

#include <iostream>
#include <cmath>
#include <string>
#include <vector>

#include "physical_consts.h"
#include "Recfast++.h"

using namespace std;

//========================================================================================
// first added 08.06.2012
//========================================================================================
class Hubble
{
    
private:
    
    int mem_index_Hz;
    double zmax, zmin;
    double (* ext_Hptr)(const double *);

    bool loaded;
    bool pointer;
    
public:

    Hubble();
    Hubble(const vector<double> &z, const vector<double> &Hz);
    Hubble(double (* Hptr)(const double *));
    ~Hubble();
    
    void init(const vector<double> &z, const vector<double> &Hz);
    void init(const double *z, const double *Hz, const int nz);
    void init(double (*Hptr)(const double *));
    void clear();
    
    bool check_limits(double z);
    
    double H(double z);
};


//========================================================================================
class Cosmos
{

private:
// --------------- Cosmology constants ------------

    double Nnu;                                 // number of neutrino species: 3 
    double g_rel;                               // 1.0+3.0*(7.0/8.0)*pow(4.0/11.0, 4.0/3.0)
    double zsRe;                                // start recombination

    //====================================================================================
    // initial data to be set 
    //====================================================================================
    double h100;                                // h=H0/100 km s^-1 Mpc^-1
    double T_CMB0;                              // Temperature of the CMB today in [K]
    double Y_p;                                 // Helium energy density 
    double O_T;                                 // total energy density 
    double O_mh2;                               // matter energy density 
    double O_bh2;                               // baryon energy density 
    double O_L;                                 // cosmological constant

    //====================================================================================
    // dependent consts
    //====================================================================================
    double h2, H0;                              // Hubble constant today 
    double rho_c, rho_g, rho_g_gr;              // gr cm^-3  
    double O_Th2, O_kh2, O_Lh2, O_m, O_b, O_k;
    double T27, Theta_CMB0, Theta27;
    double O_g, O_gh2;                          // photon density
    double O_rel, O_relh2;                      // relativistic particle density 
    double z_eq;                                // Matter radiation equality (Wu 1995)
    double Nb0;                                 // Baryon number density today

    //====================================================================================
    // auxilliary variables
    //====================================================================================
    void set_initials(double Nnuval);
    
    //====================================================================================
    // auxilliary functions
    //====================================================================================
    void calc_dependent_consts();

    //====================================================================================
    // setup for Xe & X1s etc
    //====================================================================================
    CosmoPars CosmoRF;
    RecPars RecRF;

    int recombination_history(int nzpts, double zstart, double zend,
                              vector<double> &zarr,
                              vector<double> &Xe_H, vector<double> &Xe_He,
                              vector<double> &Xe, vector<double> &TM);
    
    int recombination_history(int nzpts, double zstart, double zend,
                              vector<double> &zarr,
                              vector<double> &Xe_H, vector<double> &Xe_He,
                              vector<double> &Xe, vector<double> &dXe,
                              vector<double> &dX_H, vector<double> &TM);

    //====================================================================================
    // spline setup
    //====================================================================================
    int calc_coeff_X_spline(int M, double zstart, double zend);
    
    int n_Xe;
    int memindex_Xe;
    int memindex_dXe; 
    int memindex_XeH; 
    int memindex_dXeH; 
    int memindex_rho; 
    int memindex_XeHe; 
    bool splines_are_charged;
    bool spline_memory_allocated;

    void check_splines();
    double calc_spline(double z, int memindex);
    double calc_spline_log(double z, int memindex);
    //====================================================================================

    double RF_hPlanck, RF_kBoltz, RF_mElect; // local physical constants     
    
    //====================================================================================
    // to load an external hubble function from a table (08.06.2012)
    //====================================================================================
    Hubble loaded_Hz;
    
    //====================================================================================
    // to set models for varying consts [JC]
    //====================================================================================
    double ascale, mscale, pval;

public:

    //====================================================================================
    // structure to communicate parameters. Default values are set
    //====================================================================================
    struct parameters
    {
        //================================================================================
        // usual cosmological parameters
        //================================================================================
        double h100, T0, Yp, Nnu;
        double Omega_m, Omega_b; // Omega_m == Omega_cdm + Omega_b
        double Omega_L, Omega_k; // if Omega_L <= 0 set internally assuming Omega_tot=1
        
        //================================================================================
        // recombination specific
        //================================================================================
        double fudge;                 // old recfast recombination fudge factor
        double fDM_annihilation;      // DM annihilation efficiency
        bool CR_correction;           // switch to include CosmoRec correction functions
        
        double A2s1s;                 // Hydrogen 2gamma rate
        bool eval_ion_Tr;             // temperature evaluation of beta_phot
        bool do_separate_H_He;        // for correction function settings
        
        //================================================================================
        // default values set through constructor
        //================================================================================
        parameters()
        {
            h100=0.71; T0=2.726; Yp=0.24, Nnu=3.046;
            Omega_m=0.26, Omega_b=0.044, Omega_L=0.0; Omega_k=0.0;
            fudge=1.14; fDM_annihilation=0; CR_correction=1;
            A2s1s=8.22458; do_separate_H_He=eval_ion_Tr=0;
        }
    };

    //====================================================================================
    // Konstructors and Destructors
    //====================================================================================
    Cosmos();
    Cosmos(const parameters &pars);
    
    ~Cosmos();

    void init(const parameters &pars);
    
    void init_splines();
    
    //====================================================================================
    // to load an external hubble function from a table (08.06.2012)
    //====================================================================================
    void init_Hubble(const vector<double> &z, const vector<double> &Hz);
    void init_Hubble(const double *z, const double *Hz, const int nz);
    void init_Hubble(double (* Hptr)(const double *));
    void clear_Hubble(){ loaded_Hz.clear(); }
    
    //====================================================================================
    void display_cosmology();
    void dump_cosmology(string filename);

    double Get_Y_p() const { return Y_p;} 
    double Get_Nnu() const { return Nnu;}
    //
    double Get_O_bh2() const { return O_bh2;} 
    double Get_O_b() const   { return O_b;  } 
    //
    double Get_O_DMh2() const { return O_mh2-O_bh2;} 
    double Get_O_DM() const   { return O_m-O_b;    } 
    //
    double Get_O_mh2() const { return O_mh2;} 
    double Get_O_m() const   { return O_m;  } 
    //
    double Get_O_Lh2() const { return O_Lh2;  } 
    double Get_O_L() const   { return O_L;    } 
    // 
    double Get_O_g() const  { return O_g;   } 
    double Get_rho_c() const { return rho_c; } 
    double Get_rho_g() const { return rho_g; } 
    double Get_T27()   const { return T27;   } 
    double Get_T_CMB0() const{ return T_CMB0;} 

    //====================================================================================
    // CMB temperature
    //====================================================================================
    double TCMB(double z)    {return T_CMB0*(1.0+z);} 
    double Ttoz(double T_g)  {return T_g/T_CMB0-1.0;} 

    //====================================================================================
    // CMB temperature Theta = k Tg / me c^2
    //====================================================================================
    double ThetaCMB(double z){return Theta_CMB0*(1.0+z);} 
    double Thetatoz(double Theta_g){return Theta_g/Theta_CMB0-1.0;} 

    //====================================================================================
    // conversion Theta <--> Tg
    //====================================================================================
    double TtoTheta(double T_g)    {return const_kb_mec2*T_g;} 
    double ThetatoT(double Theta_g){return Theta_g/const_kb_mec2;} 

    double calc_Orel(double TCMB0, double Nnu, double h100)
    { 
        double H0=100.0*h100*(1.0e+3*1.0e+2)/const_Mpc;
        double a=PI*PI/15.0*const_kB*TCMB0*pow(const_kB*TCMB0/(const_hbar*const_cl), 3)
                      /pow(const_cl, 2); 
        double b=3.0*pow(H0, 2)/(8.0*PI*const_G);  
        return a/b*(1.0+Nnu*(7.0/8.0)*pow(4.0/11.0, 4.0/3.0));
    }
    
    //====================================================================================
    // compute Recfast system with a given initial solution at low z
    //====================================================================================
    int recombine_using_Recfast_system(int nzpts, double zi, double ze, double fDM, 
                                       double Xe_Hi, double Xe_Hei, double Xei, double TMi, 
                                       double dXei, const vector<double> &zarr,
                                       vector<double> &Xe_H, vector<double> &Xe_He,
                                       vector<double> &Xe, vector<double> &TM);

    //====================================================================================
    // Xe & X1s
    //====================================================================================
    double Xe_b(double z)              // this is Ne/Nb/[Ne(inf)/Nb(inf)] ~ Ne/Nb/(1-Yp/2)
    { 
        if(z>=zsRe) return 1.0; 
        // --- Xe (1-Y_p/2) =ne/nb ! but in recombination Xe=ne/nh 
        // --> factor of (1.0-Y_p)/(1.0-Y_p/2.0) 
        return (1.0-Y_p)/(1.0-Y_p/2.0*(2.0-1.0/fac_mHemH))*calc_spline_log(z, memindex_Xe);
    } 
    
    double Xe_Seager(double z)
    {                                              // xe as given by Seager 2000. Xe=ne/nh
        if(z>=zsRe) return (1.0-Y_p/2.0*(2.0-1.0/fac_mHemH))/(1.0-Y_p); 
        return calc_spline_log(z, memindex_Xe); 
    } 
    
    double dXe_dz_Seager(double z)
    {                                              // xe as given by Seager 2000. Xe=ne/nh
        if(z>=zsRe) return 0.0; 
        return calc_spline(z, memindex_dXe); 
    } 
    
    double dXe_dt_Seager(double z)
    {                                              // xe as given by Seager 2000. Xe=ne/nh
        if(z>=zsRe) return 0.0; 
        return -H(z)*(1.0+z)*calc_spline(z, memindex_dXe); 
    }
    
    //====================================================================================
    // X1s = NH1s/NH
    //====================================================================================
    double dX1s_dz(double z)
    { 
        if(z>=zsRe) return 0.0; 
        return -calc_spline(z, memindex_dXeH);
    } 
    
    double dX1s_dt(double z)
    { 
        if(z>=zsRe) return 0.0; 
        return +H(z)*(1.0+z)*calc_spline(z, memindex_dXeH);
    } 

    //====================================================================================
    // 30.06.2009 & Dec 2017 [JC]
    //====================================================================================
    double SahaBoltz_HeIII(double T);  // Xe^Saha during HeIII --> HeII recombination
    double SahaBoltz_HeII(double T);   // Xe^Saha during HeII  --> HeI  recombination

    double SahaBoltz_HeII1s(double T);
    double SahaBoltz_HeI1s(double T);
    double SahaBoltz_HI1s(double T);
    
    double XHeII1s(double z);
    double XHeI1s(double z);
    double X1s(double z);
    double XeSaha_corr(double z);     // correction to Saha problem

    //====================================================================================
    // (approximate) free proton fraction from Recfast
    //====================================================================================
    double Xe(double z){ return Xe_Seager(z); }                                   // Ne/NH
    double Xp(double z){ return 1.0-X1s(z);}
    
    //====================================================================================
    // Te_Tg from Recfast
    //====================================================================================
    double Te_Tg(double z){ return z>=zsRe ? 1.0/(1.0+1.0/kappa_cool(z)):calc_spline(z, memindex_rho); }
    double Te(double z){return TCMB(z)*Te_Tg(z);}
    double dTe_dz(double z);
    // added 10.01.2018 [JC]
    double Ntot(double z){ return NH(z)*(1.0+fHe()+Xe(z)); } // Ntot = Ne + Np + NHI + helium
    double kappa_cool(double z){ return 8.0/3.0*sigT(z)*const_cl/H(z)*rho_g_1_cm3(z)*Ne(z)/Ntot(z); }
    double DTe_high(double z){return -TCMB(z)/(1.0+kappa_cool(z));}

    //====================================================================================
    // different number densities
    //====================================================================================
    double Nb(double z){return Nb0*pow(1.0+z,3);}                                   // Nb==rho_b/mH
    double Ne(double z){return Xe_b(z)*(1.0-Y_p/2.0*(2.0-1.0/fac_mHemH))*Nb(z);}
    double N1s(double z){return X1s(z)*NH(z);}
    double NHeI1s(double z){return XHeI1s(z)*NH(z);}
    double NHeII1s(double z){return XHeII1s(z)*NH(z);}
    double NH(double z){return (1.0-Y_p)*Nb(z); }                                   // number density of hydrogen nuclei 1/cm^3
    double Np(double z){return Xp(z)*NH(z); }                                       // number density of free protons    1/cm^3
    //
    double NHe(double z){return Y_p/(4.0*fac_mHemH)*Nb(z);}                         // number density of helium nuclei   1/cm^3
    double NHeI(double z);                                                          // number density of HeI    1/cm^3
    double NHeII(double z);                                                         // number density of HeII   1/cm^3
    double NHeIII(double z);                                                        // number density of HeIII  1/cm^3
    //
    double fHe(double Ypv){ return Ypv/(4.0*fac_mHemH)/(1.0-Ypv); }                 // nHe/nH
    double fHe()          { return fHe(Y_p); }                                      // nHe/nH
    //
    double rho_g_g_cm3(double zz){return rho_g*pow(1.0+zz, 4);}                     // photon energy density g/cm^3
    double rho_g_1_cm3(double zz){return rho_g_gr*pow(1.0+zz, 4);}                  // photon energy density 1/cm^3 (mec2)
    double n_gamma(double zz){return 2.0*1.20206*1.7595e+30*pow(ThetaCMB(zz), 3);}  // photon number density in 1/cm^3

    double rho_tot(double z){ return 3.0*pow(H(z), 2)/(8.0*PI*const_G); }           // total energy density g/cm^3
    double rho_tot_eV(double z){ return rho_tot(z)*const_me*1.0e+3/const_me_gr; }   // total energy density eV/cm^3

    //====================================================================================
    // Hubble parameter 
    //====================================================================================
    double H(double z);
    
    //====================================================================================
    // different derivatives 
    //====================================================================================
    double dt_dz(double z){ return -1.0/H(z)/(1.0+z);}                        
    double dz_dt(double z){ return -H(z)*(1.0+z);}                        

    //====================================================================================
    // expansion time
    //====================================================================================
    double t_exp(double z){ return  1.0/H(z); } 
    
    //====================================================================================
    // cosmological time
    //====================================================================================
    double t(double z);                 

    //====================================================================================
    // times for energy exchange
    //====================================================================================
    double t_eg(double z){ return 1.0/(sigT(z)*const_cl*Ne(z));}                      // 1/(sigT ne c)
    double t_C(double z){ return t_eg(z);}                                            // 1/(sigT ne c)
    double t_K(double z){return t_eg(z)/(4.0*ThetaCMB(z));}                           // 1/(4 Theta sigT ne c)
    double dtau_C_dz(double z){ return dt_dz(z)/t_eg(z);}                             // dtau_C/dz 
    double t_ep(double z);                                                            // Spitzer formula
    double tau_C(double z);                                                           // Thomson optical depth
    double g_visi(double z);                                                          // Thomson visibility fcn (LH 2017)
    //-----------------------------------------------------------------------------------
    // 16.05.2017: JC
    //-----------------------------------------------------------------------------------
    double sigT(double z);                                                            // to allow Thomson cross sec to vary
    void reset_sigT_function();                                                       // sigT(z)=sigT_0=const
    void set_sigT_function(double aS, double mS, double p);                           // to set models for varying consts

    //====================================================================================
    // Integral^zh_z (1+z)^pow |dt_dz| dz 
    //====================================================================================
    double E_Int(double zh, double z, double pow);
    
    //====================================================================================
    // functions to compute energy release caused by acoustic damping (added 06.05.2017, JC)
    //====================================================================================
    double eta(double z){ return const_cl*E_Int(1.0e+16, z, 1); }                     // conformal time
    double deta_dz(double z){ return -const_cl/H(z); }                                // d eta/ dz 
    
    double Get_O_r() const  { return O_g*g_rel;   } 
    double Get_H0() const  { return H0;   } 
    
    double rho_rel(double z){return rho_c*O_g*g_rel*pow(1.0+z, 4); }                  // energy density of relativistic species g/cm^3
    double rho_gamma(double z){return rho_c*O_g*pow(1.0+z, 4); }                      // photon energy density g/cm^3
    double rho_Nu(double z){return rho_c*O_g*(g_rel-1.0)*pow(1.0+z, 4); }             // Neutrino energy density g/cm^3
    double rho_baryon(double z){return rho_c*O_b*pow(1.0+z,3); }                      // baryon energy density g/cm^3
    double rho_dm(double z){return rho_c*(O_m-O_b)*pow(1.0+z,3); }                    // dark matter energy density g/cm^3
    double rho_m(double z){return rho_c*O_m*pow(1.0+z,3); }                           // matter energy density g/cm^3
    
    double R_bg(double z){ return 3.0*rho_baryon(z)/4.0/rho_gamma(z); }               // 'baryon loading'
    double R_Nu(double z){ return rho_Nu(z)/rho_rel(z); }                             // 'Neutrino loading'
    double R_b(double z){ return rho_baryon(z)/rho_m(z); }                            // 'baryon fraction'
    double R_c(double z){ return rho_dm(z)/rho_m(z); }                                // 'CDM fraction'
    double c_s(double z){ return 1.0/sqrt(3.0*(1.0+R_bg(z))); }                       // 'sound speed' of photon + baryon fluid
    double c_sb(double z);                                                            // 'sound speed' of baryons
    double r_s(double z);                                                             // 'sound horizon'
    double Nnu_func(double Rnu){ return 11.0/7.0*pow(22.0, 1.0/3.0)*Rnu/(1.0-Rnu); }  // 'Neff'
    
    double dtau_deta(double z){ return sigT(z)*Ne(z)/(1.0+z); }                       // 'dtau_deta'
    double k_D(double z);                                                             // 'damping wavenumber'
    double k_lambda(double z);                                                        // 'damping wavenumber' for magnetic heating
    double dkm2_D_deta(double z);                                                     // eta-derivative of 1/kD^2
    double dkm2_D_dz(double z);                                                       // z-derivative of 1/kD^2
    double k_mfp(double z){ return sigT(z)*Ne(z)/(1.0+z); }                           // scale of comoving 'mean free path'
    double k_eq();                                                                    // k_eq=(2 Om H0^2 a0/aeq)^0.5

    //====================================================================================
    // y-parameters && free-free functions
    //====================================================================================
    double y_gamma(double zlow, double zup);
    double y_elect(double zlow, double zup);
    double tau_ff(double x, double z);
    double compute_z_for_given_y(double y_val);
    double compute_x_abs(double z, double val=1.0);
    double compute_x_h_low(double z);
    double compute_x_h_high(double z);
};

#endif
