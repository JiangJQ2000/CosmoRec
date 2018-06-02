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

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>
#include <stdlib.h>

#include "Cosmos.h"
#include "Recfast++.h"
#include "cosmology.Recfast.h"
#include "physical_consts.h"
#include "routines.h"
#include "Integration_routines.h"

using namespace std;

//========================================================================================
// pointer for GSL integrals (avoid type-casting)
//========================================================================================
Cosmos *GSL_cosmosptr=NULL;

//========================================================================================
// level of communication
//========================================================================================
int Cosmos_mess_flag=0;

//========================================================================================
//
// Class Hubble
//
//========================================================================================

//========================================================================================
// Konstructors and Destructors
//========================================================================================
Hubble :: Hubble(){ loaded=pointer=0; mem_index_Hz=-10; ext_Hptr=NULL; }

Hubble :: Hubble(const vector<double> &z, const vector<double> &Hz)
{ 
    loaded=pointer=0; mem_index_Hz=-10; ext_Hptr=NULL;
    init(z, Hz); 
}

Hubble :: Hubble(double (*Hptr)(const double *))
{
    loaded=pointer=1; mem_index_Hz=-10; ext_Hptr=Hptr;
}

Hubble :: ~Hubble(){ if(loaded) clear(); }

//========================================================================================
void Hubble :: init(const double *z, const double *Hz, const int nz)
{
    vector<double> lz(nz), lHz(nz);  
    
    // check ordering
    if(z[0]<z[1]) 
    {
        if(z[0]==0) lz[0]=log(1.0e-10);
        else lz[0]=log(z[0]);
        
        zmin=exp(lz[0]); 
        zmax=z[nz-1];
        
        for(int k=1; k<nz; k++) 
        {
            lz[k]=log(z[k]);
            lHz[k]=log(Hz[k]);            
        }
    }
    else 
    {
        if(z[nz-1]==0) lz[0]=log(1.0e-10);
        else lz[0]=log(z[nz-1]);
        
        zmin=exp(lz[0]); 
        zmax=z[0];
        
        for(int k=1; k<nz; k++) 
        {
            lz[k]=log(z[nz-1-k]);
            lHz[k]=log(Hz[nz-1-k]);            
        }
    }
    
    if(mem_index_Hz==-10) 
        mem_index_Hz=calc_spline_coeffies_JC(nz, &lz[0], &lHz[0], "Hubble:: Hz");
    
    else update_spline_coeffies_JC(mem_index_Hz, nz, &lz[0], &lHz[0], "Hubble:: Hz");
    
    loaded=1;
    lz.clear();
    lHz.clear();
    
    return;
}

void Hubble :: init(const vector<double> &z, const vector<double> &Hz)
{
    int nz=z.size();
    
    init(&z[0], &Hz[0], nz);
    
    return;
}

void Hubble :: init(double (*Hptr)(const double *))
{
    loaded=pointer=1; ext_Hptr=Hptr;
    return;
}

//========================================================================================
void Hubble :: clear()
{
    if(mem_index_Hz!=-10) free_spline_JC(mem_index_Hz);
    loaded=pointer=0; 
    mem_index_Hz=-10;
    ext_Hptr=NULL;
    
    return;
}

bool Hubble :: check_limits(double z)
{
    if(pointer) return 1;
    else if(loaded==1 && z>zmin && z<zmax) return 1;
    
    return 0;
}

//========================================================================================
double Hubble :: H(double z)
{
    if(pointer) return ext_Hptr(&z); 
    else if(loaded) return exp(calc_spline_JC(log(z), mem_index_Hz, " Hubble:: H(z)"));
        
    cerr << " Hubble:: Hubble was not initialized " << endl;
    
    return 0;
}
//========================================================================================



//========================================================================================
//
// Class Cosmos
//
//========================================================================================

//========================================================================================
// Konstructors and Destructors
//========================================================================================
void Cosmos :: init(const parameters &pars)
{
    h100=pars.h100; T_CMB0=pars.T0; Y_p=pars.Yp; Nnu=pars.Nnu;
    
    O_m=pars.Omega_m; O_b=pars.Omega_b; O_L=pars.Omega_L, O_k=pars.Omega_k;

    //====================================================================================
    // compute Omega_L assuming Omega_tot=1.0 if requested
    //====================================================================================
    if(O_L<=0.0) O_L=1.0-O_m-O_k-calc_Orel(T_CMB0, Nnu, h100);
    
    set_initials(Nnu);

    //====================================================================================
    // for Xe-spline interpolations...
    //====================================================================================
    calc_dependent_consts();
 
    splines_are_charged=0;
    
    //RF_hPlanck  =6.62606876E-34;      // Planck's constant [Js]    
    //RF_kBoltz   =1.3806503E-23;       // Boltzman's constant [J/K] 
    //RF_mElect   =9.10938188E-31;      // Electron mass [kg]  
    RF_hPlanck  =const_h*1.0e-7;      // Planck's constant [Js]    
    RF_kBoltz   =const_kB_J;          // Boltzman's constant [J/K] 
    RF_mElect   =const_me_gr*1.0e-3;  // Electron mass [kg]        

    //====================================================================================
    // to set default model for varying consts [JC]
    //====================================================================================
    reset_sigT_function();

    // Recfast++ cosmology
    CosmoRF.Yp=Y_p; CosmoRF.T0=T_CMB0; CosmoRF.h100=h100; CosmoRF.Neff=Nnu;
    CosmoRF.Omega_m=O_m; CosmoRF.Omega_b=O_b; CosmoRF.Omega_k=O_k; CosmoRF.Omega_L=O_L;

    // Recfast++ recombination settings
    RecRF=pars.Rec;
    // overwrite a few settings
    RecRF.verbosity=Cosmos_mess_flag;
    RecRF.npts=n_Xe; RecRF.zstart=zsRe; RecRF.zend=0.0;

    return;
}

void Cosmos :: init_Hubble(const vector<double> &z, const vector<double> &Hz)
{
    loaded_Hz.init(z, Hz);
    
    return;
}

void Cosmos :: init_Hubble(const double *z, const double *Hz, const int nz)
{
    loaded_Hz.init(z, Hz, nz);
    
    return;
}

void Cosmos :: init_Hubble(double (* Hptr)(const double *))
{
    loaded_Hz.init(Hptr);
    
    return;
}

//========================================================================================
Cosmos :: Cosmos(){ spline_memory_allocated=0; reset_sigT_function(); }
Cosmos :: Cosmos(const parameters &pars)
{ 
    spline_memory_allocated=0; 
    init(pars);
}

Cosmos :: ~Cosmos(){}

//========================================================================================
// private functions
//========================================================================================
void Cosmos :: set_initials(double Nnuval)
{
    g_rel=1.0+Nnuval*(7.0/8.0)*pow(4.0/11.0, 4.0/3.0);
    
    // parameters for recombination history
    n_Xe=6000;
    zsRe=2.5e+4;

    return;
}

void Cosmos :: init_splines()
{
    calc_coeff_X_spline(n_Xe, zsRe, 0.0);
    splines_are_charged=1;

    return;
}

void Cosmos :: calc_dependent_consts()
{
    h2=h100*h100; H0=100.0*h100*(1.0e+3*1.0e+2)/const_Mpc;  // H = h *100 * km/s/Mpc 
                                                            //   = h *100 * 10^5  cm/s/Mpc  
    T27=T_CMB0/2.7; Theta_CMB0=TtoTheta(T_CMB0); Theta27=T27;

    rho_c=3.0*pow(H0, 2)/(8.0*PI*const_G);                  //  1.87882e-29 h2    g cm^-3  

    rho_g=PI*PI/15.0*const_kB*T_CMB0                        //  4.47378e-34 T27^4 g cm^-3
          *pow(const_kB*T_CMB0/(const_hbar*const_cl), 3)/pow(const_cl, 2);

    rho_g_gr=rho_g/const_me_gr;
    
    // relativistic contributions
    O_g=rho_g/rho_c; O_gh2=O_g*h2;
    O_rel=O_g*g_rel; O_relh2=O_rel*h2;
    // matter contributions
    O_mh2=O_m*h2; O_bh2=O_b*h2;
    z_eq=O_mh2/O_relh2 -1.0;                     //  4.199628e+4/g_rel*pow(T27,-4)*O_mh2;
    // total and curvature
    O_T=O_m+O_L+O_k+O_rel;
    O_Th2=O_T*h2; O_Lh2=O_L*h2;
    O_kh2=O_k*h2;

    //====================================================================================
    // 28.05.2008
    // we changed the definition of Nb to Nb=rho_b/mH, where mH is the neutral HI mass
    //====================================================================================
    Nb0=O_b*rho_c/const_mH_gr;                   //  1.123xxe-5 cm^-3
    return;
}

//========================================================================================
int Cosmos :: recombination_history(int nzpts, double zs, double ze,
                                    vector<double> &zarr,
                                    vector<double> &Xe_H, vector<double> &Xe_He,
                                    vector<double> &Xe, vector<double> &TM)
{
    vector<double> dXe(nzpts), dX_H(nzpts);
    return recombination_history(nzpts, zs, ze, zarr, Xe_H, Xe_He, Xe, dXe, dX_H, TM);
}

//========================================================================================
double Hfunction(double z){ return GSL_cosmosptr->H(z); }

int Cosmos :: recombination_history(int nzpts, double zstart, double zend, 
                                    vector<double> &zarr,
                                    vector<double> &Xe_H, vector<double> &Xe_He,
                                    vector<double> &Xe, vector<double> &dXe,
                                    vector<double> &dX_H, vector<double> &TM)
{
    RecPars Rec=RecRF;
    Rec.npts=nzpts; Rec.zstart=zstart; Rec.zend=zend;

    if(Cosmos_mess_flag>=2)
    {
        cout << "\n Cosmos :: recombination_history: " << endl;
        cout << " running recombination_history with the following parameters: " << endl;   
        
        cout << " zs: " << Rec.zstart << "\t ze: " << Rec.zend
             << "\t nzpts: " << nzpts << endl;
        cout << " Y_p: " << Y_p << "\t TCMB: " << T_CMB0 << endl;
        cout << " Omega_m: " << O_m << "\t Omega_b: " << O_b
             << "\n Omega_L: " << O_L << "\t Omega_k: " << O_k << endl;
        cout << " Omega_gamma h2: " << O_gh2 << "\t Omega_rel h2: " << O_relh2 
             << " Nnu= " << Nnu << endl;
        cout << " Hubble constant in units of 100 km s-1 Mpc-1: " << h100 << endl;
        cout << " Fudge factor for H recombination: " << Rec.F << endl;
        cout << " DM-annihilation " << Rec.f_ann << endl;
        cout << " Chluba & Thomas correction function: " << Rec.include_CF << endl;
    }
    
    // call with hubble function of cosmology object!
    GSL_cosmosptr=this;
    Recfastpp_Cosmology::set_H_pointer(Hfunction);

    Xe_frac(CosmoRF, Rec, zarr, Xe_H, Xe_He, Xe, dXe, dX_H, TM);

    Recfastpp_Cosmology::reset_H_pointer();
    GSL_cosmosptr=NULL;
    
    return 0;
}

//========================================================================================
int Cosmos :: calc_coeff_X_spline(int M, double zstart, double zend)
{
    vector<double> zarray(M), zGSL(M), YGSL(M);
    vector<double> Xe_H(M), dXe_H(M), Xe_He(M), Xe(M), dXe(M), TM(M);
    
    if(Cosmos_mess_flag>=1)
        cout << "\n Cosmos::calc_coeff_X_spline:" 
             << "\n calculating coeffients for X spline interpolation..." << endl;
    
    //====================================================================================
    // here arrays are filled with values
    //====================================================================================
    recombination_history(M, zstart, zend, zarray, Xe_H, Xe_He, Xe, dXe, dXe_H, TM);
    
    //====================================================================================
    // resort zarray and Xe array so that ordering is zarray[0]<zarray[1]< ...
    //====================================================================================
    for(int j=0, i=M-1; i>=0; i--, j++) 
    {
        zGSL[j]=zarray[i];
        YGSL[j]=log(Xe[i]);
    }
    if(!spline_memory_allocated) 
        memindex_Xe=calc_spline_coeffies_JC(M, &zGSL[0], &YGSL[0],
                                            "Cosmos :: calc_coeff_X_spline: Xe");
    else 
        update_spline_coeffies_JC(memindex_Xe, M, &zGSL[0], &YGSL[0],
                                  "Cosmos :: calc_coeff_X_spline: Xe");

    //====================================================================================
    // dXe
    //====================================================================================
    for(int j=0, i=M-1; i>=0; i--, j++) YGSL[j]=dXe[i];
    if(!spline_memory_allocated) 
        memindex_dXe=calc_spline_coeffies_JC(M, &zGSL[0], &YGSL[0],
                                             "Cosmos :: calc_coeff_X_spline: dXe");
    else 
        update_spline_coeffies_JC(memindex_dXe, M, &zGSL[0], &YGSL[0],
                                  "Cosmos :: calc_coeff_X_spline: dXe");

    //====================================================================================
    // X1s
    //====================================================================================
    for(int j=0, i=M-1; i>=0; i--, j++) YGSL[j]=log(Xe_H[i]);
    if(!spline_memory_allocated) 
        memindex_XeH=calc_spline_coeffies_JC(M, &zGSL[0], &YGSL[0],
                                             "Cosmos :: calc_coeff_X_spline: XeH");
    else 
        update_spline_coeffies_JC(memindex_XeH, M, &zGSL[0], &YGSL[0],
                                  "Cosmos :: calc_coeff_X_spline: XeH");
                                   
    //====================================================================================
    // dXe_H
    //====================================================================================
    for(int j=0, i=M-1; i>=0; i--, j++) YGSL[j]=dXe_H[i];
    if(!spline_memory_allocated) 
        memindex_dXeH=calc_spline_coeffies_JC(M, &zGSL[0], &YGSL[0],
                                              "Cosmos :: calc_coeff_X_spline: dXeH");
    else 
        update_spline_coeffies_JC(memindex_dXeH, M, &zGSL[0], &YGSL[0],
                                  "Cosmos :: calc_coeff_X_spline: dXeH");
                                   
    //====================================================================================
    // rho=Te/Tg 
    //====================================================================================
    for(int j=0, i=M-1; i>=0; i--, j++) YGSL[j]=TM[i]/TCMB(zarray[i]);
    if(!spline_memory_allocated) 
        memindex_rho=calc_spline_coeffies_JC(M, &zGSL[0], &YGSL[0],
                                             "Cosmos :: calc_coeff_X_spline: rho");
    else 
        update_spline_coeffies_JC(memindex_rho, M, &zGSL[0], &YGSL[0],
                                  "Cosmos :: calc_coeff_X_spline: rho");
                                   
    //====================================================================================
    // XHe1s
    //====================================================================================
    for(int j=0, i=M-1; i>=0; i--, j++) YGSL[j]=log(max(Xe_He[i], 1.0e-20));
    if(!spline_memory_allocated) 
        memindex_XeHe=calc_spline_coeffies_JC(M, &zGSL[0], &YGSL[0],
                                              "Cosmos :: calc_coeff_X_spline: Xe_He");
    else 
        update_spline_coeffies_JC(memindex_XeHe, M, &zGSL[0], &YGSL[0],
                                  "Cosmos :: calc_coeff_X_spline: Xe_He");
    
    //====================================================================================
    // clean up
    //====================================================================================
    zarray.clear(); zGSL.clear(); YGSL.clear();
    Xe_H.clear(); dXe_H.clear(); Xe_He.clear(); Xe.clear(); dXe.clear(); TM.clear();

    spline_memory_allocated=1;

    return 0;
}

//========================================================================================
// local spline-functions
//========================================================================================
void Cosmos :: check_splines()
{   
    if(splines_are_charged) return; 
    else
    { 
        cerr << " Cosmos::check_splines: recombination splines are not armed."
             << " Exiting. " << endl; 
        exit(0); 
    }
}

double Cosmos :: calc_spline(double z, int memindex)    
{   
    check_splines(); 
    return calc_spline_JC(z, memindex, " Cosmos::calc_spline "); 
}
double Cosmos :: calc_spline_log(double z, int memindex)
{ 
    check_splines(); 
    return exp(calc_spline_JC(z, memindex, " Cosmos::calc_spline_log ")); 
}

//========================================================================================
// public functions
//========================================================================================
void Cosmos :: display_cosmology()
{
    cout << "\n*----------------------------------------------"
         << "--------------------------------------------*" << endl;
    cout << "* Cosmos :: These are the most important cosmological"
         << " parameters used in this computation: *" << endl;
    cout << "*------------------------------------------------"
         << "------------------------------------------*\n" << endl;

    cout << " Y_p: " << Y_p << "\t\t fHe: " << fHe() << "\t\t TCMB: " << T_CMB0 << endl;
    cout << " Omega_m   : " << O_m   << "\t Omega_b   : " << O_b
         << "\t Omega_L   : " << O_L   << "\t Omega_k   : " << O_k << endl;
    cout << " Omega_m h2: " << O_mh2 << "\t Omega_b h2: " << O_bh2
         << "\t Omega_L h2: " << O_Lh2 << "\t Omega_k h2: " << O_kh2 << endl;
    cout << " Omega_gamma h2: " << O_gh2 << "\t Omega_rel h2: " << O_relh2 
         << "\t Nnu= " << Nnu << endl;
    cout << " Hubble constant in units of 100 km s^-1 Mpc^-1: " << h100
         << " 1Mpc= " << const_Mpc << " cm " << endl;
    cout << " radiation-matter equality redshift z_eq: " << z_eq << endl;

    cout << "\n*----------------------------------------------"
         << "--------------------------------------------*" << endl;
    cout << "*------------------------------------------------"
         << "------------------------------------------*\n" << endl;
    
    return;
}

//========================================================================================
void Cosmos :: dump_cosmology(string filename)
{
    ofstream f;
    f.open(filename.c_str());
    f.precision(8);
    
    f << " Y_p: " << Y_p << "\t fHe: " << fHe() << "\t TCMB: " << T_CMB0 << endl;
    f << " Omega_m   : " << O_m   << "\t Omega_b   : " << O_b
      << "\t Omega_L   : " << O_L   << "\t Omega_k   : " << O_k << endl;
    f << " Omega_m h2: " << O_mh2 << "\t Omega_b h2: " << O_bh2
      << "\t Omega_L h2: " << O_Lh2 << "\t Omega_k h2: " << O_kh2 << endl;
    f << " Omega_gamma h2: " << O_gh2 << "\t Omega_rel h2: " << O_relh2 
      << "\t Nnu= " << Nnu << endl;
    f << " Hubble constant in units of 100 km s^-1 Mpc^-1: " << h100
      << " 1Mpc= " << const_Mpc << " cm " << endl;
    f << " radiation-matter equality redshift z_eq: " << z_eq << endl;

    f.close();

    return;
}

//========================================================================================
// Hubble-function
//========================================================================================
double Cosmos :: H(double z)
{
    if(loaded_Hz.check_limits(z)) return loaded_Hz.H(z); 
    
    //====================================================================================
    // This is equivalent with Seager 2000 ApJSS 128
    // return H0*sqrt(O_L+zp1*zp1*(O_k+zp1*O_m*(1.0+zp1/(1.0+z_eq)))); 
    // (Seager 2000 ApJSS 128)
    //====================================================================================
    double zp1=1.0+z;
    return H0*sqrt(O_L+zp1*zp1*(O_k+zp1*(O_m+O_rel*zp1)));       
}

//========================================================================================
// cosmological time
//========================================================================================
double dtdz_GSL(double logz, void *p)
{ 
    Cosmos *c =(Cosmos *) p;
    double zz=exp(logz); 
    return -zz*c->dt_dz(zz-1.0);
}

double Cosmos :: t(double z)
{
    double epsInt=1.0e-12, epsabs=0.0, a=log(1.0+z), b=log(1.0+1.0e+15), r=0.0;
    void *p=this;
    
    r=Integrate_gk_GSL(6, a, b, epsInt, epsabs, dtdz_GSL, p);  

    return r;
}

//========================================================================================
// Spitzer formula
//========================================================================================
double Cosmos :: t_ep(double z)                                                
{ 
    return t_eg(z)*sqrt(PI/2.0)*const_mp/(20.0*const_me)*(1.0-Y_p/2.0)
                  /(1.0-Y_p)*pow(TCMB(z)*(1.0+const_me/const_mp), 1.5); 
}

//========================================================================================
// integrals over Hubble; Integral^zh_z (1+z)^pow |dt_dz| dz
//========================================================================================
struct Integral_Info{ Cosmos *c; double var; };

double dE_GSL(double logz, void *p)
{ 
    Integral_Info *d =(Integral_Info *) p;
    double z1=exp(logz); 
    return -z1*pow(z1, d->var)*d->c->dt_dz(z1-1.0);
}

double Cosmos :: E_Int(double zh, double z, double pow)
{
    Integral_Info V; V.c=this; V.var=pow;
    double epsInt=1.0e-12, epsabs=0.0, a=log(z+1.0), b=log(zh+1.0);
    void *p=&V;

    return Integrate_gk_GSL(6, a, b, epsInt, epsabs, dE_GSL, p);
}

//========================================================================================
// Saha equations and solutions [restructured Dec 2017, JC]
// main assumption: only ground-state is signifantly populated
//========================================================================================
double Cosmos :: SahaBoltz_HeIII(double T)
{
    // HeII ionization energy in [J] --> nu = 1.315817072 10^16 Hz
    //double RF_EionHeII =8.7186944e-18;
    double RF_EionHeII =4.0*const_EH_inf_ergs*1.0e-7/(1.0+const_me_malp);
    
    double c1 = 2.0*PI*RF_mElect*RF_kBoltz*T/(RF_hPlanck*RF_hPlanck);
    double NHm3=NH(Ttoz(T))*1.0e+6;
    double g = pow(c1, 1.5)*exp(-RF_EionHeII/(RF_kBoltz*T))/NHm3, d=1.0+fHe()-g;
    double Xe= 0.5*(d + sqrt(d*d + 4.0*g*(1.0+2.0*fHe())));
    
    return Xe;
}

double Cosmos :: SahaBoltz_HeII(double T)
{
    double RF_EionHeI  =3.93933e-18;         // HeI ionization energy in [J]
    
    double c1 = 2.0*PI*RF_mElect*RF_kBoltz*T/(RF_hPlanck*RF_hPlanck);
    double NHm3=NH(Ttoz(T))*1.0e+6;
    double g = 4.0*pow(c1, 1.5)*exp(-RF_EionHeI/(RF_kBoltz*T))/NHm3, d=1.0-g;
    double Xe= 0.5*(d + sqrt(d*d + 4.0*g*(1.0+fHe())));
    
    return Xe;
}

//========================================================================================
// these solutions just assume that Xe is given
//========================================================================================
double Cosmos :: SahaBoltz_HeII1s(double T)
{
    double RF_EionHeI  =3.93933e-18;         // HeI ionization energy in [J]
    double RF_EionHeII =4.0*const_EH_inf_ergs*1.0e-7/(1.0+const_me_malp);
    
    double c1 = 2.0*PI*RF_mElect*RF_kBoltz*T/(RF_hPlanck*RF_hPlanck);
    double NHm3=NH(Ttoz(T))*1.0e+6;
    double gHeI  = 4.0*pow(c1, 1.5)*exp(-RF_EionHeI /(RF_kBoltz*T))/NHm3;
    double gHeII =     pow(c1, 1.5)*exp(-RF_EionHeII/(RF_kBoltz*T))/NHm3;
    double Xe=Xe_Seager(Ttoz(T));
    double XHeI1s=fHe()*gHeI*Xe/(gHeI*(Xe+gHeII)+Xe*Xe);
    
    return XHeI1s;
}

double Cosmos :: SahaBoltz_HeI1s(double T)
{
    double RF_EionHeI  =3.93933e-18;         // HeI ionization energy in [J]
    double RF_EionHeII =4.0*const_EH_inf_ergs*1.0e-7/(1.0+const_me_malp);

    double c1 = 2.0*PI*RF_mElect*RF_kBoltz*T/(RF_hPlanck*RF_hPlanck);
    double NHm3=NH(Ttoz(T))*1.0e+6;
    double gHeI  = 4.0*pow(c1, 1.5)*exp(-RF_EionHeI /(RF_kBoltz*T))/NHm3;
    double gHeII =     pow(c1, 1.5)*exp(-RF_EionHeII/(RF_kBoltz*T))/NHm3;
    double Xe=Xe_Seager(Ttoz(T));
    double XHeI1s=fHe()*Xe*Xe/(gHeI*(Xe+gHeII)+Xe*Xe);
    
    return XHeI1s;
}

double Cosmos :: SahaBoltz_HI1s(double T)
{
    // HeI ionization energy in [J]
    //double RF_EionHI   =2.17869e-18;
    double RF_EionHI   =const_EH_inf_ergs*1.0e-7/(1.0+const_me_mp);
    
    double c1 = 2.0*PI*RF_mElect*RF_kBoltz*T/(RF_hPlanck*RF_hPlanck);
    double NHm3=NH(Ttoz(T))*1.0e+6;
    double g = pow(c1, 1.5)*exp(-RF_EionHI/(RF_kBoltz*T))/NHm3;
    double Xe=Xe_Seager(Ttoz(T));
    double X1s=Xe/(Xe+g);
    
    return X1s;
}

//========================================================================================
// correction to Saha problem:
// full equation: Xe == gH/(gH+Xe) + fHe*gHeI*(2.0*gHeII+Xe)/(gHeI*(gHeII+Xe)+Xe*Xe)
// zeroth order solution is given by Seager solution. Then correction is obtained here
//========================================================================================
double Cosmos :: XeSaha_corr(double z)
{
    if(z<3000.0) return 0.0;
    
    double T=TCMB(z);
    double RF_EionHI   =const_EH_inf_ergs*1.0e-7/(1.0+const_me_mp);
    double RF_EionHeI  =3.93933e-18;         // HeI ionization energy in [J]
    double RF_EionHeII =4.0*const_EH_inf_ergs*1.0e-7/(1.0+const_me_malp);

    double c1 = 2.0*PI*RF_mElect*RF_kBoltz*T/(RF_hPlanck*RF_hPlanck);
    double NHm3=NH(z)*1.0e+6;
    
    double gHI   =     pow(c1, 1.5)*exp(-RF_EionHI  /(RF_kBoltz*T))/NHm3;
    double gHeI  = 4.0*pow(c1, 1.5)*exp(-RF_EionHeI /(RF_kBoltz*T))/NHm3;
    double gHeII =     pow(c1, 1.5)*exp(-RF_EionHeII/(RF_kBoltz*T))/NHm3;
    
    double Xe=Xe_Seager(z);  // this is Xe^(0)
    
    // compute correction terms
    double h1=gHI/(gHI+Xe);
    double h2=fHe()*gHeI*(2.0*gHeII+Xe)/(gHeI*(gHeII+Xe)+Xe*Xe);
    double D1=Xe-h1-h2;
    
    double g1=h1/(gHI+Xe);
    double g2=fHe()*gHeI*(gHeI*gHeII+(4.0*gHeII+Xe)*Xe)/pow(gHeI*(gHeII+Xe)+Xe*Xe, 2);
    
    return D1*(1.0-(g1+g2)*(1.0-(g1+g2)));
}

//========================================================================================
double Cosmos :: XHeII1s(double z)
{
    if(z>=3500.0) return SahaBoltz_HeII1s(TCMB(z));
    return calc_spline_log(z, memindex_XeHe);         // assumption hardly any NHeIII left
}

double Cosmos :: XHeI1s(double z)
{
    if(z>=3500.0) return SahaBoltz_HeI1s(TCMB(z));
    return fHe()-calc_spline_log(z, memindex_XeHe);
}

double Cosmos :: X1s(double z)
{
    if(z>=3500.0) return SahaBoltz_HI1s(TCMB(z));
    return 1.0-calc_spline_log(z, memindex_XeH);
}

//========================================================================================
// helium nuclei (assumption: only ground-state is signifantly populated)
//========================================================================================
double Cosmos :: NHeIII(double z)
{ return (z>=3500.0 ? (fHe()-XHeI1s(z)-XHeII1s(z))*NH(z) : 1.0e-20*NH(z)); }

double Cosmos :: NHeII(double z){ return XHeII1s(z)*NH(z); }
double Cosmos :: NHeI(double z){ return XHeI1s(z)*NH(z); }

//========================================================================================
// compute Recfast system with a given initial solution at low z
//========================================================================================
int Cosmos :: recombine_using_Recfast_system(int nzpts, double zi, double ze, double fDM,
                                             double Xe_Hi, double Xe_Hei, double Xei, double TMi,
                                             double dXei, const vector<double> &zarr,
                                             vector<double> &Xe_H, vector<double> &Xe_He,
                                             vector<double> &Xe, vector<double> &TM)
{
    RecPars Rec=RecRF;
    Rec.npts=nzpts; Rec.zstart=zi; Rec.zend=ze;
    Rec.f_ann=fDM;
    Rec.include_CF=0;
    
    if(Cosmos_mess_flag>=2)
    {
        cout << "\n Cosmos :: recombination_history: " << endl;
        cout << " running recombination_history with the following parameters: " << endl;
        
        cout << " zs: " << Rec.zstart << "\t ze: " << Rec.zend
             << "\t nzpts: " << nzpts << endl;
        cout << " Y_p: " << Y_p << "\t TCMB: " << T_CMB0 << endl;
        cout << " Omega_m: " << O_m << "\t Omega_b: " << O_b
             << "\n Omega_L: " << O_L << "\t Omega_k: " << O_k << endl;
        cout << " Omega_gamma h2: " << O_gh2 << "\t Omega_rel h2: " << O_relh2
             << " Nnu= " << Nnu << endl;
        cout << " Hubble constant in units of 100 km s-1 Mpc-1: " << h100 << endl;
        cout << " Fudge factor for H recombination: " << Rec.F << endl;
        cout << " DM-annihilation " << Rec.f_ann << endl;

        cout << "\n and initial values " << Xe_Hi << " " << Xe_Hei << " " 
             << Xei << " " << TMi << endl << endl;
    }
    
    // call with hubble function of cosmology object!
    GSL_cosmosptr=this;
    Recfastpp_Cosmology::set_H_pointer(Hfunction);
    
    Xe_frac_rescaled(CosmoRF, Rec, zarr, Xe_H, Xe_He, Xe, TM, Xe_Hi, Xe_Hei, Xei, TMi, dXei);

    Recfastpp_Cosmology::reset_H_pointer();
    GSL_cosmosptr=NULL;
    
    return 0;
}

//====================================================================================
// functions to compute energy release caused by acoustic damping (added 15.10.2011)
//====================================================================================
//double eff_fact_pol=8.0/9.0;
double eff_fact_pol=16.0/15.0;

double dr_s_GSL(double logz, void *p)
{ 
    Cosmos *c =(Cosmos *) p;
    double z1=exp(logz), z=z1-1.0; 
    return z1*const_cl/c->H(z)*c->c_s(z);
}

double Cosmos :: r_s(double z) // 'sound horizon'
{
    double epsInt=1.0e-7, epsabs=0.0, a=log(z+1.0), b=log(1.0e+16);
    void *p=this;

    return Integrate_using_Patterson_adaptive(a, b, epsInt, epsabs, dr_s_GSL, p);
}

//====================================================================================
double dk_D_GSL(double logz, void *p)
{ 
    Cosmos *c =(Cosmos *) p;
    double z1=exp(logz), z=z1-1.0;
    double R=c->R_bg(z);
    double tau_d=c->dtau_deta(z);
    return z1*const_cl/c->H(z)/(6.0*(1.0+R))*(R*R/(1.0+R)+eff_fact_pol)/tau_d;
}

double Cosmos :: k_D(double z) // 'damping wavenumber'
{
    double epsInt=1.0e-7, epsabs=0.0, a=log(z+1.0), b=log(1.0e+16), r=0.0;
    void *p=this;

    r=Integrate_using_Patterson_adaptive(a, b, epsInt, epsabs, dk_D_GSL, p);
    
    return sqrt(1.0/r);
}

//====================================================================================
double dk_lambda_GSL(double logz, void *p)
{
    Cosmos *c =(Cosmos *) p;
    double z1=exp(logz), z=z1-1.0;
    double tau_d=c->dtau_deta(z);
    //double R=c->R_bg(z);
    //if(z>3000) return z1*const_cl/c->H(z)/tau_d*R/(1.0+R); // == int dt / Ne sigT a^2 x R/(1+R)
    return z1*const_cl/c->H(z)/tau_d;           // == int dt / Ne sigT a^2
}

double Cosmos :: k_lambda(double z) // 'damping wavenumber' for magnetic stuff
{
    double epsInt=1.0e-7, epsabs=0.0, a=log(z+1.0), b=log(1.0e+16), r=0.0;
    void *p=this;
    
    r=Integrate_using_Patterson_adaptive(a, b, epsInt, epsabs, dk_lambda_GSL, p);
    
    return sqrt(1.0/r);
}

//====================================================================================
double Cosmos :: dkm2_D_deta(double z) // eta-derivative of 1/kD^2
{
    double R=R_bg(z);
    double tau_d=dtau_deta(z);
    
    return 1.0/(6.0*(1.0+R))*(R*R/(1.0+R)+eff_fact_pol)/tau_d;
}

double Cosmos :: dkm2_D_dz(double z) // z-derivative of 1/kD^2
{
    double R=R_bg(z);
    double tau_d=dtau_deta(z);
    
    return const_cl/H(z)/(6.0*(1.0+R))*(R*R/(1.0+R)+eff_fact_pol)/tau_d;
}

//====================================================================================
// 'sound speed' of baryons
//====================================================================================
double Cosmos :: dTe_dz(double z)
{ 
    double dz=z*1.0e-6;
    return ( Te(z+dz)-Te(z-dz) )/(2.0*dz); 
}

double Cosmos :: c_sb(double z)
{
    double Xev=Xe_Seager(z);   // Xe_b --> Xe_Seager, JC Dec 2017
    double mu_mol=const_mp_gr/(1.0-Y_p)/(1.0+fHe()+Xev);
    double Tb=Te(z);
    double cs2=const_kB*Tb/(const_cl*const_cl)/mu_mol;
    
    cs2*=1.0+z/Tb/3.0*dTe_dz(z);  // Ma & Bertschinger 1995
//    cs2*=1.0;                     // Lewis 2007 (10.1103/PhysRevD.76.063001), with Delta_Tg ~ 0
//    cs2*=5.0/3.0;
//    cout << " Cosmos :: c_sb " << 1.0+z/Tb/3.0*dTe_dz(z) << " " << z/Tb/3.0*dTe_dz(z) << endl;

    return sqrt(cs2);
}
    
double Cosmos :: k_eq()
{
    // k_eq=(2 Om H0^2 a0/aeq)^0.5
    double a0_aeq=Get_O_m()/Get_O_r()-1.0;
    
    return sqrt(2.0*Get_O_m()*a0_aeq)*Get_H0()/const_cl;
}

//========================================================================================
// y-parameters
//========================================================================================
double dy_gamma_GSL(double logz, void *p)
{
    Cosmos *c =(Cosmos *) p;
    double z1=exp(logz), z=z1-1.0;
    double dtau_dt=1.0/c->t_eg(z);
    double Thg=c->ThetaCMB(z);
    return -z1*Thg*dtau_dt*c->dt_dz(z);
}

double Cosmos :: y_gamma(double zlow, double zup)
{
    double epsInt=1.0e-5, epsabs=0.0, a=log(zlow+1.0), b=log(zup+1.0);
    void *p=this;

    return Integrate_gk_GSL(6, a, b, epsInt, epsabs, dy_gamma_GSL, p);
}

double dy_elect_GSL(double logz, void *p)
{
    Cosmos *c =(Cosmos *) p;
    double z1=exp(logz), z=z1-1.0;
    double dtau_dt=1.0/c->t_eg(z);
    double The=c->ThetaCMB(z)*c->Te_Tg(z);
    return -z1*The*dtau_dt*c->dt_dz(z);
}

double Cosmos :: y_elect(double zlow, double zup)
{
    double epsInt=1.0e-5, epsabs=0.0, a=log(zlow+1.0), b=log(zup+1.0);
    void *p=this;

    return Integrate_gk_GSL(6, a, b, epsInt, epsabs, dy_elect_GSL, p);
}

//========================================================================================
// y-parameter has certain value
//========================================================================================
struct Search_Data{ double z, val; Cosmos *c; };

double func_compute_y(double *lz, void *p)
{
    Search_Data *d=((Search_Data *)p);
    double z=exp(*lz);
    double G=d->c->y_gamma(0.0, z)/d->val;

    return G-1.0;
}

double Cosmos :: compute_z_for_given_y(double y_val)
{
    Search_Data d={0, y_val, this};
    void *p=&d;
    double lz1=log(1.0), lz2=log(1.0e+6), lzacc=1.0e-8;
    return exp(find_root_brent(func_compute_y, p, lz1, lz2, lzacc));
}

//========================================================================================
// free-free optical depth
//========================================================================================
double dtau_ff_GSL(double logz, void *p)
{
    Integral_Info *d =(Integral_Info *) p;
    double z1=exp(logz), z=z1-1.0;
    double x=d->var;
    
    double Thg=d->c->ThetaCMB(z);
    double Np=d->c->Np(z), NHeII=d->c->NHeII(z), NHeIII=d->c->NHeIII(z);
    double gff=( x<0.37 ? sqrt(3.0)/PI*log(2.25/x) : 1.0);
    double KBR=const_alpha/TWOPI*pow(const_lambdac, 3)/sqrt(3.0*TWOPI)/pow(Thg, 3.5)*gff*(Np+NHeII+4.0*NHeIII);
    double dtau_dt=1.0/d->c->t_eg(z);
    
    return -z1*dtau_dt*d->c->dt_dz(z)*KBR*(1.0-exp(-x))/pow(x, 3);
}

double dtau_ff_GSL_e(double logz, void *p)
{
    Integral_Info *d =(Integral_Info *) p;
    double z1=exp(logz), z=z1-1.0;
    double x=d->var;
    
    double Thg=d->c->ThetaCMB(z);
    double The=Thg*d->c->Te_Tg(z);
    double Np=d->c->Np(z), NHeII=d->c->NHeII(z), NHeIII=d->c->NHeIII(z);
    //double yg=d->c->y_gamma(z, glob_ygzs);
    //double xe=x*Thg/The*exp(yg);
    double xe=x*Thg/The;
    
    //if(z<=10.0){ Np=Ne=d->c->NH(z); The=d->c->ThetaCMB(3600.0); }
    
    //double gff=( xe<0.37 ? sqrt(3.0)/PI*log( 2.25/xe ) : 1.0);                    // Hu1995
    double gff=1.0+log(1.0+exp(sqrt(3.0)/PI*( log(2.25/xe)+0.5*log(The) )+1.425));  // Draine 2011
    //double gff2=1.0+log(1.0+exp(sqrt(3.0)/PI*( log(1.125/xe)+0.5*log(The) )+1.425));// Draine 2011
    //double gffN=(Np+NHeII)*gff+4.0*NHeIII*gff2;
    
    double KBR=const_alpha/TWOPI*pow(const_lambdac, 3)/sqrt(3.0*TWOPI)/pow(The, 3.5)*gff*(Np+NHeII+4.0*NHeIII);
    //double KBR=const_alpha/TWOPI*pow(const_lambdac, 3)/sqrt(3.0*TWOPI)/pow(The, 3.5)*gffN;
    double dtau_dt=1.0/d->c->t_eg(z);
    
    //return -z1*dtau_dt*d->c->dt_dz(z1-1.0)*KBR/pow(xe, 2);
    return -z1*dtau_dt*d->c->dt_dz(z1-1.0)*KBR*(1.0-exp(-xe))/pow(xe, 3);
}

double Cosmos :: tau_ff(double x, double z)
{
    Integral_Info V; V.c=this; V.var=x;
    double epsInt=1.0e-5, epsabs=0.0, a=log(1.0), b=log(z+1.0);
    void *p=&V;

    return Integrate_gk_GSL(6, a, b, epsInt, epsabs, dtau_ff_GSL_e, p);
}

//========================================================================================
// absorption frequency
//========================================================================================
double func_compute_x_abs(double *lxabs, void *p)
{
    Search_Data *d=((Search_Data *)p);
    double x=exp(*lxabs);
    double G=d->c->tau_ff(x, d->val)+1.2e-21*pow(1.0+d->val, 3)*exp(-2.0*x)/x/x;
    
    //double sr2yg=sqrt(2.0*d->c->y_gamma(0.0, d->val));
    //double F=Dawson_Int(sr2yg);
    //double G=0*d->c->tau_ff(x, d->val)+1.9e-6*(sr2yg-F)/x/x;
    
    return G/d->z-1.0;
}

double Cosmos :: compute_x_abs(double z, double val)
{
    Search_Data d={z, val, this};
    void *p=&d;
    double lx1=log(1.0e-10), lx2=log(1.0), lxacc=1.0e-8;
    
    return exp(find_root_brent(func_compute_x_abs, p, lx1, lx2, lxacc));
}

//========================================================================================
// absorption + 2y
//========================================================================================
double func_compute_x_h_low(double *lxabs, void *p)
{
    Search_Data *d=((Search_Data *)p);
    double x=exp(*lxabs);
    //double G=d->c->tau_ff(x, d->z)+1.2e-21*pow(1.0+d->z, 3)*exp(-2.0*x)/x/x;
    double Gt=d->c->tau_ff(x, d->z);
    double Gy=d->c->y_gamma(0.0, d->z);
    
    return 2.0*Gy/Gt-1.0;
}

double Cosmos :: compute_x_h_low(double z)
{
    Search_Data d={z, 0, this};
    void *p=&d;
    double lx1=log(1.0e-10), lx2=log(1.0), lxacc=1.0e-8;
    
    return exp(find_root_brent(func_compute_x_h_low, p, lx1, lx2, lxacc));
}

//========================================================================================
// heating at high frequencies
//========================================================================================
double func_compute_x_h_high(double *lxi, void *p)
{
    Search_Data *d=((Search_Data *)p);
    double x=exp(*lxi);
    double y=d->c->y_gamma(0.0, d->z);
    double alpha=(3.0-2.0/(1.0+x))/sqrt(1.0+x*y);
    double beta=1.0/(1.0+x*x*y/(1.0+x*y));
    //    double alpha=(3.0-2.0*exp(-x))/sqrt(1.0+x*y);
    //    double beta=1.0/(1.0+x*y*( 1.0-exp(-x) ));
    
    return exp(y*(alpha+beta))/(1.0+x*y)-1.0;
}

double Cosmos :: compute_x_h_high(double z)
{
    Search_Data d={z, 0, this};
    void *p=&d;
    double lx1=log(1.0e-10), lx2=log(100.0), lxacc=1.0e-8;
    
    return exp(find_root_brent(func_compute_x_h_high, p, lx1, lx2, lxacc));
}

//========================================================================================
// Thomson visibility function (LH 2017)
//========================================================================================
// to allow Thomson cross section to vary
double Cosmos :: sigT(double z){ return const_sigT_0*pow(ascale/mscale, 2)*pow((1.0+z)/1100.0, pval);}

void Cosmos :: reset_sigT_function(){ ascale=mscale=1.0; pval=0.0; return; } // --> sigT(z)=sigT_0=const
void Cosmos :: set_sigT_function(double aS, double mS, double p){ ascale=aS; mscale=mS; pval=p; return; }

double dtau_dz_GSL(double logz,void *p)
{
    Cosmos *c =(Cosmos *) p;
    double z1=exp(logz);
    return -z1*c->dtau_C_dz(z1-1.0);
}

double Cosmos :: tau_C(double z)
{
    double epsInt=1.0e-5, epsabs=0.0, i=log(1.0), f=log(1.0+z);
    void *p=this;

    return Integrate_gk_GSL(6, i, f, epsInt, epsabs, dtau_dz_GSL, p);
}

double Cosmos::g_visi(double z)
{
    // g(z) = dtau/dz * exp(-tau(z))
    return -dtau_C_dz(z)*exp(-tau_C(z));
}
//========================================================================================
//========================================================================================
