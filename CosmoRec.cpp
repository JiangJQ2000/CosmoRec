//===========================================================================================================
//
// Author: Jens Chluba with contributions from Rajat Thomas, Luke Hart, Vivian Poulin
// first implementation: April 2010
// last modification   : May   2018
// purpose: to access main CosmoRec code
//
//===========================================================================================================
// 20.05.2018: added parser for parameters reading [JC]
// 01.03.2018: added decaying particle cases [JC]
// 31.03.2017: added VFC control for batch mode / CAMB runs [LH & JC]
// 15.02.2017: added controll for magnetic field modes [JC]
// 22.01.2017: added correction function computation as standard output [LH]
// 22.01.2017: added Recfast++ variation of fundamental constants setup [LH & JC]
// 14.11.2014: control of magnetic field parameters is now all in CosmoRec.cpp (also Recfast++ part...)
// 17.10.2014: conversion from B_lambda to B0 when calling in batch mode
// 11.10.2014: added parameters for magnetic field heating
// 09.10.2014: added magnetic field heating module
// 23.07.2014: added possibility to change A2s1s two-photon decay rate
// 22.07.2014: new He-diffusion runmodes (runpars[1]== 4, 5, 6)
// 17.07.2014: included HeI-diffusion correction (most of the coding done 02-05.2011)
// 16.07.2014: default runmode for batch version with nShellsHeI=2
// 17.06.2013: small cosmetic changes
// 18.06.2012: added C-version to communicate Hubble pointer
// 08.06.2012: added version to communicate Hubble
// 23.03.2012: added runmode with helium radiative transfer module;
// 16.02.2012: added module to load precomputed recombination history (useful for runs with Camb)

//===========================================================================================================
// Standards
//===========================================================================================================
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <limits.h>
#include <vector>

//===========================================================================================================
// GSL
//===========================================================================================================
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort_double.h>

//===========================================================================================================
// several libs
//===========================================================================================================
#include "CosmoRec.h"
#include "Atom.h"
#include "HeI_Atom.h"

#include "Photoionization_cross_section.h"
#include "Sobolev.h"
#include "Cosmos.h"

#include "physical_consts.h"
#include "routines.h"
#include "parser.h"

#include "Quadrupole_lines_HI.h" 

//===========================================================================================================
// net rates for Helium case
//===========================================================================================================
#include "ODEdef_Xi.HeI.h"

//===========================================================================================================
// PDE routines (Raman-scattering, two-gamma corrections & Diffusion stuff)
//===========================================================================================================
#include "Solve_PDEs.h"
#include "Solve_PDEs_HeI.h"

//===========================================================================================================
// ODE definitions for effective HI & HeI atom and the CosmoRec ODE solver
//===========================================================================================================
#include "ODE_effective.h"
#include "ODE_solver_Rec.h"

//===========================================================================================================
// access to effective rates for HI & HeI atom
//===========================================================================================================
#include "get_effective_rates.HI.h"
#include "get_effective_rates.HeI.h"

//===========================================================================================================
// namespaces that are used
//===========================================================================================================
using namespace std;
using namespace ODE_HI_effective;
using namespace ODE_solver_Rec;

//===========================================================================================================
//
// define default HI A2s1s two-photon decay rate. This was added to allow changing it globally.
// The corresponding variables are defined as extern in './Development/Definitions/physical_consts.h'
// and './Development/Cosmology/Recfast++/include/constants.Recfast.h'. A new value is read from
// the parameter input file (see './Modules/setup_functions.cpp').
//
// (added 23.07.2014 JC)
//===========================================================================================================
double const_HI_A2s_1s = const_HI_A2s_1s_0;                                  // initial setting of the values
double const_sigT      = const_sigT_0;
double const_HeI_A2s_1s = const_HeI_A2s_1s_0;

//===========================================================================================================
// this is where the main parameters are handled and structures with settings are defined. If you want to
// add parameters or settings, this is the place to start [JC, May 2018]
//===========================================================================================================
#include "./Modules/global_variables.cpp"
#include "./Modules/global_functions.cpp"

//===========================================================================================================
// these are several Modules (i.e. plugins) that are directly loaded when compiling.
//===========================================================================================================
#include "./Modules/HI_routines.cpp"
#include "./Modules/HeI_routines.cpp"
#include "./Modules/HeI_feedback.cpp"
#include "./Modules/DPesc_HI_abs_tabulation.cpp"
#include "./Modules/Diffusion_correction.cpp"   
#include "./Modules/Diffusion_correction.HeI.cpp" 

//===========================================================================================================
// functions to read parameter file and set up memory etc
//===========================================================================================================
#include "./Modules/setup_functions.cpp"

//===========================================================================================================
// DM annihilation module
//===========================================================================================================
#include "./Modules/DM_annihilation_decay.cpp"

//===========================================================================================================
// Magnetic field heating module (added Oct. 2014)
//===========================================================================================================
#include "./Modules/Magnetic_fields.cpp"

//===========================================================================================================
// variation of fundamental constants [Oct 2016 - Jan 2017: LH & JC]
//===========================================================================================================
#include "./Modules/Variations_constants.cpp"

//===========================================================================================================
// PCA functions [JC, Aug 2017]
//===========================================================================================================
#include "./Modules/PCA_functions.cpp"

//===========================================================================================================
// main definitions for CosmoRec-ODE system
//===========================================================================================================
#include "./Modules/ODEdef_CosmoRec.cpp"

//===========================================================================================================
// writing the output
//===========================================================================================================
#include "./Modules/output_methods.cpp"


//===========================================================================================================
//
// Normal Recfast-run
//
//===========================================================================================================
void display_Recfast_pp_parameters(string filename, const Parameters_form_initfile &inp)
{
    if(filename!="") cout << "\n Using parameters corresponding to " << filename << endl;
    cout << "\n zs: " << inp.zstart << "\t ze: " << inp.zend << "\t nzpts: " << inp.nz << endl;
    cout << " Y_p: " << inp.cosm.Yp << "\t TCMB: " << inp.cosm.T0 << endl;
    cout << " Omega_m: " << inp.cosm.Omega_m << "\t Omega_b: " << inp.cosm.Omega_b
         << "\n Omega_L: " << inp.cosm.Omega_L << "\t Omega_k: " << inp.cosm.Omega_k
         << "\t Omega_rel: " << cosmos.calc_Orel(inp.cosm.T0, inp.cosm.Nnu, inp.cosm.h100)
         << " Nnu= " << inp.cosm.Nnu << endl;
    cout << " Hubble constant in units of 100 km s-1 Mpc-1: " << inp.cosm.h100 << endl;
    cout << " Fudge factor for H recombination: " << inp.cosm.Rec.F << endl;
    cout << " HI A2s1s rate: " << inp.cosm.Rec.A2s1s << endl;

    return;
}

//===========================================================================================================
// Just Recfast++
//===========================================================================================================
void run_Recfast_pp(string filename)
{
    //=======================================================================================================
    // initialize cosmology and the general startup
    //=======================================================================================================
    read_entries_from_parameter_file(filename, parameters);
    check_and_set_params(parameters, parameters.verbosity);

    if(parameters.verbosity>=1) display_Recfast_pp_parameters(filename, parameters);
    
    set_cosmology(parameters.cosm, cosmos);
    
    //=======================================================================================================
    // initialise zarr
    //=======================================================================================================
    double *za=new double[parameters.nz];
    init_xarr(parameters.zstart, max(1.0e-8, parameters.zend), za, parameters.nz, 0, 0);
    
    if(parameters.verbosity>=0) cout << " writing data to " << parameters.path << endl;
    
    ofstream ofile;
    string name=parameters.path;
    if(parameters.cosm.Rec.include_CF==1) name+="Xe_Recfast.corr_fac"+parameters.addname;
    else name+="Xe_Recfast"+parameters.addname;
    ofile.open(name.c_str());
    ofile.precision(10); 
    
    ofile << "# Recfast++ output. The columns are z, Xe=Ne/NH, Te\n#" << endl;

    for(int l=0; l<parameters.nz; l++) ofile << scientific << za[l] << " "
                                             << cosmos.Xe_Seager(za[l]) << " " 
                                             << cosmos.Te(za[l]) << endl;
    
    if(parameters.verbosity>=0) cout << " finished Recfast++ run." << endl;

    //=======================================================================================================
    // write visibility function given the ionization history (LH)
    //=======================================================================================================
    if (write_RF_visi)
    {
        if(parameters.verbosity>=0) cout << " writing visibility function output ..." << flush;
        
        ofstream vfile;
        string vname=parameters.path;
        if (parameters.cosm.Rec.include_CF == 1) vname += "Visiblity_Recfast.corr_fac.dat";
        else vname+="Visiblity_Recfast.dat";
        vfile.open(vname.c_str());
        vfile.precision(10);
        
        vfile << "# Visibility function output. Columsn are z, g(z)\n#" << endl;
        
        for (int l=0; l<parameters.nz; l++) vfile << scientific << za[l] << " "
                                                  << cosmos.g_visi(za[l]) << " " << endl;
        
        if (parameters.verbosity>=0) cout << " done. " << endl;
    }
    
    //=======================================================================================================
    // compute Compton y-parameter
    //=======================================================================================================
    double r=0.0, r2=0.0, dt, zeval;
    for(int l=0; l<parameters.nz; l++)
    {
        zeval=(za[l]+za[l-1])*0.5;
        dt=-(za[l]-za[l-1])/cosmos.H(zeval)/(1.0+zeval);
        
        //if(zeval<800.0 && zeval>10.0)
        r+=const_sigT*const_cl*cosmos.Ne(zeval)*const_kb_mec2*cosmos.TCMB(zeval)*(cosmos.Te_Tg(zeval)-1.0)*dt;
        
        if(B0_CosmoRec>0.0)
            r2+=Gamma_decay(zeval, cosmos.H(zeval), B0_CosmoRec, nB_CosmoRec, kD_CosmoRec)
                            /(0.26*const_e*pow(1.0+zeval, 4))*dt;
    }
    
    if(parameters.verbosity>=0)
    {
        cout << " y= " << r << " y_max= " << r2/4.0 << endl;
        cout << " Now exiting program ..." << endl << endl;
    }
    
    delete [] za;
    
    exit(0);        
    
    return;
}


//===========================================================================================================
//
// extend CosmoRec computation at low redshifts using Recfast++
//
//===========================================================================================================
void compute_Recfast_part(double zi, double ze, double Xe_Hi, double Xe_Hei, 
                          double Xei, double rhoi, double dXei,
                          ofstream &ofile, ofstream &ofileRecfast)
{
    if(parameters.verbosity>=0)
        cout << " compute_Recfast_part:: the solution will be extended to z= " << ze 
             << " solving the Recfast system" << endl;
    
    int nz=nz_Recfast_part;
    vector<double> zarr(nz), Xe_H(nz), Xe_He(nz), Xe(nz), TM(nz);
    
    init_xarr(zi, max(1.0e-5, ze), &zarr[0], nz, 1, 0);
    
    double TMi=cosmos.TCMB(zi)*rhoi;
    cosmos.recombine_using_Recfast_system(nz, zi, ze, parameters.cosm.Rec.f_ann,
                                          Xe_Hi, max(1.0e-30, Xe_Hei),
                                          Xei, TMi, dXei, zarr, Xe_H, Xe_He, Xe, TM);
    
    //=======================================================================================================
    // only Xe and Tm
    //=======================================================================================================
    for(int i=1; i<nz; i++)
    {
        write_Xe_Te_to_file(ofile, zarr[i], Xe[i], TM[i]);
        write_X_Recfast_z_to_file(ofileRecfast, zarr[i], cosmos);
    }
    
    //========================================================================================
    // Save results for output to CosmoMC
    //========================================================================================
    if(!Diffusion_correction || (Diffusion_correction && iteration_count==Diff_iteration_max) )
    {
        vector<double> dum(3);
        for(int i=1; i<nz; i++)
        {
            dum[0]=zarr[i]; 
            dum[1]=Xe[i]; 
            dum[2]=TM[i];
        
            output_CosmoRec.push_back(dum);
        }
    }
    
    zarr.clear(); Xe_H.clear(); Xe_He.clear(); Xe.clear(); TM.clear();
    
    return;
}


//===========================================================================================================
//
// for main calculation for CosmoRec
//
//===========================================================================================================
#include "./Modules/main.CosmoRec.cpp"


//===========================================================================================================
//
// CosmoRec part after setup
//
//===========================================================================================================
int CosmoRec()
{
    output_CosmoRec.clear();
    
    string addname_in=parameters.addname;
    
    //==============================================================================
    // setup for the iteration
    //==============================================================================
    bool iteration=(Diffusion_correction ? 1 : 0);
    int err=0;
    
    //==============================================================================
    // set rescale factor for rates etc
    //==============================================================================
    if(FC_CosmoRec.runmode>0) {
        rescale_variables_all(FC_CosmoRec.a_scale, FC_CosmoRec.m_scale, 
                              FC_CosmoRec.runmode, FC_CosmoRec.f_a, FC_CosmoRec.f_b,
                              FC_CosmoRec.f_t, FC_He_switch);
    }

    //==============================================================================
    // enter the main run of ODE --> PDE-solver --> ODE
    //==============================================================================
    do
    {
        if(Diffusion_correction && iteration_count>0)
        { 
            print_message("computing diffusion corrections from iteration # " + int_to_string(iteration_count-1)); 
            
            if(parameters.verbosity>=0) cout << " Entering HI-diffusion part " << endl;
            
            //======================================================================
            // reset parameters
            //======================================================================
            restart_CosmoRec();
                        
            //======================================================================
            // file from previous output
            //======================================================================
            string addnamePesc=".it_"+int_to_string(iteration_count-1)+addname_in;
            if(iteration_count==1) addnamePesc=addname_in;
            
            vector<double> DF_vec_z;
            vector<vector<double> > DF_Ly_n;
            vector<vector<double> > DF_2_gamma_vec;
            vector<vector<double> > DF_Raman_vec;
            vector<vector<double> > DF_nD1s_vec;
            vector<double> DI1_2s_vec;  
            
            //======================================================================
            // solver for the diffusion correction
            // HI populations are communicated; the effective rates have to be
            // preloaded.
            //======================================================================            
            compute_DPesc_with_diffusion_equation_effective(DF_vec_z, DF_Ly_n, 
                                                            DF_2_gamma_vec, DF_Raman_vec, 
                                                            DF_nD1s_vec, DI1_2s_vec, 
                                                            min(Diff_corr_zmax*FC_CosmoRec.f_t*1.25, parameters.zstart/1.002),
                                                            max(Diff_corr_zmin*FC_CosmoRec.f_t, parameters.zend*1.002),
                                                            parameters.CR.nS_2gamma, parameters.CR.nS_Raman, cosmos, 
                                                            Hydrogen_Atoms, 
                                                            pass_on_the_Solution_CosmoRec,
                                                            iteration_count-1);
            
            setup_DF_interpol_data(DF_vec_z, DF_Ly_n, DF_2_gamma_vec,
                                   DF_Raman_vec, DF_nD1s_vec, DI1_2s_vec);
            
            Diffusion_correction_is_on=1;
            
            //======================================================================
            // HeI-diffusion corrections
            // added Feb-May 2011 (Jens Chluba & Jeffrey Fung)
            //======================================================================            
            if(Diffusion_correction_HeI)
            //if(Diffusion_correction_HeI && iteration_count==Diff_iteration_max)
            {
                if(parameters.verbosity>=0) cout << "\n Entering HeI-diffusion part " << endl;

                vector<int> HeI_level_indices;

                if(DI1_2s_correction_on_HeI) switch_on_1s_2s_HeI_correction();
                else switch_off_1s_2s_HeI_correction();
                
                bool last_it=0;
                if(iteration_count==Diff_iteration_max) last_it=1;
                
                compute_HeI_DPesc_with_diffusion_equation_effective(HeI_Atoms.Get_nShells(), DF_vec_z,
                                                                    DF_Ly_n, DI1_2s_vec,
                                                                    HeI_level_indices,
                                                                    min(3200.0*FC_CosmoRec.f_t, parameters.zstart/1.002),
                                                                    HeI_was_switched_off_at_z*1.002,
                                                                    cosmos, 
                                                                    HeI_Atoms,
                                                                    pass_on_the_Solution_CosmoRec_HeI,
                                                                    iteration_count-1, last_it);
                
                setup_DF_interpol_data_HeI(DF_vec_z, DF_Ly_n, DI1_2s_vec, HeI_level_indices);
                
                Diffusion_correction_HeI_is_on=1;
            }
            
            //======================================================================
            // change the name for output
            //======================================================================
            if(iteration_count<Diff_iteration_max) 
                parameters.addname=".it_"+int_to_string(iteration_count)+addname_in;
            else parameters.addname=".trans"+addname_in;
            
            print_message("starting iteration # " + int_to_string(iteration_count)); 
        }
        
        //==========================================================================
        // computation of recombination history with EMLA
        //==========================================================================
        err=compute_history_with_effective_rates();  // defined in main.CosmoRec.cpp
        
        if(err!=0) iteration=0;
            
        iteration_count++;      
        
        if(iteration_count>Diff_iteration_max) iteration=0;
    }
    while(iteration);

    //==============================================================================
    // restore initial state before run (important for batch mode run)
    //==============================================================================
    iteration_count=Diff_iteration_min;
    flag_He=1;              // has to be reset, since after each run it will be == 0
    parameters.addname=addname_in;
    Diffusion_correction_is_on=0;
    Diffusion_correction_HeI_is_on=0;
    
    return err;
}


//===========================================================================================================
//
// CosmoRec for run from console giving a parameter file
//
//===========================================================================================================
int CosmoRec(int narg, char *args[])
{
    string filename=args[narg-1], mode="";
    
    if(narg==2){}
    else if(narg==3)
    {
        mode=args[narg-2];
        if(mode=="REC") print_message(" Running Recfast++ module.");
        else{ mode=""; print_message(" Runmode does not exist. Exiting"); exit(0); }
    }
    else{ cerr << " Too many/few parameters " << endl; exit(1); }
    
    //========================================================================
    // just run history using the recfast module in the cosmology class
    //========================================================================
    if(mode=="REC") run_Recfast_pp(filename);

    //========================================================================
    // initialize cosmology and the general startup
    //========================================================================
    print_message(" You are entering CosmoRec "+CosmoRec_version+".");
    read_entries_from_parameter_file(filename, parameters);
    set_startup_data_CR(filename, parameters, cosmos);

    set_array_dimensions_and_parameters();
    allocate_memory();
    
    //========================================================================
    // initialise zarr (linear)
    //========================================================================
    init_xarr(parameters.zstart, parameters.zend, zarr, parameters.nz, 0, 0);
    
    //========================================================================
    // enter CosmoRec
    //========================================================================
    int err=CosmoRec();
    if(err!=0){ cerr << " CosmoRec:: run was not completed. Exiting. " << endl; }
    
    //========================================================================
    // clean up memory
    //========================================================================
    cleanup_CosmoRec();

    if(parameters.verbosity>=0) cout << " CosmoRec finished normally. " << endl;

    return err;
}



//===========================================================================================================
//
// to call CosmoRec in batch mode, i.e. avoiding initialization & loading tables
//
//===========================================================================================================
static bool CosmoRec_is_initialized=0;

//===========================================================================================================
int return_solution_mem_Xe;
int return_solution_mem_Te;

void return_solution_to_calling_program(const int nz, double *z_arr, double *Xe_arr, double *Te_arr)
{
    int nz_loc=output_CosmoRec.size();
    vector<double> za(nz_loc);
    vector<double> ya(nz_loc);
    
    //=======================================================================================================
    // Xe
    //=======================================================================================================
    for(int i=0; i<nz_loc; i++)
    { 
        za[nz_loc-1-i]=output_CosmoRec[i][0]; 
        ya[nz_loc-1-i]=log(output_CosmoRec[i][1]); 
    }

    double zmin=min(za[0], za[nz_loc-1]);
    double zmax=max(za[0], za[nz_loc-1]);
    
    if(!CosmoRec_is_initialized)
        return_solution_mem_Xe=calc_spline_coeffies_JC(nz_loc, &za[0], &ya[0],
                                                       "return_solution_to_calling_program :: Xe");      
    else update_spline_coeffies_JC(return_solution_mem_Xe, nz_loc, &za[0], &ya[0]);      
        
    //=======================================================================================================
    // Te
    //=======================================================================================================
    for(int i=0; i<nz_loc; i++)
        ya[nz_loc-1-i]=log(output_CosmoRec[i][2]); 
    
    if(!CosmoRec_is_initialized)
        return_solution_mem_Te=calc_spline_coeffies_JC(nz_loc, &za[0], &ya[0],
                                                       "return_solution_to_calling_program :: Te");      
    else update_spline_coeffies_JC(return_solution_mem_Te, nz_loc, &za[0], &ya[0]);      
    
    //=======================================================================================================
    // copy output
    //=======================================================================================================
    for(int i=0; i<nz; i++)
    { 
        if(z_arr[i]>=zmax || z_arr[i]<=zmin) 
        {
            Xe_arr[i]=cosmos.Xe_Seager(z_arr[i]);
            Te_arr[i]=cosmos.Te(z_arr[i]);
        }
        else 
        {
            Xe_arr[i]=exp(calc_spline_JC(z_arr[i], return_solution_mem_Xe));
            Te_arr[i]=exp(calc_spline_JC(z_arr[i], return_solution_mem_Te));                
        }
    }
    
    return;
}

//===========================================================================================================
// Just return Recfast++ output
//===========================================================================================================
void run_Recfast_pp(const int nz, double *z_arr,
                    double *Xe_arr, double *Te_arr,
                    double mu1, double mu2, double mu3, double mu4, int CT_corr=0)
{
    check_and_set_params(parameters, parameters.verbosity);
        
    if(parameters.verbosity>=1) display_Recfast_pp_parameters("", parameters);
    
    parameters.cosm.Rec.include_CF=CT_corr;  // overwrite just this one setting if requested.
    set_cosmology(parameters.cosm, cosmos);
    
    for(int l=0; l<nz; l++)
    {
        Xe_arr[l]=cosmos.Xe_Seager(z_arr[l]);
        Te_arr[l]=cosmos.Te(z_arr[l]);
    }
    
    //=======================================================================================================
    // Add PCA terms if mode amplitudes are non-zero
    //=======================================================================================================
    add_PCA_functions(mu1, mu2, mu3, mu4, nz, z_arr, Xe_arr);

    //=======================================================================================================
    // write outputs of Xe and Te
    //=======================================================================================================
    if(write_Xe_Te_sol>0)
    {
        if(parameters.verbosity>=0) cout << " writing data to " << parameters.path << endl;
    
        ofstream ofile;
        string name=COSMORECDIR+"./outputs/";
        if(CT_corr==1) name+="Xe_Recfast.corr_fac"+parameters.addname+".dat";
        else name+="Xe_Recfast"+parameters.addname+".dat";
        ofile.open(name.c_str());
        ofile.precision(10);
        
        ofile << "# Recfast++ output. The columns are z, Xe=Ne/NH, Te\n#" << endl;
        
        for(int l=0; l<nz; l++) ofile << scientific << z_arr[l] << " "
                                      << Xe_arr[l] << " "
                                      << Te_arr[l] << endl;
        
        if(parameters.verbosity>=0)
            cout << " finished Recfast++ run. Now exiting program ... " << endl << endl;
    }
    
    return;
}

//===========================================================================================================
// to call CosmoRec with a list of parameters in batch mode (e.g. when calling it from CosmoMC).
//
//----------- Runmodes --------------------------------------------------------------------------------------
//
// runmode == 0: CosmoRec run with diffusion
//            1: CosmoRec run without diffusion
//            2: Recfast++ run (equivalent of the original Recfast version)
//            3: Recfast++ run with correction function of Chluba & Thomas, 2010
//
// On entry, the array z_arr should contain the redshifts at which Xe and Te are required. nz
// determines the number of redshift points. Xe_arr and Te_arr will contain the solution on exit.
//
// Furthermore, runpars[0] defines the dark matter annihilation efficiency in eV/s.
// runpars[1] switches the accuracy of the recombination model:
//
//----------- Output parameters -----------------------------------------------------------------------------
//
// runpars[1]==-1: closest equivalent of 'HyRec' case (Haimoud & Hirata, 2010)
// runpars[1]== 0: default setting
// runpars[1]== 1: 2g for n<=4 & Raman for n<=3
// runpars[1]== 2: 2g for n<=8 & Raman for n<=7
// runpars[1]== 3: 2g for n<=8 & Raman for n<=7 + Helium feedback up to n=5
// runpars[1]== 4: default setting              + Helium radiative transfer
// runpars[1]== 5: 2g for n<=4 & Raman for n<=3 + Helium radiative transfer up to n=3
// runpars[1]== 6: 2g for n<=4 & Raman for n<=3 + Helium radiative transfer up to n=5 (full setting)
//
// The value of runpars[1] is only important for runmode 0 & 1.
//
// CosmoRec can be set to output the ionization history & electron temperature, solutions for
// the populations of the resolved levels, and some information about the Cosmology (usually
// directed into the directory "./outputs/"). This can cause problems when using CosmoRec
// with a parallel CosmoMC, however, it possible to control the outputs through runpars[2]:
//
//----------- Output parameters -----------------------------------------------------------------------------
//
// runpars[2]==0: don't write out anything (default).
// runpars[2]==1: write out only the recombination history.
// runpars[2]==2: write out the recombination history, and the cosmology.
// runpars[2]==3: write out the recombination history, populations, and the cosmology.
//
//----------- A2s1s two-photon rate -------------------------------------------------------------------------
//
// runpars[3]>=0: use given number to set const_HI_A2s_1s. For ==0 default value is used.
//                Do not vary too much [A2s1s ~ (6 - 10) 1/sec]!
//
//-----------------------------------------------------------------------------------------------------------
// heating by magnetic fields mainly following Kunze & Komatsu 2013
//
// runpars[4]> 0: B0 in nG. For ==0 no heating is added.
// runpars[5]   : nB spectral index for magnetic fields
// runpars[6]  0: turbulent decay; 1: ambipolar diffusion; 2: both
//            <0: same but using Paoletti & Finelli norm for B0
//
//-----------------------------------------------------------------------------------------------------------
// runpars[7]!=1: Fine structure constant rescaled
// runpars[8]!=1: Mass of electron rescaled
// runpars[9]!=0: Redshift dependent power law of fundamental constants scaling
// runpars[10] 0: no rescaling;    1: Boltzmann/T rescaling; 2: sigma_T rescaling
//             3: A2s1s rescaling; 4: Phot rates rescaled; 5: Pesc rescaled (only for Recfast++ runmodes)
//             6: All components rescaled
//
//-----------------------------------------------------------------------------------------------------------
// PCA terms EM1-EM4 from Farhang et al., 2012 & 2013 [only when calling in batch mode]
//
// runpars[11-14]: EM1..EM4
//
//-----------------------------------------------------------------------------------------------------------
// CAMB overrides. Negative values keep the default or selected accuracy preset.
//
// runpars[15]>=0: override Diff_iteration_max.
// runpars[16]>=0: override nShells.
// runpars[17]>=0: override nShellsHeI.
// runpars[18]>=0: override HI_absorption.
// runpars[19]>=0: override nS_2gamma.
// runpars[20]>=0: override nS_Raman.
//
//===========================================================================================================
int CosmoRec(const int runmode, const double runpars[21],
             const double omegac, const double omegab, 
             const double omegak, const double Nnu,  
             const double h0, const double tcmb, const double yhe,
             const int nz, double *z_arr, double *Xe_arr, double *Te_arr, 
             const int label)
{
    //=======================================================================================================
    // switch off messages 
    //=======================================================================================================
    parameters.verbosity=-1;
    parameters.path=COSMORECDIR+"./outputs/";
    parameters.addname=".batch."+int_to_string(label, 6)+".dat";
    
    //=======================================================================================================
    // switch on/off outputs
    //=======================================================================================================
    write_Xe_Te_sol=0;
    write_populations=0;
    bool write_cosmology=0;
    
    if((int)runpars[2] > 0) write_Xe_Te_sol=1;
    if((int)runpars[2] > 1) write_cosmology=1;
    if((int)runpars[2] > 2) write_populations=1;

    //=======================================================================================================
    // general setup
    //=======================================================================================================
    parameters.nz=3000;     // number of redshift points for CosmoRec computation;
                            // for higher accuracy this number should be increased;
                            // a smaller number than ~750 is not recommended.
    parameters.zstart=3000; // starting redshift
    parameters.zend  =0;    // ending redshift
    
    //=======================================================================================================
    // cosmological parameters
    //=======================================================================================================
    parameters.cosm.Yp= yhe;                   // Helium mass fraction
    parameters.cosm.T0= tcmb;                  // CMB temperature today
    //
    parameters.cosm.Omega_m= omegac + omegab;  // Omega_matter
    parameters.cosm.Omega_b= omegab;           // Omega_b
    parameters.cosm.Omega_L= 0.0;              // Omega_Lambda (computed internally)
    parameters.cosm.Omega_k= omegak;           // Omega_curv
    //
    parameters.cosm.h100 = h0;                 // h100
    parameters.cosm.Nnu  = Nnu;                // N_nu_eff
    
    //=======================================================================================================
    // Recfast++ settings [some are default values]
    //=======================================================================================================
    parameters.cosm.Rec.F= 1.14;               // Recfast++ fudge-factor
    parameters.cosm.Rec.A2s1s=runpars[3];      // explicitly set value of A2s1s rate (added July 23.07.2014)
    parameters.cosm.Rec.f_ann=runpars[0];
    
    //=======================================================================================================
    // magnetic field heating (added 11.10.2014)
    //=======================================================================================================
    int PMFmode=runpars[6];
    bool turb=0, ambi=0, KKLor=0;
    if(PMFmode<0){ KKLor=1; PMFmode*=-1.0; }
    else if(PMFmode==0) turb=1;
    else if(PMFmode==1) ambi=1;
    else if(PMFmode==2){ turb=1; ambi=1; }
    else { cerr << " PMF mode not existing... " << endl; exit(1); }
    
    parameters.cosm.Rec.B0=runpars[4];
    parameters.cosm.Rec.nB=runpars[5];
    parameters.cosm.Rec.turb=turb;
    parameters.cosm.Rec.ambi=ambi;
    parameters.cosm.Rec.KK_or_PF_Lorentz=KKLor;
    
    //=======================================================================================================
    // scaling of fundamental constants using input variables from camb
    //=======================================================================================================
    if(runpars[10] != 0 && runpars[9] != 0.0 && runmode < 2)
        throw_error("CosmoRec", "Cannot currently run CosmoRec with the power law!");
    
    parameters.cosm.Rec.aS=runpars[7];
    parameters.cosm.Rec.mS=runpars[8];
    parameters.cosm.Rec.pS=runpars[9];
    parameters.cosm.Rec.RM=runpars[10];

    //=======================================================================================================
    // runmode 2: return Recfast++ result
    //=======================================================================================================
    if(runmode==2)
    {
        run_Recfast_pp(nz, z_arr, Xe_arr, Te_arr, runpars[11], runpars[12], runpars[13], runpars[14], 0);
        return 0;
    }

    //=======================================================================================================
    // runmode 3: return Recfast++ result with CT2010 correction function
    //=======================================================================================================
    if(runmode==3)
    {
        run_Recfast_pp(nz, z_arr, Xe_arr, Te_arr, runpars[11], runpars[12], runpars[13], runpars[14], 1);
        return 0;
    }
    
    //=======================================================================================================
    // parameters for the runmodes & output
    //=======================================================================================================
    parameters.CR.nShells=3;
    parameters.CR.nS_effective=500;
    
    parameters.CR.nShellsHeI=2;
    parameters.CR.HI_absorption=2;
    parameters.CR.spin_forbidden=1;
    parameters.CR.HeI_Feedback=0;
    
    parameters.CR.Diffusion_flag=1;
    parameters.CR.induced_flag=2;
    parameters.CR.nS_2gamma=3;
    parameters.CR.nS_Raman=2;
    
    //=======================================================================================================
    // number of iterations over HI & HeI diffusion problem
    //=======================================================================================================
    Diff_iteration_max=2;   // Match CosmoRec's internal default.
                            // CAMB can explicitly override this through
                            // runpars[15] / cosmorec_diff_iteration_max.
    if((int)runpars[15] >= 0) Diff_iteration_max=(int)runpars[15];
    
    //=======================================================================================================
    // settings for non-default cases
    //=======================================================================================================
    if((int)runpars[1]==-1)
    {
        parameters.CR.nShells=4;
        parameters.CR.nS_2gamma=4; parameters.CR.nS_Raman=3;
        parameters.CR.nShellsHeI=2;
        parameters.CR.HeI_Feedback=-2; parameters.CR.HI_absorption=1;
    }
    else if((int)runpars[1]==1)
    {
        parameters.CR.nShells=4;
        parameters.CR.nS_2gamma=4; parameters.CR.nS_Raman=3;
    }
    else if((int)runpars[1]==2)
    {
        parameters.CR.nShells=10;
        parameters.CR.nS_2gamma=8; parameters.CR.nS_Raman=7;
        parameters.CR.nShellsHeI=3;
        parameters.CR.HeI_Feedback=0;
    }
    else if((int)runpars[1]==3)
    {
        parameters.CR.nShells=10;
        parameters.CR.nS_2gamma=8; parameters.CR.nS_Raman=7;
        parameters.CR.nShellsHeI=5;
        parameters.CR.HeI_Feedback=-5;
    }
    
    //=======================================================================================================
    // runmode with helium diffusion module
    //=======================================================================================================
    else if((int)runpars[1]==4){ parameters.CR.HeI_Feedback=-2; parameters.CR.HI_absorption=3; }
    else if((int)runpars[1]==5)
    {
        parameters.CR.nShells=4;
        parameters.CR.nS_2gamma=4; parameters.CR.nS_Raman=3;
        parameters.CR.nShellsHeI=3;
        parameters.CR.HeI_Feedback=-3; parameters.CR.HI_absorption=3;
    }
    else if((int)runpars[1]==6)
    {
        parameters.CR.nShells=4;
        parameters.CR.nS_2gamma=4; parameters.CR.nS_Raman=3;
        parameters.CR.nShellsHeI=5;
        parameters.CR.HeI_Feedback=-5; parameters.CR.HI_absorption=3;
    }
    
    if((int)runpars[16] >= 0) parameters.CR.nShells=(int)runpars[16];
    if((int)runpars[17] >= 0) parameters.CR.nShellsHeI=(int)runpars[17];
    if((int)runpars[18] >= 0) parameters.CR.HI_absorption=(int)runpars[18];
    if((int)runpars[19] >= 0) parameters.CR.nS_2gamma=(int)runpars[19];
    if((int)runpars[20] >= 0) parameters.CR.nS_Raman=(int)runpars[20];
    
    //=======================================================================================================
    if(runmode==1) parameters.CR.Diffusion_flag=0;

    //=======================================================================================================
    // initialize cosmology and the general startup
    //=======================================================================================================
    if(!CosmoRec_is_initialized)
    {
        cout << "\n ====================================================="
             << "===================================================== " << endl;
        cout << " || You are entering CosmoRec "+CosmoRec_version+" (batch mode)." 
             << " Several initializations will be executed now." << endl;
        cout << " ====================================================="
             << "===================================================== " << endl;
    }
    
    int err=0;
    
    do
    {
        set_startup_data_CR("", parameters, cosmos);
        set_array_dimensions_and_parameters();  
        
        if(!CosmoRec_is_initialized && err==0) allocate_memory();
        
        //===================================================================================================
        // write info about cosmology
        //===================================================================================================
        if(write_cosmology==1)
            cosmos.dump_cosmology(COSMORECDIR+"./outputs/batch.cosmology."+int_to_string(label, 6)+".dat");
        
        //===================================================================================================
        // initialise zarr (only has to be done once...)
        //===================================================================================================
        init_xarr(parameters.zstart, parameters.zend, zarr, parameters.nz, 0, 0);
        
        //===================================================================================================
        // enter CosmoRec
        //===================================================================================================
        err=CosmoRec();
        
        if(err!=0)
        { 
            cout << " CosmoRec:: run was not completed. Restarting Solver. " << err << endl; 
            wait_f_r(); 
        }
    }   
    while(err!=0);

    //=======================================================================================================
    // copy solution to output arrays
    //=======================================================================================================
    return_solution_to_calling_program(nz, z_arr, Xe_arr, Te_arr);

    //=======================================================================================================
    // Add PCA terms if mode amplitudes are non-zero
    //=======================================================================================================
    add_PCA_functions(runpars[11], runpars[12], runpars[13], runpars[14], nz, z_arr, Xe_arr);
    
    //=======================================================================================================
    // exiting CosmoRec
    //=======================================================================================================
    if(parameters.verbosity>=0) cout << " finished run... " << endl;
    
    if(!CosmoRec_is_initialized)
    {
        cout << "\n ====================================================="
             << "===================================================== " << endl;
        cout << " || The first run was successful." 
             << " Will now switch to 'silent' mode." << endl;
        cout << " ====================================================="
             << "===================================================== " << endl;
    }

    CosmoRec_is_initialized=1;

    return err;
}

//===========================================================================================================
// Wrap the C++ Fortran routine to be allow calling from Fortran. Arguments are as above.
// Added 06.03.2011 (Richard Shaw)
//===========================================================================================================
// 18.06.2012: Added version to communicate Hubble pointer
// 08.06.2012: Added version to communicate Hubble 
//===========================================================================================================
extern "C" {
    
    void cosmorec_calc_cpp_(const int * runmode, const double * runpars, 
                            const double * omega_c, const double * omega_b, const double * omega_k, 
                            const double * num_nu, const double * h0, 
                            const double * t_cmb, const double * y_he, 
                            double * za_in, double * xe_out, double * tb_out, 
                            const int * len, const int* label) 
    {
        //======================================================================
        // Call standard CosmoRec routine
        //======================================================================
        CosmoRec(*runmode, runpars, *omega_c, *omega_b, *omega_k, *num_nu, 
                 *h0 / 100.0, *t_cmb, *y_he, *len, za_in, xe_out, tb_out, *label);
        
        return;
    } 
    
    void cosmorec_calc_h_cpp_(const int * runmode, const double * runpars, 
                              const double * omega_c, const double * omega_b, const double * omega_k, 
                              const double * num_nu, const double * h0, 
                              const double * t_cmb, const double * y_he, 
                              const double * z_Hz, const double * Hz, const int * nz,
                              double * za_in, double * xe_out, double * tb_out, 
                              const int * len, const int* label) 
    {
        cosmos.init_Hubble(z_Hz, Hz, *nz);
        
        //======================================================================
        // Call standard CosmoRec routine
        //======================================================================
        CosmoRec(*runmode, runpars, *omega_c, *omega_b, *omega_k, *num_nu, 
                 *h0 / 100.0, *t_cmb, *y_he, *len, za_in, xe_out, tb_out, *label);
        
        cosmos.clear_Hubble();
        
        return;
    } 
    
    void cosmorec_calc_hptr_cpp_(const int * runmode, const double * runpars, 
                                 const double * omega_c, const double * omega_b, const double * omega_k, 
                                 const double * num_nu, const double * h0, 
                                 const double * t_cmb, const double * y_he, 
                                 double (* Hptr)(const double * z), 
                                 double * za_in, double * xe_out, double * tb_out, 
                                 const int * len, const int* label) 
    {
        cosmos.init_Hubble(Hptr);
        
        //======================================================================
        // Call standard CosmoRec routine
        //======================================================================
        CosmoRec(*runmode, runpars, *omega_c, *omega_b, *omega_k, *num_nu, 
                 *h0 / 100.0, *t_cmb, *y_he, *len, za_in, xe_out, tb_out, *label);
        
        cosmos.clear_Hubble();
        
        return;
    } 
    
    //=======================================================================================================
    //
    // electron temperature equation (rho=Te/Tg) [added by Vivian Poulin, summer 2017]
    //
    //=======================================================================================================
    void evaluate_TM(double z, double Xe, double fHe, double rho, double Tg, double Hz, double * drho_dt)
    {
        ODE_effective::evaluate_TM(z, Xe, fHe, rho, Tg, Hz, *drho_dt);
        return;
    }
}

//===========================================================================================================
