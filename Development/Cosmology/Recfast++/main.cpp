//========================================================================================
// Authors: Jens Chluba & Luke Hart
// First implementation: Oct 2010
// CITA, University of Toronto
// All rights reserved.
//========================================================================================
// 05.12.2016: added simple redshift scaling of constants ~(1+z)^p
// 01.12.2016: changed setup of Tr flag and added setup for variation of fundamental
//             constants (developed by Luke Hart)
// 14.11.2014: added setup for magnetic field heating modes

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "recombination.Recfast.h"
#include "cosmology.Recfast.h"
#include "constants.Recfast.h"
#include "evalode.Recfast.h"
#include "Magnetic_fields.Recfast.h"
#include "Variation_constants.Recfast.h"

using namespace std;

bool verbose=1;

//========================================================================================
// define default HI A2s1s two-photon decay rate. Added to allow changing it globally.
// (added 23.07.2014 JC)
//========================================================================================
double RECFAST_atomic_constants::RF_Lam2s1sH = 8.22458;         // Recfast++ default value

//========================================================================================
void show_cosmological_parameters(vector<double> &param)
{
    cout << "\n Recfast++ :: show_cosmological_parameters: " << endl;
    cout << " running Recfast++ history for the following parameters: " << endl;   
    
    cout << "\n zs: " << param[1] << "\t ze: " << param[2] << "\t nzpts: " << (int)param[0] << endl;
    cout << " Y_p: " << param[3] << "\t TCMB: " << param[4] << endl;
    cout << " OmegaM: " << param[5] << "\t OmegaB: " << param[6] 
         << " OmegaL: " << param[7] << "\t OmegaK: " << param[8] << endl;
    cout << " Omega_rel (photons & neutrinos): " << calc_Orel(param[4], param[10], param[9]) 
         << " Nnu= " << param[10] << endl;
    
    cout << " Hubble constant in units of 100 km s-1 Mpc-1: " << param[9] << endl;
    cout << " Fudge factor for H recombination: " << param[11] << endl;
    
    cout << " Using HI A2s1s = " << RECFAST_atomic_constants::RF_Lam2s1sH << " s^-1";
    if(param[14]==0.0) cout << " (default value)";
    cout << endl;

    if(param[12]!=0.0)
        cout << " DM-annihilations with efficiency " << param[12] << " will be included " << endl;
    if(param[13]==1) 
        cout << " Will include recombination corrections following Chluba & Thomas 2010 / Hart & Chluba 2017" << endl;
    
    if(Get_T_eval()) cout << " Using T=Tr for photo-ionization rate" << endl;
    else cout << " Using T=Te for photo-ionization rate (Recfast default)" << endl;

    if(Get_B0_PMF()>0.0) cout << " Will include heating by magnetic fields "
                              << "(decay==" << Get_decay_setting()
                              << ", ambi==" << Get_ambi_setting()
                              << ", Lorentz==" << Get_Lorentz_setting() << ")" << endl;

    if(Get_aScale()!=1.0 || Get_mScale()!=1.0  || Get_pScale()!=0.0)
        cout << " Variation of fundamental constants: alp/alp_r== " << Get_aScale()
             << " me/me_r== " << Get_mScale()
             << " pow== " << Get_pScale()
             << " mode== " << Get_rMode() << endl;

    cout << endl;
    
    return;
}

//========================================================================================
void read_parameter_file(string fname, vector<double> &params)
{
    if(verbose) 
        cout << " Recfast++ :: read_parameter_file : reading parameters from: " 
             << fname << endl;
    
    ifstream ifile(fname.c_str());
    if(!ifile)
    { 
        cerr << " Recfast++ :: read_parameter_file :" 
             << " error reading parameter-file. Exiting. " << endl; 
        exit(0); 
    }
    
    //====================================================================================
    ifile >> params[3];  // Yp
    ifile >> params[4];  // Temperature of CMB at z=0
    ifile >> params[5];  // Omega matter 
    ifile >> params[6];  // Omega Baryons 
    ifile >> params[7];  // Omega Lambda
    ifile >> params[8];  // Omega Curvature 
    ifile >> params[9];  // h100
    ifile >> params[10]; // effective number of neutrinos 
    ifile >> params[11]; // fudge-factor; normally F=1.14
    
    ifile >> params[12]; // fDM [eV/s] which gives annihilation efficiency;
                         // typical value fDM=2.0e-24 eV/s (see Chluba 2010 for definitions)
    ifile >> params[13]; // switch on/off recombination corrections (Chluba & Thomas 2010)

    ifile >> params[14]; // HI A2s1s
    if(params[14]>0.0) RECFAST_atomic_constants::RF_Lam2s1sH=params[14];

    //====================================================================================
    // magnetic field heating
    //====================================================================================
    double B0, nB; bool turb, ambi;
    ifile >> B0;   // B0 in nG
    ifile >> nB;   // nB
    ifile >> turb; // decaying magnetic turb on / off
    ifile >> ambi; // ambipolar diff on / off
    
    if(B0>0.0 && (turb || ambi) )
    {
        double kD=286.91/B0;         // in Mpc^-1; cos(theta)=1
        Set_PMF(turb, ambi, 0, B0, nB, kD);
        Set_Tr_eval();
    }
    else Switch_PMF_off();           // switch PMF heating off
    
    //====================================================================================
    // variation of alpha and me/mp
    //====================================================================================
    double aS, mS, pS; int RM;
    ifile >> aS; // alp/alp_r
    ifile >> mS; // me/me_r
    ifile >> pS; // power for (1+z)^p
    ifile >> RM; // mode
    
    if(RM>0) Set_VFC_params(aS, mS, pS, RM);
    else Switch_VFC_off(); // switch VFC off
    
    //====================================================================================
    // compute Omega_L according to Omega_K
    //====================================================================================
    if(params[7]<=0.0) 
        params[7]=1.0-params[5]-params[8]-calc_Orel(params[4], params[10], params[9]);
    
    ifile.close();
    
    if(verbose) show_cosmological_parameters(params);
    
    return;
}

//========================================================================================
void write_Recfast_output(string oname, vector<double> &zarr, vector<double> &Xe_H, 
                          vector<double> &Xe_He, vector<double> &Xe, vector<double> &TM)
{
    if(verbose) 
        cout << " Recfast++ :: write_Recfast_output : writing output into: " 
             << oname << endl;

    ofstream ofile(oname.c_str());
    ofile.precision(8);
    
    for(int i=0; i<(int)zarr.size(); i++)
        ofile << zarr[i] << " " << Xe[i] << " " << Xe_H[i] << " " 
              << Xe_He[i] << " " << TM[i] << endl;
    
    return;
}

//========================================================================================
// main program
//========================================================================================
int main(int narg, char *args[])
{
    //====================================================================================
    // System call to clear screen 
    //====================================================================================
    system("clear");
    
    if(narg < 2) 
    {
        cout << " Recfast-error: usage $ ./Recfast++ parameters.dat" << endl;
        exit(0);
    }
    
    string fname=args[narg-1];

    //====================================================================================
    // reading parameters
    //====================================================================================
    vector<double> params(15);
    int npz=10000;
    double zstart=2.5e+4, zend=0.001;
    
    params[0]=npz;
    params[1]=zstart;
    params[2]=zend;
    
    read_parameter_file(fname, params);
    
    //====================================================================================
    // running Recfast part
    //====================================================================================
    if(verbose) cout << " Recfast++ :: Entering computation. " << endl;

    vector<double> zarr(npz), Xe_H(npz), Xe_He(npz), Xe(npz), TM(npz);
    
    //====================================================================================
    // main code
    //====================================================================================
    Xe_frac(&params[0], &zarr[0], &Xe_H[0], &Xe_He[0], &Xe[0], &TM[0]);
    
    string oname="./outputs/Xe_Recfast";
    if(params[12]!=0.0) oname+=".DM_annihilations";
    if(params[13]==1) oname+=".Rec_corrs_CT2010";
    if(Get_B0_PMF()>0.0)
    {
        oname+=".PMF_heating";
        if(Get_decay_setting()) oname+=".turb";
        if(Get_ambi_setting()) oname+=".ambi";
    }

    if((Get_aScale()!=1.0 || Get_pScale()!=0.0) && Get_aScale()>0.0) oname+=".alp";
    if((Get_mScale()!=1.0 || Get_pScale()!=0.0) && Get_mScale()>0.0) oname+=".me";
    if(Get_pScale()!=0.0) oname+=".power";
    if(Get_rMode()==1) oname+=".Boltz";
    else if(Get_rMode()==2) oname+=".sigT";
    else if(Get_rMode()==3) oname+=".A2s1s";
    else if(Get_rMode()==4) oname+=".phot";
    else if(Get_rMode()==5) oname+=".Ly-a";
    else if(Get_rMode()==6) oname+=".all";
    if(Get_T_eval()) oname+=".Tr";
    oname+=".dat";
    
    write_Recfast_output(oname, zarr, Xe_H, Xe_He, Xe, TM);

    if(verbose) cout << " Recfast++ :: Run completed. Exiting. " << endl;

    return 0;
}
