//========================================================================================
// Author: Jens Chluba & Luke Hart
// Last modification: Jan 2017
// CITA, University of Toronto
// All rights reserved.
//========================================================================================
// 20.01.2017: added switch for selection of correction function (LH & JC)
// 04.01.2017: added new routines for separate H and HeI correction functions (LH)

//========================================================================================
// This module allows to add the recombination corrections as described by 
// Chluba & Thomas 2010. Simple linear interpolation of DNe/Ne is used, as explained by
// Rubino-Martin et al. 2010.
//========================================================================================

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

#include "constants.Recfast.h"
#include "Rec_corrs_CT.Recfast.h"

using namespace std;
using namespace RECFAST_physical_constants; // defined in constants.Recfast.h

//========================================================================================
// data-file with correction
//========================================================================================
// if separate correction functions requested --> new file & scheme
bool file_sect=(RECFAST_atomic_constants::RF_separate_corr ? 1 : 0); // new file for fH & fHe but old for ftot
//bool file_sect=(RECFAST_atomic_constants::RF_separate_corr ? 1 : 1); // using new file for both

// if Tr evaluation wanted --> new file with respective evaluation
bool use_Tr=(RECFAST_atomic_constants::RF_eval_ion_Tr ? 1 : 0);

string Rec_corrs_CT_file_old=(string)RECFASTPPPATH+"./src/Data/DXe_Xe.CT2010.dat";
string Rec_corrs_CT_file_Tr=(string)RECFASTPPPATH+"./src/Data/Corr_LH_Tr.dat";
string Rec_corrs_CT_file_noTr=(string)RECFASTPPPATH+"./src/Data/Corr_LH.dat";

string Rec_corrs_CT_file_new=(use_Tr ? Rec_corrs_CT_file_Tr : Rec_corrs_CT_file_noTr);
string Rec_corrs_CT_file=(file_sect ? Rec_corrs_CT_file_new : Rec_corrs_CT_file_old);

//========================================================================================
// some global variables for compution
//========================================================================================
int Rec_corrs_CT_start_index=0;
vector<vector<double> > Rec_corrs_CT_Data;

//========================================================================================
// reading the correction factor
//========================================================================================
void read_recombination_correction(int mflag)
{
    Rec_corrs_CT_start_index=0;
    if(Rec_corrs_CT_Data.size()>0) return;
    
    ifstream ifile(Rec_corrs_CT_file.c_str());
    if(!ifile)
    { 
        cerr << " Recfast++ :: read_recombination_correction:" 
             << " file with correction factor ( " << Rec_corrs_CT_file 
             << " ) not found. Exiting. " << endl; 
        exit(0);
    }
    
    int nCol=(file_sect ? 4 : 2);
    vector<double> dum(nCol);
    Rec_corrs_CT_Data.clear();
    
    if(mflag>=1) 
    {
        cout << " Recfast++ :: read_recombination_correction:";
    
        if(!file_sect)
            cout << " reading recombination correction according to Chluba & Thomas 2010 for ftot " << endl;
        else
        {
            if(RECFAST_atomic_constants::RF_separate_corr)
                cout << " reading recombination correction according to Hart & Chluba 2017 for fH & fHe " << endl;
            else cout << " reading recombination correction according to Hart & Chluba 2017 for ftot " << endl;
        }
    }
    
    do
    {
        ifile >> dum[0];
        ifile >> dum[1];
        if(file_sect){ ifile >> dum[2]; ifile >> dum[3]; }
        
        if(ifile) Rec_corrs_CT_Data.push_back(dum);
    }
    while(ifile);

    return;
}

//========================================================================================
double interpolate_data(double z,int cflag)
{
    //====================================================================================
    // find next k index with z>z(k)
    //====================================================================================
    if(z<=Rec_corrs_CT_Data[Rec_corrs_CT_start_index][0]) Rec_corrs_CT_start_index++;
    
    //====================================================================================
    // simple linear interpolation [added LH]
    //====================================================================================
    // since we have three different sets to interpolate now 
    // 1 -- basic correction f_H+f_He
    // 2 -- hydrogen correction f_H
    // 3 -- helium correction f_He
    // Add an error flag in case the values aren't the ones we want
    //====================================================================================
    if (cflag != 1 && cflag != 2 && cflag != 3) {
        cerr << " Recfast++ :: interpolate_data:"
             << " invalid selection of correction function flag. exitting. " << endl;
        exit(0);
    }
    
    double dz=z-Rec_corrs_CT_Data[Rec_corrs_CT_start_index-1][0];
    double df_dz=(Rec_corrs_CT_Data[Rec_corrs_CT_start_index][cflag]
                 -Rec_corrs_CT_Data[Rec_corrs_CT_start_index-1][cflag])
                /(Rec_corrs_CT_Data[Rec_corrs_CT_start_index][0]
                 -Rec_corrs_CT_Data[Rec_corrs_CT_start_index-1][0]);
    
    double fcorr=Rec_corrs_CT_Data[Rec_corrs_CT_start_index-1][cflag]+df_dz*dz; 
    return 1.0+fcorr;
}

//========================================================================================
// correction factor 1+DNe/Ne as explained by Rubino-Martin et al. 2010.
//----------------------------------------------------------------------------------------
// Modification according to Hart & Chluba, 2017
// cflag == 1: fH+fHeI (default)
// cflag == 2: fH
// cflag == 3: fHeI
//========================================================================================
double Recombination_correction_factor(double z, int cflag)
{
    int nz=Rec_corrs_CT_Data.size();
    
    if(Rec_corrs_CT_Data[0][0]<=z) return 1.0;
    if(Rec_corrs_CT_Data[nz-1][0]>=z) return 1.0+Rec_corrs_CT_Data[nz-1][cflag];
    return interpolate_data(z, cflag);
}

//========================================================================================
//========================================================================================
