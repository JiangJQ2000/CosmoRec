//===========================================================================================================
// Author: Jens Chluba 
// last modification: Nov 2010
// purpose: to access the main CosmoRec code
//===========================================================================================================

#ifndef COSMOREC_H
#define COSMOREC_H

//===========================================================================================================
// CosmoRec for run from console
//===========================================================================================================
int CosmoRec(int narg, char *args[]);

//===========================================================================================================
// to call CosmoRec with a list of parameters in batch mode (e.g. when calling it from CosmoMC).
// 
// runmode == 0: CosmoRec run with diffusion
//            1: CosmoRec run without diffusion
//            2: Recfast++ run
//            3: Recfast++ run with correction function of Chluba & Thomas, 2010
//
// On entry, the array z_arr should contain the redshifts at which Xe and Te are required. nz 
// determines the number of redshift points. Xe_arr and Te_arr will contain the solution on exit.
//
// Furthermore, runpars[0] defines the dark matter annihilation efficiency in eV/s.
// runpars[1] switches the accuracy of the recombination model:
//
// runpars[1]==-1: closest equivalent of 'HyRec' case (Haimoud & Hirata, 2010)
// runpars[1]== 0: default setting
// runpars[1]== 1: 2g for n<=4 & Raman for n<=3
// runpars[1]== 2: 2g for n<=8 & Raman for n<=7
// runpars[1]== 3: 2g for n<=8 & Raman for n<=7 + Helium feedback up to n=5 (full setting)
//
// The value of runpars[1] is only important for runmode 0 & 1.
//
// By default CosmoRec writes into a file the ionisation history, solutions for the
// populations of the resolved levels, and some information about the Cosmology (usually 
// directed into the directory "./outputs/"). This can cause problems when using CosmoRec 
// with a parallel CosmoMC, however, it possible to controll the outputs through runpars[2]:
// 
// runpars[2]==0: write out the recombination history, populations, and the cosmology.
// runpars[2]==1: write out the recombination history, and the cosmology.
// runpars[2]==2: write out only the recombination history.
// runpars[2]==3: don't write out anything.
//
//===========================================================================================================
int CosmoRec(const int runmode, const double runpars[5], 
             const double omegac, const double omegab, 
             const double omegak, const double Nnu,  
             const double h0, const double tcmb, const double yhe, 
             const int nz, double *z_arr, double *Xe_arr, double *Te_arr,
             const int label);

//===========================================================================================================
// Wrap the C++ Fortran routine to be allow calling from Fortran. Arguments are as above.
// Added 06.03.2011 (Richard Shaw)
//===========================================================================================================
extern "C" {
    
    void cosmorec_fortran_(int * runmode, double * runpars, double * omega_c, double * omega_b, double * omega_k, 
                           double * num_nu, double * h0, double * t_cmb, double * y_he, 
                           double * za_in, double * xe_out, double * tb_out, int * len, int* label);
}

//===========================================================================================================
// to cleanup after finishing with CosmoRec
//===========================================================================================================
void cleanup_CosmoRec();

//===========================================================================================================
#endif
