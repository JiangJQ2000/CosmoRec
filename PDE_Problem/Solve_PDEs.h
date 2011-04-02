//======================================================================================
// Author Jens Chluba (July 2010)
//======================================================================================

#ifndef SOLVE_PDES_H
#define SOLVE_PDES_H

#include "Cosmos.h"
#include "Atom.h"

//======================================================================================
// this allows to switch the 1s-2s absorption correction off; in the radiative transfer 
// everything is included, but the integral is only taken over the emission term
//======================================================================================
void switch_off_1s_2s_absorption_correction();

int compute_DPesc_with_diffusion_equation_effective(vector<double> &DF_vec_z, 
                                                    vector<vector<double> > &DF_Ly_n, 
                                                    vector<vector<double> > &DF_2_gamma_vec, 
                                                    vector<vector<double> > &DF_Raman_vec, 
                                                    vector<double> &DI1_2s_vec, 
                                                    double zs, double ze, 
                                                    int nmax_2g_corrs, int nmax_R_corrs,
                                                    Cosmos &cos, Gas_of_Atoms &HIA, 
                                                    vector<vector<double> > HI_Solution, 
                                                    int it_num);

#endif
