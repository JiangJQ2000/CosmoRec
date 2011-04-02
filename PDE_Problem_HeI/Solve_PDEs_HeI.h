//==================================================================================================
// Authors Jeffrey Fung & Jens Chluba Feb/March 2010
//==================================================================================================

#ifndef SOLVE_PDES_HEI_H
#define SOLVE_PDES_HEI_H

#include "Cosmos.h"
#include "Atom.h"
#include "HeI_Atom.h"

int compute_HeI_DPesc_with_diffusion_equation_effective(vector<double> &DF_vec_z, 
                                                        vector<vector<double> > &DF_vec,
                                                        vector<int> &HeI_level_indices,
                                                        double zs, double ze, 
                                                        Cosmos &cos, 
                                                        Gas_of_HeI_Atoms &HeIAtoms,
                                                        vector<vector<double> > HI_HeI_Solution, 
                                                        int it_num);

#endif
