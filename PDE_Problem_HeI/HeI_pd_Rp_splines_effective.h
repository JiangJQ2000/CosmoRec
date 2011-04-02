//===========================================================================================================
// Authors Jeffrey Fung & Jens Chluba Feb/March 2010
//
// comment: HeI populations have to be pre-loaded; Also the effective rates have to be initialized
//===========================================================================================================

#ifndef HEI_PD_RP_SPLINES_EFFECTIVE_H
#define HEI_PD_RP_SPLINES_EFFECTIVE_H

#include <string>
#include "HeI_Atom.h"
#include "Cosmos.h"

//===========================================================================================================
// For precalulated death probabilities and emission rates
//===========================================================================================================
void set_up_splines_for_HeI_pd_Rp_effective(double zend, double zstart, Cosmos &cosmos, Gas_of_HeI_Atoms &HeA);

void clear_HeI_pd_Rp_effective_memory();
void Set_HeI_pd_Rp_splines_effective_verbosity(int v);

//===========================================================================================================
// pd-functions (the index ires is the corresponding resolved level index)
//===========================================================================================================
double calc_HeI_pd_i_splines_effective(double z, int ires);

//===========================================================================================================
// Rp-functions (the index ires is the corresponding resolved level index)
//===========================================================================================================
double calc_HeI_Rp_Rm_i_splines_effective(double z, int ires);

//===========================================================================================================
// Dnem-functions (the index ires is the corresponding resolved level index)
//===========================================================================================================
double calc_HeI_Dnem_i_splines_effective(double z, int ires);

#endif
