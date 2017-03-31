//========================================================================================
// Author: Jens Chluba
// Last modification: Oct 2010
// CITA, University of Toronto
// All rights reserved.
//========================================================================================

//========================================================================================
// This module contains the evaluation of the RECFAST ODE-system. The code is based on the
// C-version for RECFAST from Seager 1999
//========================================================================================

#ifndef EVALODE_H
#define EVALODE_H

//========================================================================================
// Photo-ionization with T=Te (Recfast default, although incorrect) or T=Tg (phys.correct)
//========================================================================================
void Set_Te_eval();
void Set_Tr_eval();
void Set_default(); // --> T=Te
bool Get_T_eval();

//========================================================================================
// evaluation of Recfast-system
//========================================================================================
void evaluate_Recfast_System(double z, double y[], double fvec[], int n);

//========================================================================================
// wrapper for ODE-solver
//========================================================================================
void fcn(int *neq, double *z, double *y, double *f);

//========================================================================================
// wrapper for ODE-solver with rescaled derivative for Xe
// this is important to stitch the recfast solution to the output of the multi-level code
//========================================================================================
void set_rescaling(double ff);
void fcn_rescaled(int *neq, double *z, double *y, double *f);

#endif
