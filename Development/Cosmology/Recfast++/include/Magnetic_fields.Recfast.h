//========================================================================================
// Author: Jens Chluba (Oct. 2014)
// Johns Hopkins University
//========================================================================================
// 30.03.2017: added set functions for parameters related to PMF [JC]

//========================================================================================
// This simple module allows to add the effect of heating by magnetic fields.
// Expressions from Sethi & Subramanian, 2005 and Kunze & Komatsu, 2013 are used.
//========================================================================================

#ifndef MAGNETIC_FIELDS_H
#define MAGNETIC_FIELDS_H

//========================================================================================
// Settings and parameter access
//========================================================================================
void Switch_PMF_off();
void Set_PMF(bool decay_on, bool ambi_on, bool Lorentz, double B0, double nB, double kD);

double Get_B0_PMF();
double Get_nB_PMF();
double Get_kD_PMF();
bool Get_decay_setting();
bool Get_ambi_setting();
bool Get_Lorentz_setting();

//========================================================================================
// ODE terms
//========================================================================================
void evaluate_magnetic_field_heating_RF(double z, double Hz, double fHe, double NH,
                                        double xp, double xHep, double Tm, double nbar,
                                        double B0, double nB, double kD, double fvec[]);

#endif

//========================================================================================
//========================================================================================
