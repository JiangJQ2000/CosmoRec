//========================================================================================
// Author: Jens Chluba (Oct. 2014)
// Johns Hopkins University
//========================================================================================

//========================================================================================
// This simple module allows to add the effect of heating by magnetic fields.
// Expressions from Sethi & Subramanian, 2005 and Kunze & Komatsu, 2013 are used.
//========================================================================================

#ifndef MAGNETIC_FIELDS_H
#define MAGNETIC_FIELDS_H

void evaluate_magnetic_field_heating_RF(double z, double Hz, double fHe, double NH,
                                        double xp, double xHep, double Tm, double nbar,
                                        double B0, double nB, double kD, double fvec[]);

#endif

//========================================================================================
//========================================================================================
