//========================================================================================
// Author: Jens Chluba
// Last modification: May 2011
// CITA, University of Toronto
// All rights reserved.
//========================================================================================

//========================================================================================
// This simple module allows to add the effect of DM-annihilation. The important 
// equations as implemented here were given in Chluba 2010. Original works are 
// Chen & Kamionkowski 2004 and Padmanabhan & Finkbeiner 2005.
// In this module additional excitations can be switched on/off
//========================================================================================

#include <cmath>
#include <iostream>
#include <cstdlib>

#include "constants.Recfast.h"
#include "DM_annihilation.Recfast.h"

using namespace std;
using namespace RECFAST_physical_constants; // defined in constants.Recfast.h

bool RECFAST_add_excitations=0; // switch the effect of extra excitations on/off

//========================================================================================
// branching into ionizations according to Chen & Kamionkowski 2004
// (see also Chluba 2010)
//========================================================================================
double branching_ions_Chen(double x){ return (1.0-x)/3.0; }
double branching_ex_Chen  (double x){ return (1.0-x)/3.0; }
double branching_heat_Chen(double xp, double ZHeII, double fHe)
{ return (1.0 + 2.0*xp + fHe*(1.0 + 2.0*ZHeII))/3.0/(1.0+fHe); }



//========================================================================================
// branching into ionizations according to Shull & van Steenberg 1985 
// (see also Chluba 2010)
//========================================================================================
double branching_ions_SvSt(double x){ return 0.4*pow(1.0-pow(x, 0.4), 7.0/4.0); }



//========================================================================================
// evalutation of ODE terms
//========================================================================================
void evaluate_DM_annihilation_terms(double z, double Hz, double fHe, double xp, 
                                    double xHep, double CHe, double CH, double fDM, 
                                    double fvec[])
{
    //====================================================================================
    // fDM [eV/s] gives annihilation efficiency; typical value fDM=2.0e-24 eV/s
    // (see Chluba 2010 for definitions)
    //====================================================================================    
    double Fcal=(1.0+z)*(1.0+z)/Hz;
    double dE_dz=fDM*Fcal;  // here the factor 1/NH was canceled with definition of dE_dt
    //
    double bHe_ion=branching_ions_Chen(xHep/fHe);
    double bHI_ion=branching_ions_Chen(xp);
    double b_heat=branching_heat_Chen(xp, xHep/fHe, fHe);
    //
    double fheat=2.0/3.0*1.602176487e-12/(RF_kBoltz*1.0e+7)/(1.0+fHe+xp+xHep);
    
    //====================================================================================
    // return different terms
    //====================================================================================
    fvec[0]-=bHe_ion/24.6*fHe/(1.0+fHe)*dE_dz;
    fvec[1]-=bHI_ion/13.6    /(1.0+fHe)*dE_dz;
    fvec[2]-=fheat*b_heat*dE_dz;

    //====================================================================================
    // additional excitations (added May 2011)
    //====================================================================================
    if(RECFAST_add_excitations)
    {
        double bHe_ex=branching_ex_Chen(xHep/fHe);
        double bHI_xe=branching_ex_Chen(xp);

        fvec[0]-=(1.0-CHe)*bHe_ex/21.0*fHe/(1.0+fHe)*dE_dz;
        fvec[1]-=(1.0-CH )*bHI_xe/10.2    /(1.0+fHe)*dE_dz;
    }
    
    return;
}
