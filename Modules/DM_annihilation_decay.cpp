//========================================================================================
// Author: Jens Chluba (Nov. 19th 2010)
// CITA, University of Toronto
// All rights reserved.
//========================================================================================
// 21.05.2018: added control of extra excitations and reorganized function. Also added
//             a simple decaying particle module [JC]

//========================================================================================
// This simple module allows to add the ionizing effect of DM-annihilation / decay.
// The important equations as implemented here were given in Chluba 2010. Original works
// are Chen & Kamionkowski 2004 and Padmanabhan & Finkbeiner 2005.
//========================================================================================

//========================================================================================
// branching into ionizations according to Chen & Kamionkowski 2004
// (see also Chluba 2010)
//========================================================================================
double branching_ratio_ions_Chen(double x){ return (1.0-x)/3.0; }
double branching_ratio_ex_Chen  (double x){ return (1.0-x)/3.0; }
double branching_ratio_heat_Chen(double xp, double ZHeII, double fHe)
{ return (1.0 + 2.0*xp + fHe*(1.0 + 2.0*ZHeII))/3.0/(1.0+fHe); }


//========================================================================================
// branching into ionizations according to Shull & van Steenberg 1985 
// (see also Chluba 2010)
//========================================================================================
double branching_ratio_ions_SvSt(double x)
{ return 0.4*pow(1.0-pow(x, 0.4), 7.0/4.0); }


//========================================================================================
// evalutation of ODE terms
//========================================================================================
void evaluate_ionization_heating_excitation_terms(double fHe, double Tg,
                                                  double XHII, double XHeII, double dE_dt,
                                                  double &dXHI_1s_dt, double &dXHeI_1s_dt,
                                                  double &drho_dt)
{
    //------------------------------------------------------------------------------------
    // different branching ratios. Simple approximations here but this can be made more
    // comprehensive using Slatyer et al., if needed
    //------------------------------------------------------------------------------------
    double bHe_ion=branching_ratio_ions_Chen(XHeII/fHe);
    double bHI_ion=branching_ratio_ions_Chen(XHII);
    double b_heat =branching_ratio_heat_Chen(XHII, XHeII/fHe, fHe);
    // for temperature equation
    double fheat=2.0/3.0*const_e/const_kB/(1.0+fHe+XHII+XHeII);
    
    //====================================================================================
    // return different terms (HI, HeI, Te)
    //====================================================================================
    dXHI_1s_dt+=-bHI_ion/13.6/(1.0+fHe)*dE_dt;
    if(XHeII>0.0) dXHeI_1s_dt+=-bHe_ion/24.6*fHe/(1.0+fHe)*dE_dt;
    drho_dt+=fheat*b_heat*dE_dt/Tg;

    //====================================================================================
    // additional excitations
    //====================================================================================
    if(parameters.CR.add_excitations)
    {
        double bHe_ex=branching_ratio_ex_Chen(XHeII/fHe);
        double bHI_ex=branching_ratio_ex_Chen(XHII);
        
        dXHI_1s_dt+=-bHI_ex/10.2/(1.0+fHe)*dE_dt;
        if(XHeII>0.0) dXHeI_1s_dt+=-bHe_ex/21.0*fHe/(1.0+fHe)*dE_dt;
    }

    return;
}

//========================================================================================
// evalutation of ODE terms for annihilation cases
//========================================================================================
void evaluate_DM_annihilation_terms(double z, double fHe, double Tg,
                                    double XHII, double XHeII,
                                    double &dXHI_1s_dt, double &dXHeI_1s_dt,
                                    double &drho_dt)
{
    //====================================================================================
    // fDM [eV/s] gives annihilation efficiency; typical value fDM=2.0e-24 eV/s
    // (see Chluba 2010 for definitions)
    // In this module additional excitations are not included at the moment
    //====================================================================================
    double fDM=parameters.cosm.Rec.f_ann;
    double Fcal=pow(1.0+z, 3);
    double dE_dt=fDM*Fcal;  // here the factor 1/NH was canceled with definition of dE_dt
    
    evaluate_ionization_heating_excitation_terms(fHe, Tg, XHII, XHeII, dE_dt,
                                                 dXHI_1s_dt, dXHeI_1s_dt, drho_dt);
    
    return;
}

//========================================================================================
// evaluation of ODE terms for decay cases. Here things are expressed relative to the
// Dark matter abundance with fdec controling the fraction of decay.
//========================================================================================
void evaluate_DM_decay_terms(double z, double fHe, double Tg,
                             double XHII, double XHeII,
                             double &dXHI_1s_dt, double &dXHeI_1s_dt,
                             double &drho_dt)
{
    //====================================================================================
    // total energy deposition in eV per hydrogen nucleus per sec
    //====================================================================================
    double fdec=parameters.cosm.Rec.f_dec;
    double Gamma_dec=parameters.cosm.Rec.Gamma_dec;
    double rho_cdm_NH=cosmos.Get_O_DM()/cosmos.Get_O_b()*(1.0-cosmos.Get_Y_p())*938.28e+6;
    double f_eff=0.1; // energy transfer efficiency factor
    double dE_dt=rho_cdm_NH * f_eff * fdec * Gamma_dec * exp(-Gamma_dec*cosmos.t(z));

    evaluate_ionization_heating_excitation_terms(fHe, Tg, XHII, XHeII, dE_dt,
                                                 dXHI_1s_dt, dXHeI_1s_dt, drho_dt);
    
    return;
}

//========================================================================================
//========================================================================================
