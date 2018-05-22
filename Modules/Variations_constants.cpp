//===========================================================================================================
//
// 
// Variation of fundamental constants for use in CosmoRec
// Author: Luke Hart (Dec 2016)
//
//
//===========================================================================================================
void Switch_VFC_off()
{
    // set CosmoRec parameters
    FC_CosmoRec.a_scale = FC_CosmoRec.m_scale = 1.0;
    FC_CosmoRec.p_scale = 0.0;
    FC_CosmoRec.runmode = 0;

    // set Recfast++ parameters
    Cosmos::parameters *Cpars=&parameters.cosm;
    Cpars->Rec.aS=Cpars->Rec.mS=1.0; Cpars->Rec.pS=0.0; Cpars->Rec.RM=0;
    Cpars->Rec.do_separate_H_He=0;

    return;
}

void Set_VFC_pars(double aS, double mS, double pS, int RM)
{
    if(RM==0){ Switch_VFC_off(); return; }
    if(aS<0 && mS<0){ Switch_VFC_off(); return; }
    
    // set CosmoRec parameters
    FC_CosmoRec.a_scale =(aS<=0.0 ? 1.0 : aS);
    FC_CosmoRec.m_scale =(mS<=0.0 ? 1.0 : mS);
    FC_CosmoRec.p_scale = 0.0; // not really used yet in CosmoRec! [JC]
    FC_CosmoRec.runmode = RM;
    
    // final check if mode activated
    if(aS==1.0 && mS==1.0 && pS==0.0){ Switch_VFC_off(); return; }
    
    // set Recfast++ parameters
    Cosmos::parameters *Cpars=&parameters.cosm;
    Cpars->Rec.aS=FC_CosmoRec.a_scale;
    Cpars->Rec.mS=FC_CosmoRec.m_scale;
    Cpars->Rec.pS=FC_CosmoRec.p_scale;
    Cpars->Rec.RM=RM;
    Cpars->Rec.do_separate_H_He=1;

    return;
}

//===========================================================================================================
// Rescaling functions for all the rates used in CosmoRec
//===========================================================================================================
void rescale_variables_all(double alpha_scale, double me_scale, size_t run_mode,
                           double &Ascale, double &Bscale, double &Tscale, int &use_He)
{
    //=======================================================================================================
    // This function is designed to rescale all the explicit scalings we need to worry about and drop this
    // into the reading of the file. For in depth details about the omission mode please look at
    // evalode.Recfast.cpp description of the variables
    //=======================================================================================================
    // Ignore all scalings if the alpha and mu values are already set to 1
    //=======================================================================================================
    const_HI_A2s_1s = const_HI_A2s_1s_0;
    const_HeI_A2s_1s = const_HeI_A2s_1s_0;
    const_sigT = const_sigT_0;
    Ascale=Bscale=Tscale=1.0;
    
    if(run_mode == 0) return;
        
    if (alpha_scale != 1. || me_scale != 1.) {
        
        // Rescaled Thomson scattering
        if (run_mode == 2 || run_mode == 6) {
            const_sigT = const_sigT_0*pow(alpha_scale/me_scale, 2);
        }
        
        // Rescaled 2 photon rate
        if (run_mode == 3 || run_mode == 6) {
            const_HI_A2s_1s = const_HI_A2s_1s_0*pow(alpha_scale,8)*me_scale;
            if (use_He == 1) const_HeI_A2s_1s = const_HeI_A2s_1s_0*pow(alpha_scale,8)*me_scale;
        }
        
        // Rescaled Boltzmann factors
        if (run_mode == 1 || run_mode == 6) {
            Tscale = pow(alpha_scale,-2)*pow(me_scale,-1);
        }
        
        // Rescaled recomb/photoion rates
        if (run_mode == 4 || run_mode == 6) {
            Ascale = pow(alpha_scale*me_scale,-3);
            Bscale = pow(alpha_scale,5)*me_scale;
        }
        
        // Rescaled Pesc
        if (run_mode == 5)
            throw_error("rescale_variables_all", "This CosmoRec mode does not exist...");
    }
    
    // rescale hydrogen atom. This is done interally and relative to whatever state the atom is in [LH & JC]
    if (run_mode == 1 || run_mode == 6) 
    {
        Hydrogen_Atoms.rescale_gas(alpha_scale, me_scale);
        if (use_He == 1) HeI_Atoms.rescale_gas(alpha_scale, me_scale);
    }
    
    // set rescalings for PDE problem [LH & JC]
    if (run_mode == 6) 
    {
        Set_rescaling_PDE_HI(alpha_scale, me_scale);
        if (use_He == 1) Set_rescaling_PDE_HeI(alpha_scale, me_scale);
    }
    
    return;
}

//===========================================================================================================
// Rescaling the rates
// We need to rescale the Ai Bi and Rij rates by our coefficients and this function does that
// independently since we need to do this for a series of vectors
// This function is overloaded for 2 vectors and a vector of vectors as well as a single vector
//===========================================================================================================
void rescale_rates(vector<double>&A,vector<double>&B,vector<vector<double> >&R)
{
    if(FC_CosmoRec.f_a==1.0 || FC_CosmoRec.f_b==1.0) return;  // save time if not needed
    
    // For each of the rates we're going to go through the iterations and times every element by the
    // f_i prefactors
    // Leave this alone as this works
    for (vector<double>::iterator a1 = A.begin(); a1 != A.end(); ++a1) *a1 *= FC_CosmoRec.f_a;
    for (vector<double>::iterator b1 = B.begin(); b1 != B.end(); ++b1) *b1 *= FC_CosmoRec.f_b;
    for (vector<vector<double> >::iterator r1 = R.begin(); r1 != R.end(); ++r1) {
        for (vector<double>::iterator r2 = r1->begin(); r2 != r1->end(); ++r2) {
            *r2 *= FC_CosmoRec.f_b;
        }
    }
    return;
}

//===========================================================================================================
//===========================================================================================================
