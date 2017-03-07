//==================================================================================================
// Hydrogen atom
//==================================================================================================
void update_H_populations(const double *XiH, Gas_of_Atoms &H_Atoms)
{
    // XiH[0..nH] contains ALL the populations of levels
    H_Atoms.Set_population_of_level(0, XiH[0]);
    for(int i=0; i<get_number_of_resolved_levels(); i++) 
        H_Atoms.Set_population_of_level(get_HI_index(i), XiH[get_HI_index(i)]);
    
    return;
}

//==================================================================================================
// neutral Helium setup
//==================================================================================================
void update_HeI_populations(const double *XiHeI, Gas_of_HeI_Atoms &HeIA)
{
    // XiHeI[0..nH] contains ALL the populations of levels
    HeIA.Set_Xi(0, XiHeI[0]);
    for(int i=0; i<get_number_of_resolved_levels_HeI(); i++) 
        HeIA.Set_Xi(get_HeI_index(i), XiHeI[get_HeI_index(i)]);
    
    return;
}

//==================================================================================================
// compute Xe, Xp, XHeI, XHeII, XHeIII
//==================================================================================================
void compute_fractions(const Data_Level_I &LI, double &Xe, double &Xp, double &XHeII, double &XHeIII)
{
    int iHI=LI.index_HI, iHeI=LI.index_HeI;
    
    Xp=1.0-LI.X[iHI];
    XHeII=XHeIII=0.0;
    
    if(flag_He==1)
    {
        XHeII =parameters.fHe-LI.X[iHeI];
        XHeIII=0.0;
    }

    Xe=Xp+XHeII+2.0*XHeIII;
    
    return;
}

//==================================================================================================
//==================================================================================================
static double ODE_Solver_z_pre=0.0, ODE_Solver_rho_pre=0.0;

//==================================================================================================
void update_Solver_vars_eff(const double *XiH, Gas_of_Atoms &H_Atoms, 
                            const double *XiHeI, Gas_of_HeI_Atoms &HeIA, 
                            double z, double rho, int Pop_flag=1)
{
    // update populations
    if(Pop_flag==1)
    {
        update_H_populations(XiH, H_Atoms);
        if(flag_He==1) update_HeI_populations(XiHeI, HeIA);
    }
    
    if(z!=ODE_Solver_z_pre || rho!=ODE_Solver_rho_pre)
    {
        ODE_Solver_z_pre=z;
        ODE_Solver_rho_pre=rho;
        
        // effective rates
        double Tg=cosmos.TCMB(z);
        get_rates(Tg*FC_CosmoRec.f_t, Tg*rho*FC_CosmoRec.f_t, Ai_df_dx, Bi_df_dx, RijVec_df_dx);
        if(FC_CosmoRec.runmode>0) rescale_rates(Ai_df_dx, Bi_df_dx, RijVec_df_dx);
        
        // equivalent for Helium
        if(flag_He==1) {
            get_rates_HeI(Tg*FC_CosmoRec.f_t, Ai_df_dx_HeI, Bi_df_dx_HeI, RijVec_df_dx_HeI);
            if(FC_CosmoRec.runmode>0) rescale_rates(Ai_df_dx_HeI, Bi_df_dx_HeI, RijVec_df_dx_HeI);
        }
           
    }
    
    return;
}

//==================================================================================================
// aux-functions for variable mapping
//==================================================================================================
void copy_LI_to_ysol(double *LIX, double *ysol, int neq, int index_HI, int index_HeI)
{
    int nres_HI=get_number_of_resolved_levels();
    int nres_HeI=get_number_of_resolved_levels_HeI();
    
    //============================================================================
    // 1. eq is Te
    //============================================================================
    ysol[0]=LIX[neq-1];
    ysol[1]=LIX[index_HI];
    
    for(int k=0; k<nres_HI; k++) ysol[k+2]=LIX[index_HI+get_HI_index(k)];
    
    if(flag_He>=1)
    {
        ysol[nres_HI+2]=LIX[index_HeI];
        for(int k=0; k<nres_HeI; k++) ysol[nres_HI+k+3]=LIX[index_HeI+get_HeI_index(k)];
    }
    
    return;
}

//==================================================================================================
void copy_ysol_to_LI(double *ysol, Data_Level_I &LI, int neq, int index_HI, int index_HeI)
{
    int nres_HI=get_number_of_resolved_levels();
    int nres_HeI=get_number_of_resolved_levels_HeI();
    
    //============================================================================
    // 1. eq is Te
    //============================================================================
    LI.X[neq-1]=ysol[0];
    LI.X[index_HI]=ysol[1];
    
    for(int k=0; k<nres_HI; k++) LI.X[index_HI+get_HI_index(k)]=ysol[k+2];
    
    if(flag_He>=1)
    {
        LI.X[index_HeI]=ysol[nres_HI+2];
        for(int k=0; k<nres_HeI; k++) LI.X[index_HeI+get_HeI_index(k)]=ysol[nres_HI+k+3];
    }
    
    // also update Xe
    double Xe, Xp, XHeII, XHeIII;
    compute_fractions(LI, Xe, Xp, XHeII, XHeIII);
    LI.X[0]=Xe;

    return;
}

//==================================================================================================
// fcn for hydrogen
//==================================================================================================
void fcn_HI_effective(double z, double rho, double Xe, double Tg, double NH, double H_z,
                      double Xp, int iHI, Data_Level_I &LI)
{
    //==================================================================
    // effective rates
    //==================================================================
    ODE_effective::evaluate_2s_two_photon_decay(z, Tg, Hydrogen_Atoms.X(0), Hydrogen_Atoms.X(1),
                                                Hydrogen_Atoms.Level(1).Get_Dnu_1s(),
                                                const_HI_A2s_1s, LI.g[iHI], LI.g[iHI+1]);
    
    //==================================================================
    // HI Lyman-series
    //==================================================================
    int nres=get_number_of_resolved_levels();
    int nLy_max=(int)min(Lyn_max_parameter, Hydrogen_Atoms.Get_n_of_Level(get_HI_index(nres-1)));
    
    for(int n=2; n<=nLy_max; n++)
    {
        Transition_Data_A T=Hydrogen_Atoms.Level(n, 1).Get_Data_A(1, 0);
        
        ODE_effective::evaluate_Ly_n_channel(z, Tg, Hydrogen_Atoms.X(0), Hydrogen_Atoms.X(n, 1),
                                             NH, H_z, T.A21, T.lambda21, T.Dnu,
                                             LI.g[iHI],
                                             LI.g[iHI+Hydrogen_Atoms.Get_Level_index(n,1)]);
    }
    
    //==================================================================
    // HI nD-1s Quadrupole transitions (added May 2011)
    //==================================================================
    if(Atom_activate_HI_Quadrupole_lines)
    {
        for(int k=0; k<get_number_of_resolved_levels(); k++)
        {
            int l=Hydrogen_Atoms.Get_l_of_Level(get_HI_index(k));
            
            if(l==2)
            {
                int n=Hydrogen_Atoms.Get_n_of_Level(get_HI_index(k));
                
                if(n<=nD_Quadrupole_max)
                    ODE_effective::evaluate_nD_1s_Q_channel(Hydrogen_Atoms.X(n, 1),
                                                            Hydrogen_Atoms.X(n, 2),
                                                            Hydrogen_Atoms.Level(n, 2).Get_A21(1, 0),
                                                            LI.g[iHI+Hydrogen_Atoms.Get_Level_index(n, 1)],
                                                            LI.g[iHI+Hydrogen_Atoms.Get_Level_index(n, 2)]);
            }
        }
    }
    
    //==================================================================
    // these are the ones with effective rates that depend on Te
    //==================================================================
    ODE_HI_effective::evaluate_effective_Rci_Ric_terms(z, Xe, NH*Xp, Hydrogen_Atoms, Ai_df_dx,
                                                       Bi_df_dx, get_HI_index, &LI.g[iHI]);
    
    ODE_HI_effective::evaluate_effective_Rij_terms(z, Tg, Hydrogen_Atoms, RijVec_df_dx,
                                                   get_HI_index, &LI.g[iHI]);
    
    //==================================================================
    // Diffusion correction
    //==================================================================
    if(Diffusion_correction_is_on)
    {
        evaluate_HI_Diffusion_correction(z, Hydrogen_Atoms, NH, H_z, Tg, &LI.g[iHI]);
        
        if(DI1_2s_correction_on) evaluate_HI_DI1_DI2_correction(z, Hydrogen_Atoms.X(0),
                                                                Hydrogen_Atoms.X(1),
                                                                Hydrogen_Atoms.Level(1).Get_Dnu_1s(),
                                                                Tg, &LI.g[iHI]);
    }
    
    return;
}

//==================================================================================================
// fcn for helium
//==================================================================================================
void fcn_HeI_effective(double z, double Xe, double Tg, double NH, double H_z,
                       double XHeII, int iHeI, Data_Level_I &LI)
{   
    double Te=Tg;   
    
    //==================================================================
    // these are the ones with effective rates that depend on Xi
    //==================================================================
    ODE_HeI_effective::evaluate_effective_Rci_Ric_terms(z, Xe, NH*XHeII, HeI_Atoms, Ai_df_dx_HeI,
                                                        Bi_df_dx_HeI, get_HeI_index, &LI.g[iHeI]);
    
    ODE_HeI_effective::evaluate_effective_Rij_terms(z, Tg, HeI_Atoms, RijVec_df_dx_HeI,
                                                    get_HeI_index, &LI.g[iHeI]);    

    //==========================================================================
    // 2s two-photon channels (2S - 1S singlet channel)
    //==========================================================================
    ODE_effective::evaluate_2s_two_photon_decay(z, Tg, HeI_Atoms.Xi(0), HeI_Atoms.Xi(1),
                                                HeI_Atoms.Sing.Level(2, 0).Get_Dnu_1s2(),
                                                const_HeI_A2s_1s, LI.g[iHeI], LI.g[iHeI+1]);
    
    //==================================================================
    // HeI 'Lyman-series' & quadrupole lines (if switched on)
    //==================================================================
    for(int k=1; k<get_number_of_resolved_levels_HeI(); k++)
    {
        int ik=get_HeI_index(k);
        if(flag_spin_forbidden==0 && ik==(int)HeI_Atoms.Get_Level_index(2, 1, 1, 1)) continue;
        
        Transition_Data_HeI_A T=HeI_Atoms.Get_Trans_Data(ik, 1, 0, 0, 0);
        if(T.A21!=0) ODE_effective::evaluate_Ly_n_channel(z, Tg, HeI_Atoms.Xi(0), HeI_Atoms.Xi(ik),
                                                          NH, H_z, T.A21, T.lambda21, T.Dnu, 
                                                          LI.g[iHeI], LI.g[iHeI+ik], 
                                                          HeI_Atoms.Get_gw(ik)/T.gwp);
    }
    
    //==================================================================
    // HeI-diffusion correction (added 22.01.2011; JF & JC)
    //==================================================================
    if(Diffusion_correction_HeI_is_on)
    {
        evaluate_HeI_Diffusion_correction(z, HeI_Atoms, NH, H_z, Tg, 
                                          &LI.g[iHeI], LI.g[0], LI.g[LI.index_HI]);
        
        if(DI1_2s_correction_on_HeI) 
            evaluate_HeI_DI1_DI2_correction(z, HeI_Atoms.Xi(0), HeI_Atoms.Xi(1), 
                                            HeI_Atoms.Sing.Level(2, 0).Get_Dnu_1s2(), 
                                            Tg, &LI.g[iHeI]);
    }
    
    //==================================================================
    // HeI-feedback
    //==================================================================
    if(HeISTfeedback==1 && !Diffusion_correction_HeI_is_on)
    {
        if(_HeI_feedback_w_HI_abs==0) 
            evaluate_HeI_feedback(z, Hydrogen_Atoms.X(1, 0), HeI_Atoms, cosmos, HeI_Lines,
                                  Dn_HeI_feedback_no_HI_abs, &LI.g[iHeI]);
        
        if(_HeI_feedback_w_HI_abs==1) 
            evaluate_HeI_feedback(z, Hydrogen_Atoms.X(1, 0), HeI_Atoms, cosmos, HeI_Lines,
                                  Dn_HeI_feedback, &LI.g[iHeI]);
    }
    
    //==================================================================
    // absorption by neutral hydrogen
    //==================================================================
    if(flag_HI_absorption==1 && z<=zcrit_HI*FC_CosmoRec.f_t && !Diffusion_correction_HeI_is_on)
    {  
        double dXe_HeI_HI_abs_dt=0.0;
        //==============================================================
        // Ly-a line
        //==============================================================
        evaluate_HI_abs_HeI(z, Hydrogen_Atoms.X(1, 0), NH, H_z, Tg, Te,
                            HeI_Atoms, cosmos, 
                            HeI_Atoms.HILyc, DP_interpol_S, dXe_HeI_HI_abs_dt);

        // Helium levels
        LI.g[iHeI+0]+= dXe_HeI_HI_abs_dt;
        LI.g[iHeI+2]+=-dXe_HeI_HI_abs_dt;
        // electron equation
        LI.g[0]+=dXe_HeI_HI_abs_dt;
        // HI 1s
        LI.g[LI.index_HI]+=-dXe_HeI_HI_abs_dt;
        
        //==============================================================
        // from intercombination line
        //==============================================================
        if(z<=zcrit_HI_Int*FC_CosmoRec.f_t && flag_spin_forbidden==1)
        {
            dXe_HeI_HI_abs_dt=0.0;
            evaluate_HI_abs_HeI_Intercombination(z, Hydrogen_Atoms.X(1, 0), NH, H_z, Tg, Te,
                                                 HeI_Atoms, cosmos, 
                                                 HeI_Atoms.HILyc, DP_interpol_T, 
                                                 dXe_HeI_HI_abs_dt);        
            
            // Helium levels
            LI.g[iHeI+0]+= dXe_HeI_HI_abs_dt;
            LI.g[iHeI+HeI_Atoms.Get_Level_index(2, 1, 1, 1)]+=-dXe_HeI_HI_abs_dt;
            // electron equation
            LI.g[0]+=dXe_HeI_HI_abs_dt;
            // HI 1s
            LI.g[LI.index_HI]+=-dXe_HeI_HI_abs_dt;
        }
    }
    
    return;
}

//==================================================================================================
// fcn for ODE evaluation
//==================================================================================================
void fcn_effective(double z, Data_Level_I &LI)
{
    double rho=LI.X[LI.neq-1];
    double Tg=cosmos.TCMB(z);
    int iHI=LI.index_HI;
    int iHeI=LI.index_HeI;
    
    double Xe, Xp, XHeII, XHeIII;
    compute_fractions(LI, Xe, Xp, XHeII, XHeIII);
    double NH=cosmos.NH(z), H_z=cosmos.H(z);
    
    //==================================================================
    // copy Xi_H --> Atoms && update recombination coefficients
    //==================================================================
    update_Solver_vars_eff(&LI.X[iHI], Hydrogen_Atoms, 
                           &LI.X[iHeI], HeI_Atoms,
                           z, rho);
    
    //==================================================================
    // reset everything 
    //==================================================================
    for(int i=0; i<LI.neq; i++) LI.g[i]=0.0;
    
    //==================================================================
    // Te part
    //==================================================================
    ODE_effective::evaluate_TM(z, Xe, parameters.fHe, rho, Tg, H_z, LI.g[LI.neq-1]);
    
    //==================================================================
    // hydrogen part
    //==================================================================
    fcn_HI_effective(z, rho, Xe, Tg, NH, H_z, Xp, iHI, LI);

    //==================================================================
    // neutral helium part
    //==================================================================
    if(flag_He==1) fcn_HeI_effective(z, Xe, Tg, NH, H_z, XHeII, iHeI, LI);
    
    
    //==================================================================
    // DM annihilation
    //==================================================================
    if(DM_annihilation)
//        evaluate_DM_annihilation_terms(z, parameters.fHe, Tg, Xp, XHeII, fDM_CosmoRec,
//                                       LI.g[iHI], LI.g[iHeI], LI.g[LI.neq-1]);
        
        evaluate_DM_decay_terms(z, parameters.fHe, Tg, Xp, XHeII, fDM_CosmoRec, 1.0,
                                LI.g[iHI], LI.g[iHeI], LI.g[LI.neq-1]);
    
    //==================================================================
    // heating by magnetic fields
    //==================================================================
    if(magnetic_fields)
        evaluate_magnetic_field_heating(z, parameters.fHe, Tg, Tg*rho,
                                        Xp, XHeII,
                                        B0_CosmoRec, nB_CosmoRec, kD_CosmoRec,
                                        LI.g[iHI], LI.g[LI.neq-1]);

    return;
}

//==================================================================================================
//
// jac for ODE evaluation
//
//==================================================================================================
double df_dx_2point_sym(const double &fp1, const double &fm1, const double &h, const double &eps)
{
    if(fm1==0.0){  cout << " alert!!! " << endl; return fp1/(2.0*h); }  
    
    double dum=fp1/fm1-1.0;
    if(fabs(dum)<=eps) return 0.0; 
    return fm1*dum/(2.0*h);
}

//==================================================================================================
void fcn_for_jacobian_eval(int jX, int neq, double DXj, double z, double Tg, 
                           double NH, Data_Level_I &LI0, double *f)
{
    double Te=Tg;   

    //==============================================================================================
    // reset everything 
    //==============================================================================================
    for(int i=0; i<LI0.neq; i++) f[i]=1.0e-300;        
    
    //==============================================================================================
    LI0.X[jX]+=DXj;
    
    //------------------------------------------------------------------
    // HI & HeI populations
    //------------------------------------------------------------------
    if(jX>=LI0.index_HI && jX<LI0.index_HeI) 
        Hydrogen_Atoms.Set_population_of_level(jX-LI0.index_HI, LI0.X[jX]);
    
    else if(jX>=LI0.index_HeI && jX<LI0.neq-1) 
        HeI_Atoms.Set_Xi(jX-LI0.index_HeI, LI0.X[jX]);       
    
    double Xe, Xp, XHeII, XHeIII;
    compute_fractions(LI0, Xe, Xp, XHeII, XHeIII);

    //------------------------------------------------------------------
    // rho-dependent parts
    //------------------------------------------------------------------
    if(jX==LI0.neq-1)
    {
        update_Solver_vars_eff(&LI0.X[LI0.index_HI], Hydrogen_Atoms,
                               &LI0.X[LI0.index_HeI], HeI_Atoms,
                               z, LI0.X[LI0.neq-1], 0);
        
        ODE_HI_effective::evaluate_effective_Rci_Ric_terms_for_Te_derivative(z, Xe, NH*Xp, Ai_df_dx, 
                                                                             Bi_df_dx, get_HI_index, 
                                                                             &f[LI0.index_HI]);
        
        ODE_HI_effective::evaluate_effective_Rij_terms(z, Tg, Hydrogen_Atoms, RijVec_df_dx, 
                                                       get_HI_index, &f[LI0.index_HI]);        
    }
    else 
    {
        //==========================================================================
        // Diffusion correction (added 18.07.2010)
        //==========================================================================
        if(Diffusion_correction_is_on)
        {
            evaluate_HI_Diffusion_correction(z, Hydrogen_Atoms, NH, cosmos.H(z), 
                                             Tg, &f[LI0.index_HI]);
            
            if(DI1_2s_correction_on) 
                evaluate_HI_DI1_DI2_correction(z, Hydrogen_Atoms.X(0), Hydrogen_Atoms.X(1),
                                               Hydrogen_Atoms.Level(1).Get_Dnu_1s(), 
                                               Tg, &f[LI0.index_HI]);
        }      
        
        if(flag_He==1)
        {
            //==================================================================
            // HeI-diffusion correction (added 22.01.2011; JF & JC)
            //==================================================================
            if(Diffusion_correction_HeI_is_on)
            {
                double dum;
                evaluate_HeI_Diffusion_correction(z, HeI_Atoms, NH, cosmos.H(z), Tg, 
                                                  &f[LI0.index_HeI], dum, f[LI0.index_HI]);
                
                if(DI1_2s_correction_on_HeI) 
                    evaluate_HeI_DI1_DI2_correction(z, HeI_Atoms.Xi(0), HeI_Atoms.Xi(1), 
                                                    HeI_Atoms.Sing.Level(2, 0).Get_Dnu_1s2(), 
                                                    Tg, &f[LI0.index_HeI]);
            }                

            //==========================================================================
            // HeI-feedback
            //==========================================================================
            if(HeISTfeedback==1 && !Diffusion_correction_HeI_is_on)
            {
                if(_HeI_feedback_w_HI_abs==0) 
                    evaluate_HeI_feedback(z, Hydrogen_Atoms.X(1, 0), HeI_Atoms, cosmos, HeI_Lines,
                                          Dn_HeI_feedback_no_HI_abs, &f[LI0.index_HeI]);
                
                if(_HeI_feedback_w_HI_abs==1) 
                    evaluate_HeI_feedback(z, Hydrogen_Atoms.X(1, 0), HeI_Atoms, cosmos, HeI_Lines,
                                          Dn_HeI_feedback, &f[LI0.index_HeI]);
            }
            
            //==================================================================
            // absorption by neutral hydrogen
            //==================================================================
            if(flag_HI_absorption==1 && z<=zcrit_HI*FC_CosmoRec.f_t && !Diffusion_correction_HeI_is_on)
            {  
                if(jX==LI0.index_HeI || jX==LI0.index_HeI+2 || jX==LI0.index_HI)
                {
                    double dXe_HeI_HI_abs_dt=0.0;
                    //==================================================================
                    // HeI Ly-a line
                    //==================================================================
                    evaluate_HI_abs_HeI(z, Hydrogen_Atoms.X(1, 0), NH, cosmos.H(z), Tg, Te,
                                        HeI_Atoms, cosmos, HeI_Atoms.HILyc, 
                                        DP_interpol_S, dXe_HeI_HI_abs_dt);
                    
                    // Helium levels
                    f[LI0.index_HeI+0]+= dXe_HeI_HI_abs_dt;
                    f[LI0.index_HeI+2]+=-dXe_HeI_HI_abs_dt;
                    // HI 1s
                    f[LI0.index_HI]+=-dXe_HeI_HI_abs_dt;
                }
                
                //==================================================================
                // from intercombination line
                //==================================================================
                if(jX==LI0.index_HeI || jX==LI0.index_HI || 
                   jX==LI0.index_HeI+(int)HeI_Atoms.Get_Level_index(2, 1, 1, 1))
                {
                    if(z<=zcrit_HI_Int*FC_CosmoRec.f_t && flag_spin_forbidden==1)
                    {
                        double dXe_HeI_HI_abs_dt=0.0;
                        evaluate_HI_abs_HeI_Intercombination(z, Hydrogen_Atoms.X(1, 0), NH, 
                                                             cosmos.H(z), Tg, Te, 
                                                             HeI_Atoms, cosmos, HeI_Atoms.HILyc, 
                                                             DP_interpol_T, dXe_HeI_HI_abs_dt);     
                        
                        // Helium levels
                        f[LI0.index_HeI+0]+= dXe_HeI_HI_abs_dt;
                        f[LI0.index_HeI+HeI_Atoms.Get_Level_index(2, 1, 1, 1)]+=-dXe_HeI_HI_abs_dt;
                        // HI 1s
                        f[LI0.index_HI]+=-dXe_HeI_HI_abs_dt;
                    }
                }
            }       
        }
    }
        
    //==============================================================================================
    // DM annihilation
    //==============================================================================================
    if(DM_annihilation && (jX==LI0.index_HI || jX==LI0.index_HeI || jX==LI0.neq-1))
    {
        //evaluate_DM_annihilation_terms(z, parameters.fHe, Tg, Xp, XHeII, fDM_CosmoRec,
        //                               f[LI0.index_HI], f[LI0.index_HeI], f[LI0.neq-1]);

        evaluate_DM_decay_terms(z, parameters.fHe, Tg, Xp, XHeII, fDM_CosmoRec, 1.0,
                                f[LI0.index_HI], f[LI0.index_HeI], f[LI0.neq-1]);
    }

    //==============================================================================================
    // heating by magnetic fields
    //==============================================================================================
    if(magnetic_fields && (jX==LI0.index_HI || jX==LI0.index_HeI || jX==LI0.neq-1))
    {
        evaluate_magnetic_field_heating(z, parameters.fHe, Tg, Tg*LI0.X[LI0.neq-1],
                                        Xp, XHeII,
                                        B0_CosmoRec, nB_CosmoRec, kD_CosmoRec,
                                        f[LI0.index_HI], f[LI0.neq-1]);
    }

    //==============================================================================================
    // restore old values
    //==============================================================================================
    LI0.X[jX]-=DXj;
    
    //------------------------------------------------------------------
    // RESET HI & HeI populations
    //------------------------------------------------------------------
    if(jX>=LI0.index_HI && jX<LI0.index_HeI) 
        Hydrogen_Atoms.Set_population_of_level(jX-LI0.index_HI, LI0.X[jX]);
    
    else if(flag_He==1 && jX>=LI0.index_HeI && jX<LI0.neq-1) 
        HeI_Atoms.Set_Xi(jX-LI0.index_HeI, LI0.X[jX]);     

    //===========================================================================
    // derivative with respect to rho --> restore recombination rates after calls
    //===========================================================================
    if(jX==LI0.neq-1)
    {
        update_Solver_vars_eff(&LI0.X[LI0.index_HI], Hydrogen_Atoms,
                               &LI0.X[LI0.index_HeI], HeI_Atoms,
                               z, LI0.X[LI0.neq-1], 0);
    }
    //==============================================================================================
    //==============================================================================================

    return;
}

//==================================================================================================
void evaluate_derives_atomic_part(double *r, double DXj, double eps, int iY, int neqY, int (*get_index_func)(int))
{
    r[iY]+=df_dx_2point_sym(Level_I_temp_JAC.g[iY], Level_I_temp_JAC.X[iY], DXj, eps);
    
    for(int m=0; m<neqY; m++)
    {
        int kl=iY+get_index_func(m);
        r[kl]+=df_dx_2point_sym(Level_I_temp_JAC.g[kl], Level_I_temp_JAC.X[kl], DXj, eps);
    }

    return;
}


void compute_numerical_derivative(int colX, int neq, int neq_HI, int neq_HeI, double z,
                                  double Tg, double NH, Data_Level_I &LI0, double *r)
{ 
    //===========================================================================
    // r[i] is reset here
    //===========================================================================   
    double X0=LI0.X[colX], DXj=X0*1.0e-7, eps=1.0e-12;
    if(X0==0.0) DXj=1.0e-7;
    
    //===========================================================================
    // numerical derivatives with respect to Xj: here part of the populations and 
    // also recombination rates could be changed
    //===========================================================================
    // get f(Xj+DXj)
    fcn_for_jacobian_eval(colX, neq,  DXj, z, Tg, NH, LI0, Level_I_temp_JAC.g);
    // get f(Xj-DXj)
    fcn_for_jacobian_eval(colX, neq, -DXj, z, Tg, NH, LI0, Level_I_temp_JAC.X);
    
    //===========================================================================
    // define numerical derivative
    //===========================================================================
    r[0]+=df_dx_2point_sym(Level_I_temp_JAC.g[0], Level_I_temp_JAC.X[0], DXj, eps);
    evaluate_derives_atomic_part(r, DXj, eps, LI0.index_HI, neq_HI, get_HI_index);
    if(flag_He==1) evaluate_derives_atomic_part(r, DXj, eps, LI0.index_HeI, neq_HeI, get_HeI_index);
    
    return;
}

//==================================================================================================
//
//
//==================================================================================================
static double jac_effective_z=-1e+300;

//==================================================================================================
//
// jac for ODE evaluation
//
//==================================================================================================
void jac_effective(double z, Data_Level_I &LI, int col)
{
    //==================================================================
    // do numerical derivative if it was set that way
    //==================================================================
    if(ODE_Solver_info_CR.num_jac) return fcn_effective(z, LI);

    //==================================================================
    // common data
    //==================================================================
    int iHI=LI.index_HI;
    int iHeI=LI.index_HeI;
    int neq_HI=get_number_of_resolved_levels();
    int neq_HeI=get_number_of_resolved_levels_HeI();
    int neq=1+neq_HI+neq_HeI;
    
    //==================================================================
    // first decide if one should update z-dependent functions
    //==================================================================
    if(jac_effective_z!=z)
    {
        jac_effective_z=z;
        
        //==================================================================
        // copy Xi_H --> Atoms && update recombination coefficients
        //==================================================================
        update_Solver_vars_eff(&LI.X[iHI], Hydrogen_Atoms,
                               &LI.X[iHeI], HeI_Atoms,
                               z, LI.X[LI.neq-1]);
    }   
    
    //==================================================================
    // local variables
    //==================================================================
    double rho=LI.X[LI.neq-1];
    double Tg=cosmos.TCMB(z);
    double Xe, Xp, XHeII, XHeIII;
    compute_fractions(LI, Xe, Xp, XHeII, XHeIII);
    double NH=cosmos.NH(z), H_z=cosmos.H(z);
    
    //==================================================================
    // reset everything 
    //==================================================================
    for(int i=0; i<LI.neq; i++) LI.g[i]=0.0;
    
    //==================================================================
    // Te derivatives
    //==================================================================
    if(col==0)
    {
        ODE_effective::df_drho_evaluate_TM(z, Xe, parameters.fHe, Tg, H_z, 
                                           LI.g[LI.neq-1]);
        
        compute_numerical_derivative(LI.neq-1, neq, neq_HI, neq_HeI, 
                                     z, Tg, NH, LI, LI.g);
    }
    //==================================================================
    // HI X1s derivatives
    //==================================================================
    else if(col==1)
    {
        ODE_effective::df_dXHI1s_evaluate_TM(z, Xe, parameters.fHe, rho, Tg, H_z, 
                                             LI.g[LI.neq-1]);
        
        ODE_effective::df_dX1s_evaluate_2s_two_photon_decay(z, Tg, 
                                                            Hydrogen_Atoms.Level(1).Get_Dnu_1s(), 
                                                            const_HI_A2s_1s, LI.g[iHI], 
                                                            LI.g[iHI+1]); 
        
        ODE_HI_effective::df_dXHI1s_evaluate_effective_Rci_Ric_terms(z, NH*Xe, NH*Xp, Ai_df_dx,
                                                                     Bi_df_dx, get_HI_index, 
                                                                     &LI.g[iHI]);
        
        if(flag_He==1) 
            ODE_HeI_effective::df_dXHI1s_evaluate_effective_Rci_Ric_terms(z, NH*XHeII, Ai_df_dx_HeI,
                                                                          Bi_df_dx_HeI, 
                                                                          get_HeI_index, 
                                                                          &LI.g[iHeI]);
        
        int nLy_max=(int)min(Lyn_max_parameter, 
                             Hydrogen_Atoms.Get_n_of_Level(get_HI_index(neq_HI-1)));
        
        for(int n=2; n<=nLy_max; n++)
        {
            Transition_Data_A T=Hydrogen_Atoms.Level(n, 1).Get_Data_A(1, 0);
            
            int nly=Hydrogen_Atoms.Get_Level_index(n, 1);
            ODE_effective::df_dX1s_evaluate_Ly_n_channel(z, Tg, Hydrogen_Atoms.X(0), 
                                                         Hydrogen_Atoms.X(n, 1), NH, H_z, 
                                                         T.A21, T.lambda21, T.Dnu, 
                                                         LI.g[iHI], 
                                                         LI.g[iHI+nly]);
        }
        
        compute_numerical_derivative(iHI, neq, neq_HI, neq_HeI, z, Tg, NH, LI, LI.g);
    }
    //==================================================================
    // HI Xi derivatives (i>1s)
    //==================================================================
    else if(col>1 && col<2+neq_HI)
    {
        int m=col-2, ik=get_HI_index(m);
        ///
        if(col==2) 
            ODE_effective::df_dX2s_evaluate_2s_two_photon_decay(z, Tg, 
                                                                Hydrogen_Atoms.Level(1).Get_Dnu_1s(), 
                                                                const_HI_A2s_1s, LI.g[iHI], 
                                                                LI.g[iHI+1]);
        
        ODE_HI_effective::df_dXi_evaluate_effective_Rci_Ric_terms(m, z, Ai_df_dx, Bi_df_dx, 
                                                                  get_HI_index, 
                                                                  &LI.g[iHI]);
        
        ODE_HI_effective::df_dXi_evaluate_effective_Rij_terms(m, z, Tg, Hydrogen_Atoms, 
                                                              RijVec_df_dx, get_HI_index, 
                                                              &LI.g[iHI]);
        
        int n=Hydrogen_Atoms.Get_n_of_Level(ik);
        int l=Hydrogen_Atoms.Get_l_of_Level(ik);
        
        if(n<=Lyn_max_parameter && l==1)
        {
            Transition_Data_A T=Hydrogen_Atoms.Level(ik).Get_Data_A(1, 0);
            
            ODE_effective::df_dXnp_evaluate_Ly_n_channel(z, Tg, Hydrogen_Atoms.X(0), 
                                                         Hydrogen_Atoms.X(ik), NH, H_z, T.A21, 
                                                         T.lambda21, T.Dnu, LI.g[iHI], 
                                                         LI.g[iHI+ik]);
            
            //==================================================================
            // nD-1s quadrupole lines (added May 2011)
            //==================================================================
            if(Atom_activate_HI_Quadrupole_lines && n>2 && n<=nD_Quadrupole_max)
                ODE_effective::df_dXnp_evaluate_nD_1s_Q_channel(Hydrogen_Atoms.X(n, 2), 
                                                                Hydrogen_Atoms.Level(n, 2).Get_A21(1, 0), 
                                                                LI.g[iHI+Hydrogen_Atoms.Get_Level_index(n, 1)], 
                                                                LI.g[iHI+Hydrogen_Atoms.Get_Level_index(n, 2)]);
        }
        
        //==================================================================
        // nD-1s quadrupole lines (added May 2011)
        //==================================================================
        if(Atom_activate_HI_Quadrupole_lines && n>2 && n<=nD_Quadrupole_max && l==2)
            ODE_effective::df_dXnd_evaluate_nD_1s_Q_channel(Hydrogen_Atoms.X(n, 1), 
                                                            Hydrogen_Atoms.Level(n, 2).Get_A21(1, 0), 
                                                            LI.g[iHI+Hydrogen_Atoms.Get_Level_index(n, 1)], 
                                                            LI.g[iHI+Hydrogen_Atoms.Get_Level_index(n, 2)]);
        
        
        compute_numerical_derivative(iHI+ik, neq, neq_HI, neq_HeI, z, Tg, NH, LI, LI.g);
    }
    //==================================================================
    // HeI X1s derivatives
    //==================================================================
    else if(flag_He==1 && col==2+neq_HI)
    {
        ODE_effective::df_dXHeI1s_evaluate_TM(z, Xe, parameters.fHe, rho, Tg, H_z, 
                                              LI.g[LI.neq-1]);
        
        ODE_effective::df_dX1s_evaluate_2s_two_photon_decay(z, Tg, 
                                                            HeI_Atoms.Sing.Level(2, 0).Get_Dnu_1s2(), 
                                                            const_HeI_A2s_1s, 
                                                            LI.g[iHeI], 
                                                            LI.g[iHeI+1]); 
        
        ODE_HI_effective ::df_dXHeI1s_evaluate_effective_Rci_Ric_terms(z, NH*Xp, Ai_df_dx, Bi_df_dx,
                                                                       get_HI_index, 
                                                                       &LI.g[iHI]);
        
        ODE_HeI_effective::df_dXHeI1s_evaluate_effective_Rci_Ric_terms(z, NH*Xe, NH*XHeII,
                                                                       Ai_df_dx_HeI, 
                                                                       Bi_df_dx_HeI, 
                                                                       get_HeI_index, 
                                                                       &LI.g[iHeI]);
        
        for(int k=1; k<neq_HeI; k++)
        {
            int ik=get_HeI_index(k);
            if(flag_spin_forbidden==0 && ik==(int)HeI_Atoms.Get_Level_index(2, 1, 1, 1)) continue;
            
            Transition_Data_HeI_A T=HeI_Atoms.Get_Trans_Data(ik, 1, 0, 0, 0);
            if(T.A21!=0) 
                ODE_effective::df_dX1s_evaluate_Ly_n_channel(z, Tg, HeI_Atoms.Xi(0), 
                                                             HeI_Atoms.Xi(ik), NH, H_z, T.A21, 
                                                             T.lambda21, T.Dnu, 
                                                             LI.g[iHeI], 
                                                             LI.g[iHeI+ik], 
                                                             HeI_Atoms.Get_gw(ik)/T.gwp);
        }
        
        compute_numerical_derivative(iHeI, neq, neq_HI, neq_HeI, z, Tg, NH, LI, LI.g);
    }
    //==================================================================
    // HeI Xi derivatives (i>1s)
    //==================================================================
    else if(flag_He==1 && col>2+neq_HI)
    {
        int m=col-(3+neq_HI), ik=get_HeI_index(m);
        
        //if(3+neq_HI) // JC: 23rd July 2015: this was a stupid bug (seemed to not make a difference...)
        if(col==3+neq_HI)
            ODE_effective::df_dX2s_evaluate_2s_two_photon_decay(z, Tg, 
                                                                HeI_Atoms.Sing.Level(2,0).Get_Dnu_1s2(), 
                                                                const_HeI_A2s_1s, 
                                                                LI.g[iHeI], 
                                                                LI.g[iHeI+1]); 
        
        ODE_HeI_effective::df_dXi_evaluate_effective_Rci_Ric_terms(m, z, Ai_df_dx_HeI, 
                                                                   Bi_df_dx_HeI, get_HeI_index, 
                                                                   &LI.g[iHeI]);
        
        ODE_HeI_effective::df_dXi_evaluate_effective_Rij_terms(m, z, Tg, HeI_Atoms, 
                                                               RijVec_df_dx_HeI, 
                                                               get_HeI_index, 
                                                               &LI.g[iHeI]);        
        
        // JC: 23rd July 2015: l=0 and J=1 transition was neglected with previous version.
        //                     Very small rate though so this did not have a large effect.
        if(!( flag_spin_forbidden==0 && ik==(int)HeI_Atoms.Get_Level_index(2, 1, 1, 1) ) )
        {
            Transition_Data_HeI_A T=HeI_Atoms.Get_Trans_Data(ik, 1, 0, 0, 0);
            if(T.A21!=0)
                ODE_effective::df_dXnp_evaluate_Ly_n_channel(z, Tg, HeI_Atoms.Xi(0),
                                                             HeI_Atoms.Xi(ik), NH, H_z,
                                                             T.A21, T.lambda21, T.Dnu,
                                                             LI.g[iHeI],
                                                             LI.g[iHeI+ik],
                                                             HeI_Atoms.Get_gw(ik)/T.gwp);
        }
        
        compute_numerical_derivative(iHeI+ik, neq, neq_HI, neq_HeI, z, Tg, NH, LI, LI.g);
    }
    else if(flag_He==0 && col>=2+neq_HI){}
    else{ cerr << " This should not happen. Exiting. " << endl; exit(0); }
    
    return;
}

//==================================================================================================
void fcn_effective(int *neq, double *z, double *y, double *f, int col)
{    
    int index_HI=Level_I_temp.index_HI;
    int index_HeI=Level_I_temp.index_HeI;
    int nres_HI=get_number_of_resolved_levels();
    int nres_HeI=get_number_of_resolved_levels_HeI();
        
    //==================================================================
    copy_ysol_to_LI(y, Level_I_temp, Level_I_temp.neq, index_HI, index_HeI);
    
    //==================================================================
    // col<0 : evaluate the rhs of system  dy/dz == g(t,y); 
    // col>=0: evaluate derivative with respect to j==col
    //==================================================================
    if(col<0) fcn_effective(*z, Level_I_temp); 
    else jac_effective(*z, Level_I_temp, col); 
    
    //==================================================================
    // remember that Level_I_temp.g is g(t, y) --> need dz/dt factor!!!
    //==================================================================
    double dz_dt=-cosmos.H( (*z) )*( 1.0+ (*z) );
    
    //==================================================================
    // reset everything
    //==================================================================
    for(int i=0; i<*neq; i++) f[i]=0.0;

    //==================================================================
    f[0]=Level_I_temp.g[Level_I_temp.neq-1]/dz_dt;
    f[1]=Level_I_temp.g[index_HI]/dz_dt;
    for(int k=0; k<nres_HI; k++) f[k+2]=Level_I_temp.g[index_HI+get_HI_index(k)]/dz_dt;
    
    if(flag_He==1)
    {
        f[nres_HI+2]=Level_I_temp.g[index_HeI]/dz_dt;
        
        for(int k=0; k<nres_HeI; k++) 
            f[nres_HI+k+3]=Level_I_temp.g[index_HeI+get_HeI_index(k)]/dz_dt;
    }
    
    return;
}   

//==================================================================================================
//==================================================================================================
