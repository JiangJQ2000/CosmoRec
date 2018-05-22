//===========================================================================================================
// Author Jens Chluba 
// here functions which are common for all the programs
// connected with the recombination problem are defined
//===========================================================================================================

//===========================================================================================================
void set_cosmology(Cosmos::parameters &cosm, Cosmos &cosmos)
{
    if(parameters.verbosity>=1) print_message(" Initializing cosmology object ");

    //=======================================================================================================
    // initializing cosmos object...
    //=======================================================================================================
    cosmos.init(cosm);
    cosmos.init_splines();
    cosmos.set_sigT_function(cosm.Rec.aS, cosm.Rec.mS, cosm.Rec.pS); // 16.05.2017: JC
    
    if(parameters.verbosity>=1) cosmos.display_cosmology();
    
    return;
}

//===========================================================================================================
// initialization
//===========================================================================================================
void set_array_dimensions_and_parameters(int wait=0)
{
    //=======================================================================================================
    // setup Level I and start indices of HI and HeI
    //=======================================================================================================
    Level_I.nHIeq=Hydrogen_Atoms.Get_total_number_of_Levels();
    Level_I.nHeIeq=HeI_Atoms.Get_total_number_of_Levels();
    Level_I.neq=Level_I.nHIeq+Level_I.nHeIeq+2;
    Level_I.index_HI=1;
    Level_I.index_HeI=Level_I.nHIeq+1;

    //-------------------------------------------------------------------------------------------------------
    // copy info to other arrays
    //-------------------------------------------------------------------------------------------------------
    Level_I_temp.nHIeq=Level_I_temp_JAC.nHIeq=Level_I.nHIeq;
    Level_I_temp.nHeIeq=Level_I_temp_JAC.nHeIeq=Level_I.nHeIeq;
    Level_I_temp.neq=Level_I_temp_JAC.neq=Level_I.neq;
    Level_I_temp.index_HI=Level_I_temp_JAC.index_HI=Level_I.index_HI;
    Level_I_temp.index_HeI=Level_I_temp_JAC.index_HeI=Level_I.index_HeI;

    //=======================================================================================================
    if(parameters.verbosity>=2)
    {
        cout << "\n total number of equations: " << Level_I.neq
             << "\n Hydrogen equations: " << Level_I.nHIeq
             << "\n Helium equations: " << Level_I.nHeIeq << endl;
    
        cout << " nShells: " << parameters.CR.nShells
             << "\n included number of hydrogen levels: " << Level_I.nHIeq << endl;
    
        cout << " nShellsHeI: " << parameters.CR.nShellsHeI
             << "\n included number of neutral helium levels: " << Level_I.nHeIeq << endl;
    }
    //=======================================================================================================

    if(wait==1) wait_f_r();
    
    return;
}

//===========================================================================================================
void read_cosmology_params_from_file(file_content &pfc, Cosmos::parameters &C, bool show_params)
{
    bool found; double dval;

    //==============================================================================================
    // cosmological parameters
    //==============================================================================================
    parser_read(pfc, "Yp", dval, found, show_params);
    if(found) C.Yp=dval;
    
    parser_read(pfc, "T0", dval, found, show_params);
    if(found) C.T0=dval;
    
    parser_read(pfc, "Omega_m", dval, found, show_params);
    if(found) C.Omega_m=dval;
    
    parser_read(pfc, "Omega_b", dval, found, show_params);
    if(found) C.Omega_b=dval;
    
    parser_read(pfc, "Omega_k", dval, found, show_params);
    if(found) C.Omega_k=dval;
    
    parser_read(pfc, "Omega_L", dval, found, show_params);
    if(found) C.Omega_L=dval;
    
    parser_read(pfc, "h100", dval, found, show_params);
    if(found) C.h100=dval;
    
    parser_read(pfc, "N_eff", dval, found, show_params);
    if(found) C.Nnu=dval;
    
    return;
}

//===========================================================================================================
void read_CosmoRec_params_from_file(file_content &pfc, CosmoRec_parameters &CR, bool show_params)
{
    bool found; int ival;
    
    //----------------------------------------------------------------------------------------------
    // settings for HI part
    //----------------------------------------------------------------------------------------------
    parser_read(pfc, "number of hydrogen shells", ival, found, show_params);
    if(found) CR.nShells=ival;

    parser_read(pfc, "nS for effective HI rates", ival, found, show_params);
    if(found) CR.nS_effective=ival;

    parser_read(pfc, "run HI PDE part", ival, found, show_params);
    if(found) CR.Diffusion_flag=ival;

    parser_read(pfc, "correction to 2s-1s channel", ival, found, show_params);
    if(found) CR.induced_flag=ival;

    parser_read(pfc, "nS two-photon decays", ival, found, show_params);
    if(found) CR.nS_2gamma=ival;

    parser_read(pfc, "nS Raman-scattering", ival, found, show_params);
    if(found) CR.nS_Raman=ival;

    //----------------------------------------------------------------------------------------------
    // settings for HeI part
    //----------------------------------------------------------------------------------------------
    parser_read(pfc, "number of helium shells", ival, found, show_params);
    if(found) CR.nShellsHeI=ival;

    parser_read(pfc, "HI absorption", ival, found, show_params);
    if(found) CR.HI_absorption=ival;

    parser_read(pfc, "spin forbidden transitions", ival, found, show_params);
    if(found) CR.spin_forbidden=ival;

    parser_read(pfc, "HeI feedback setting", ival, found, show_params);
    if(found) CR.HeI_Feedback=ival;

    return;
}

//===========================================================================================================
void read_Recfast_pp_params_from_file(file_content &pfc, RecPars &R, bool show_params)
{
    bool found; int ival; double dval;
    
    //----------------------------------------------------------------------------------------------
    // standard recombination physics settings
    // redshift range and number of points will be set by cosmology object, so don't bother
    //----------------------------------------------------------------------------------------------
    parser_read(pfc, "Recfast fudge factor", dval, found, show_params);
    if(found) R.F=dval;
    
    parser_read(pfc, "include correction function", ival, found, show_params);
    if(found) R.include_CF=(bool)ival;
    
    parser_read(pfc, "A2s1s", dval, found, show_params);
    if(found) R.A2s1s=dval;
    
    //==============================================================================================
    // different non-standard cases
    //----------------------------------------------------------------------------------------------
    // annihilation
    //----------------------------------------------------------------------------------------------
    parser_read(pfc, "f_ann", dval, found, show_params);
    if(found) R.f_ann=dval;

    //----------------------------------------------------------------------------------------------
    // decay
    //----------------------------------------------------------------------------------------------
    parser_read(pfc, "f_dec", dval, found, show_params);
    if(found) R.f_dec=dval;

    if(R.f_dec>0) {
        parser_read(pfc, "Gamma_dec", dval, found, show_params);
        if(found) R.Gamma_dec=dval;
    }
    
    //----------------------------------------------------------------------------------------------
    // magnetic field heating
    //----------------------------------------------------------------------------------------------
    parser_read(pfc, "B0", dval, found, show_params);
    if(found) R.B0=dval;

    if(R.B0>0.0) {
        parser_read(pfc, "include turbulent decay", ival, found, show_params);
        if(found) R.turb=(bool)ival;
        
        parser_read(pfc, "include ambipolar diffusion", ival, found, show_params);
        if(found) R.ambi=(bool)ival;
        
        parser_read(pfc, "nB", dval, found, show_params);
        if(found) R.nB=dval;

        parser_read(pfc, "Paoletti-Finelli Lorentz", ival, found, show_params);
        if(found) R.KK_or_PF_Lorentz=(bool)ival;
    }
    
    //----------------------------------------------------------------------------------------------
    // variation of alpha and me/mp
    //----------------------------------------------------------------------------------------------
    parser_read(pfc, "Variation mode", ival, found, show_params);
    if(found) R.RM=ival;

    if(R.RM>0) {
        parser_read(pfc, "alp/alp_ref", dval, found, show_params);
        if(found) R.aS=dval;
        
        parser_read(pfc, "me/me_ref", dval, found, show_params);
        if(found) R.mS=dval;

        parser_read(pfc, "power for (1+z)^p", dval, found, show_params);
        if(found) R.pS=dval;
    }
    
    return;
}

//===========================================================================================================
void read_entries_from_parameter_file(string filename, Parameters_form_initfile &inp)
{
    //==============================================================================================
    // read values (if given) from file
    //==============================================================================================
    string funcname="read_entries_from_parameter_file";
    bool found; int ival; double dval; string sval;
    file_content pfc;
    parser_read_file(filename, pfc, 0);
    
    //==============================================================================================
    // reset verbosity level
    //==============================================================================================
    parser_read(pfc, "verbosity level CosmoRec", ival, found, 0);
    if(!found) parser_read(pfc, "verbosity level Recfast++", ival, found, 0);
    if(found) inp.verbosity=ival;

    if(inp.verbosity>-1)
        cout << " CosmoRec :: " << funcname << " : reading parameters from: " << filename << endl;
    
    bool show_params=(inp.verbosity > 3 ? 1 : 0);
    // read again to show the parameter...
    if(show_params) parser_read(pfc, "verbosity level CosmoRec", ival, found, show_params);

    //==============================================================================================
    // path etc
    //==============================================================================================
    parser_read(pfc, "path for output", sval, found, show_params);
    if(found) inp.path=sval;
    
    parser_read(pfc, "addition to name for output", sval, found, show_params);
    if(found) inp.addname=sval;
    
    //----------------------------------------------------------------------------------------------
    // info for redshift output
    //----------------------------------------------------------------------------------------------
    parser_read(pfc, "nzpts", ival, found, show_params);
    if(found) inp.nz=ival;
    if(inp.nz<1000){ cerr << funcname << " : Choose nzpts > 1000 " << endl; exit(1); }
    
    parser_read(pfc, "zstart", dval, found, show_params);
    if(found) inp.zstart=dval;
    if(inp.zstart>3400){ cerr << funcname << " : Choose zstart <= 3400 " << endl; exit(1); }
    
    parser_read(pfc, "zend", dval, found, show_params);
    if(found) inp.zend=dval;
    if(inp.zend<0){ cerr << funcname << " : Choose zstart >=0 " << endl; exit(1); }

    //==============================================================================================
    // cosmological parameters
    //==============================================================================================
    read_cosmology_params_from_file(pfc, inp.cosm, show_params);

    //==============================================================================================
    // recombination parameters for CosmoRec
    //==============================================================================================
    read_CosmoRec_params_from_file(pfc, inp.CR, show_params);

    //==============================================================================================
    // recombination parameters for Recfast++
    //==============================================================================================
    read_Recfast_pp_params_from_file(pfc, inp.cosm.Rec, show_params);
    inp.cosm.Rec.verbosity=inp.verbosity;

    //==============================================================================================
    // cleanup
    //==============================================================================================
    parser_free(pfc);
    
    //==============================================================================================
    // show read parameters
    //==============================================================================================
    if(inp.verbosity>3) inp.show();
    
    return;
}

//===========================================================================================================
void check_and_set_cosmology_params(Cosmos::parameters &C, int verbosity)
{
    string funcname="check_and_set_cosmology_params";
    
    //==============================================================================================
    // pure cosmology part of parameters (T0, Yp, etc) is checked by cosmology object. but do set
    //==============================================================================================
    if(C.Omega_L<=0.0) C.Omega_L=1.0-C.Omega_m-C.Omega_k-cosmos.calc_Orel(C.T0, C.Nnu, C.h100);

    //----------------------------------------------------------------------------------------------
    // check recombination part of Recfast++
    //----------------------------------------------------------------------------------------------
    if(C.Rec.F==0) C.Rec.F=1.14;
    if(C.Rec.F<0){ cerr << funcname << " : Choose fudge >=0 " << endl; exit(1); }
    
    if(C.Rec.A2s1s<=0.0) C.Rec.A2s1s=const_HI_A2s_1s; // use values consistent with CosmoRec default
    // explicitly set value of A2s1s rate in CosmoRec (added 23.07.2014)
    if(C.Rec.A2s1s>0.0) const_HI_A2s_1s=C.Rec.A2s1s;

    //----------------------------------------------------------------------------------------------
    // annihilation
    //----------------------------------------------------------------------------------------------
    if(C.Rec.f_ann<0){ cerr << funcname << " : f_ann >=0 " << endl; exit(1); }
    
    //----------------------------------------------------------------------------------------------
    // decay
    //----------------------------------------------------------------------------------------------
    if(C.Rec.f_dec<0){ cerr << funcname << " : f_dec >=0 " << endl; exit(1); }
    if(C.Rec.Gamma_dec<0){ cerr << funcname << " : Gamma_dec >=0 " << endl; exit(1); }

    //----------------------------------------------------------------------------------------------
    // magnetic field heating
    //----------------------------------------------------------------------------------------------
    if(C.Rec.B0<0){ cerr << funcname << " : B0 >=0 " << endl; exit(1); }
    if(C.Rec.nB<-3.0){ cerr << funcname << " : nB > -3.0 " << endl; exit(1); }
   
    set_PMF_parameters(C.Rec.B0, C.Rec.nB, C.Rec.turb, C.Rec.ambi, C.Rec.KK_or_PF_Lorentz, C.h100);

    //----------------------------------------------------------------------------------------------
    // variation of alpha and me/mp
    //----------------------------------------------------------------------------------------------
    if(C.Rec.RM>6){ cerr << funcname << " : Variation mode <=6 " << endl; exit(1); }

    Set_VFC_pars(C.Rec.aS, C.Rec.mS, C.Rec.pS, C.Rec.RM);
    
    //==============================================================================================
    // show final parameters
    //==============================================================================================
    if(verbosity>2) C.show();

    return;
}

//===========================================================================================================
//
// setting the startup parameters according to param[12]
//
//===========================================================================================================
void set_startup_data(string filename,
                      Parameters_form_initfile &inp,
                      Cosmos &cos, int wait=0)
{
    if(inp.verbosity>=0) cout << endl;
    
    //=======================================================================================================
    // cosmology and Recfast++ parameter checks. This function also initializes derived CosmoRec parameters
    // from standard recombination settings.
    //=======================================================================================================
    check_and_set_cosmology_params(inp.cosm, inp.verbosity);
    
    //=======================================================================================================
    // CosmoRec specific settings
    //=======================================================================================================
    if(inp.CR.HeI_Feedback==0 || inp.CR.HeI_Feedback==1 || inp.CR.HeI_Feedback==-1)
    { 
        inp.CR.HeI_Feedback=0; 
        HeISTfeedback=0; 
        _HeI_feedback_w_HI_abs=0; 
    }
    else if(inp.CR.HeI_Feedback>0)
    { 
        inp.CR.HeI_Feedback=(int)min(inp.CR.HeI_Feedback, inp.CR.nShellsHeI);
        HeISTfeedback=1; 
        _HeI_feedback_w_HI_abs=0; 
    }
    else 
    { 
        inp.CR.HeI_Feedback=(int)min(-inp.CR.HeI_Feedback, inp.CR.nShellsHeI);
        HeISTfeedback=1;
        _HeI_feedback_w_HI_abs=1; 
    }
    
    if(inp.CR.HI_absorption==1){}
    else if(inp.CR.HI_absorption==2){ _HI_abs_appr_flag=1; inp.CR.HI_absorption=1; }
    else if(inp.CR.HI_absorption==3){ _HI_abs_appr_flag=1; inp.CR.HI_absorption=1; Diffusion_correction_HeI=1; }
    else if(inp.CR.HI_absorption>=4)
    { 
        cerr << " set_startup_data :: HI_absorption flag N/A. Exiting... " << endl;
        exit(0);
    }

    if(inp.cosm.Rec.f_ann>0) DM_annihilation=1;
    if(inp.cosm.Rec.f_dec>0) DM_decay=1;
    if(B0_CosmoRec >0) magnetic_fields=1;
    
    if(inp.verbosity>=1)
    {
        if(filename!="") cout << " Using parameters corresponding to " << filename << endl;
        cout << "\n zs: " << inp.zstart << "\t ze: " << inp.zend << "\t nzpts: " << inp.nz << endl;
        cout << " Y_p: " << inp.cosm.Yp << "\t TCMB: " << inp.cosm.T0 << endl;
        cout << " Omega_m: " << inp.cosm.Omega_m << "\t Omega_b: " << inp.cosm.Omega_b
             << "\n Omega_L: " << inp.cosm.Omega_L << "\t Omega_k: " << inp.cosm.Omega_k
             << "\t Omega_rel: " << cos.calc_Orel(inp.cosm.T0, inp.cosm.Nnu, inp.cosm.h100)
             << " Nnu= " << inp.cosm.Nnu << endl;
        cout << " Hubble constant in units of 100 km s-1 Mpc-1: " << inp.cosm.h100 << endl;
        cout << " Fudge factor for H recombination: " << inp.cosm.Rec.F << endl;
        cout << " HI A2s1s rate: " << const_HI_A2s_1s << endl;
    }
    
    if(DM_annihilation && inp.verbosity>=0)
        cout << " Including dark matter annihilation. f_DM= " << inp.cosm.Rec.f_ann << " eV/sec " << endl;

    if(DM_decay && inp.verbosity>=0)
        cout << " Including dark matter decay. f_dec= " << inp.cosm.Rec.f_dec
             << " Gamma_dec= " << inp.cosm.Rec.Gamma_dec << endl;

    if(magnetic_fields && inp.verbosity>=0)
        cout << " Including heating by magnetic fields. B0= " << B0_CosmoRec
             << " nG and nB= " << nB_CosmoRec << endl;

    if(FC_CosmoRec.runmode>0 && inp.verbosity>=0)
        cout << " Including variations of fundamental constants."
             << " Mode = " << FC_CosmoRec.runmode
             << ", a_scale = " << FC_CosmoRec.a_scale
             << ", m_scale = " << FC_CosmoRec.m_scale
             << ", p_scale = " << FC_CosmoRec.p_scale << endl;
    
    if(inp.CR.HI_absorption==1 && inp.verbosity>=0)
            cout << " HI absorption will be included during HeI-recombination " << endl;
    //
    if(Diffusion_correction_HeI==1 && inp.verbosity>=0)
        cout << " After the first iteration of the normal EMLA, the HeI-diffusion" 
             << " correction will be computed from PDE-solver " << endl;
    //
    if(inp.CR.spin_forbidden>=1)
    { 
        if(inp.verbosity>=0)
            cout << " HeI 2P-1S intercombination line is switched on " << endl; 
        inp.CR.spin_forbidden=1;
    }
    else
    { 
        if(inp.verbosity>=0)
            cout << " HeI 2P-1S intercombination line is switched off " << endl; 
        inp.CR.spin_forbidden=0;
    }

    if(Atom_activate_HI_Quadrupole_lines)
    { 
        if(inp.verbosity>=0 && inp.CR.nShells==3)
            cout << " HI quadrupole lines are included " << endl; 
    }
    
    if(inp.CR.Diffusion_flag>0)
    { 
        inp.CR.Diffusion_flag=1;
        Diffusion_correction=1;
        DI1_2s_correction_on=0;
        switch_off_1s_2s_correction();
        
        if(inp.verbosity>=0)
        {
            cout << " Diffusion correction for HI Ly-a will be included." 
                 << " The number of iterations will be "+int_to_string(Diff_iteration_max) 
                 << endl;
        }

        if(inp.CR.induced_flag>=1)
        { 
            switch_on_1s_2s_correction();
            DI1_2s_correction_on=1;
            
            if(inp.verbosity>=0)
                cout << " Switching on stimulated 2s-1s two-photon decays " << endl; 
            
            if(inp.CR.induced_flag==2)
            { 
                if(inp.verbosity>=0)
                    cout << " Ly-n feedback to the 1s-2s two-photon channel " << endl; 
            }
            else switch_off_1s_2s_absorption_correction();     // defined in Solve_PDEs.h
        }
    }
        
    //=======================================================================================================
    // below given z extend the solution not solving the whole system
    //=======================================================================================================
    if(inp.zend<z_extend)
    {
        flag_compute_Recfast_part=1;
        ze_Recfast=max(inp.zend, 0.001);
        inp.zend=z_extend;
    }   
    
    if( inp.verbosity>=0 && (write_Xe_Te_sol!=0 || write_populations!=0) )
    {
        cout << "\n output-path: " << inp.path << endl;
        cout << " addname    : " << inp.addname << endl;
    }
    
    set_detailed_name(inp, detailedname);
    
    if(inp.verbosity>=0) cout << endl;
    
    if(wait==1) wait_f_r();
    
    //=======================================================================================================
    // for CosmoRec run, use Recfast++ history without any non-standard physics
    //=======================================================================================================
    Cosmos::parameters cosm_ref;
    cosm_ref.h100=inp.cosm.h100; cosm_ref.T0=inp.cosm.T0;
    cosm_ref.Yp=inp.cosm.Yp; cosm_ref.Nnu=inp.cosm.Nnu;
    cosm_ref.Omega_m=inp.cosm.Omega_m; cosm_ref.Omega_b=inp.cosm.Omega_b;
    cosm_ref.Omega_L=inp.cosm.Omega_L; cosm_ref.Omega_k=inp.cosm.Omega_k;
    cosm_ref.Rec.F=inp.cosm.Rec.F;
    cosm_ref.Rec.A2s1s=inp.cosm.Rec.A2s1s;
    cosm_ref.Rec.include_CF=0;

    set_cosmology(cosm_ref, cos);

    //=======================================================================================================
    // set up atomic species
    //=======================================================================================================
    if(!hydrogen_and_helium_atom_are_set)
    { 
        setup_Hydrogen_atom();
        setup_Helium_atom();
        hydrogen_and_helium_atom_are_set=1;
    }

    return;
}

//===========================================================================================================
// reading startup parameters from startup file
//===========================================================================================================
void read_startup_data(string filename,
                       Parameters_form_initfile &inp,
                       Cosmos &cos, int wait=0)
{
    read_entries_from_parameter_file(filename, inp);
    
    set_startup_data(filename, inp, cos, wait);
    
    return;
}   

//===========================================================================================================
// reset to fresh startup
//===========================================================================================================
void restart_CosmoRec()
{   
    flag_He=1;
    set_array_dimensions_and_parameters();
    
    return;
}   

//===========================================================================================================
// memory (de)allocation
//===========================================================================================================
void clear_atoms()
{
    delete HAtoms;
    delete HeIAtoms;
    
    Ai_df_dx.clear();
    Bi_df_dx.clear();
    RijVec_df_dx.clear();
    //
    Ai_df_dx_HeI.clear();
    Bi_df_dx_HeI.clear();
    RijVec_df_dx_HeI.clear();
    
    output_CosmoRec.clear();
    pass_on_the_Solution_CosmoRec.clear();
    pass_on_the_Solution_CosmoRec_HeI.clear(); // JF
    
    return;
}

//===========================================================================================================
void allocate_memory()
{
    zarr =new double[parameters.nz];
    // Level I
    Level_I.X =new double[Level_I.neq];
    Level_I.g =new double[Level_I.neq];
    //
    Level_I_temp.X =new double[Level_I.neq];
    Level_I_temp.g =new double[Level_I.neq];
    //
    Level_I_temp_JAC.X =new double[Level_I.neq];
    Level_I_temp_JAC.g =new double[Level_I.neq];
    
    ODE_Solver_info_CR.Snewptr=NULL; // to make sure that the memory will be set

    return;
}

//===========================================================================================================
void deallocate_memory()
{
    delete[] zarr;
    // Level I
    delete[] Level_I.X;
    delete[] Level_I.g;
    //
    delete[] Level_I_temp.X;
    delete[] Level_I_temp.g;
    //
    delete[] Level_I_temp_JAC.X;
    delete[] Level_I_temp_JAC.g;
    
    delete[] Xfiles;
    
    clear_all_memory(Sz, tols, ODE_Solver_info_CR);
    
    return;
}

//===========================================================================================================
//
// to cleanup after finishing with CosmoRec
//
//===========================================================================================================
void cleanup_CosmoRec()
{
    //======================================================================
    // clean up memory
    //======================================================================
    deallocate_memory();
    clear_atoms();
    free_all_splines_JC();
    
    return;
}

//===========================================================================================================
//===========================================================================================================
