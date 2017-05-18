//--------------------------------------------------------------------------------------------
void set_detailed_name()
{
    detailedname="CosmoRec.";
    detailedname+="HI.nS_eff_"+int_to_string(nS_effective)+"."; 
    if(Atom_activate_HI_Quadrupole_lines) detailedname+="HI_Q.";
    //
    detailedname+="nH_"+int_to_string(nShells)+".";
    detailedname+="nHe_"+int_to_string(nShellsHeI)+".";
    //
    if(flag_spin_forbidden==0) detailedname+="no_Intercomb.";
    if(set_HeI_Quadrupole_lines) detailedname+="HeI_Q.";
    //
    if(Diffusion_correction_HeI) detailedname+="HeI_diff.";
    else
    {
        if(flag_HI_absorption==1 && _HI_abs_appr_flag==0) detailedname+="H_abs.";
        if(flag_HI_absorption==1 && _HI_abs_appr_flag==1) detailedname+="H_abs.fcorr.";
        //
        if(HeISTfeedback==1) detailedname+="HeI_feedback_nf_"+int_to_string(nHeFeedback)+".";
        if(_HeI_feedback_w_HI_abs==1) detailedname+="w_HI_abs.";
    }
    //
    if(DM_annihilation) detailedname+="DM_annihil.";

    if(magnetic_fields)
    {
        detailedname+="magn_heat.";
        if(decaying_turb_on) detailedname+="dec.";
        if(ambipolar_heat_on) detailedname+="ambi.";
        if(KK_or_PF_Lorentz) detailedname+="PF.";
    }
    
    if(FC_CosmoRec.runmode>0)
    {
        if(FC_CosmoRec.a_scale!=1.0) detailedname+="alp.";
        if(FC_CosmoRec.m_scale!=1.0) detailedname+="me.";
        if(FC_CosmoRec.a_scale!=1.0 || FC_CosmoRec.m_scale!=1.0)
        {
            if(FC_CosmoRec.runmode==1) detailedname+="Boltz.";
            else if(FC_CosmoRec.runmode==2) detailedname+="sigT.";
            else if(FC_CosmoRec.runmode==3) detailedname+="A2s1s.";
            else if(FC_CosmoRec.runmode==4) detailedname+="phot.";
            else if(FC_CosmoRec.runmode==5) detailedname+="Pesc.";  // currently not included...
            else if(FC_CosmoRec.runmode==6) detailedname+="all.";
        }
    }
    
    if(show_CosmoRec_mess>=0)
        cout << " set_detailed_name: detailedname: " << detailedname << endl;
    
    return; 
}

//--------------------------------------------------------------------------------------------
void open_files(ofstream *Xfiles, string path, int fresh=1)
{
    ios::openmode wrmode;
    if(fresh==1) wrmode=ios::out;
    else wrmode=ios::app;
    
    string prename=path+detailedname;
    if(write_Xe_Te_sol)
    {
        if(iteration_count==0) 
        {
            Xfiles[0].open((prename+"X_Recfast"+addname).c_str(), wrmode);
            Xfiles[0] << "# Recfast++ output. The columns are z, Xe=Ne/NH, Te\n#" << endl;
        }
        
        Xfiles[1].open((prename+"Xe"+addname).c_str(), wrmode);
        Xfiles[1] << "# CosmoRec "+CosmoRec_version+" output. The columns are z, Xe=Ne/NH, Te\n#" 
                  << endl;        
    }
    
    //  Xfiles[2].open((prename+"rho"+addname).c_str(), wrmode); 
    //
    if(write_populations) Xfiles[3].open((prename+"HI.X"+addname).c_str(), wrmode); 
//  Xfiles[4].open((prename+"HeI.X"+addname).c_str(), wrmode); 
//  Xfiles[5].open((prename+"HI.HeI.1s"+addname).c_str(), wrmode); 
    if(write_HI_distortion) Xfiles[5].open((prename+"HI_distortion"+addname).c_str(), wrmode); 

    // Added for correction functions LH & JC
    if(write_corrections)
    {
        if(Get_T_eval())
            Xfiles[6].open((prename+"correction_function.Tr"+addname).c_str(),wrmode);
        else Xfiles[6].open((prename+"correction_function"+addname).c_str(),wrmode);
    }
    
    // Added extra file option for separated Xe components
    if (write_split) Xfiles[7].open((prename+"Xe.H.HeI"+addname).c_str(),wrmode);

    if(write_Xe_Te_sol && iteration_count==0) Xfiles[0].precision(8);
    if(write_Xe_Te_sol) Xfiles[1].precision(8);
    if(write_populations) Xfiles[3].precision(12);
    if(write_HI_distortion) Xfiles[5].precision(12);
    if(write_corrections) Xfiles[6].precision(12); // Added for correction functions LH
    if(write_split) Xfiles[7].precision(12); // Added for split Xe LH

    
    return;
}

void close_files(ofstream *Xfiles)
{
    if(write_Xe_Te_sol && iteration_count==0) Xfiles[0].close();
    if(write_Xe_Te_sol) Xfiles[1].close();
    if(write_populations) Xfiles[3].close();
    if(write_HI_distortion) Xfiles[5].close();
    if(write_corrections) Xfiles[6].close(); // Added for correction functions LH
    if(write_split) Xfiles[7].close(); // Added for split Xe LH

    
    return;
}

//--------------------------------------------------------------------------------------------
void write_HI_HeI_1s_array_to_file(ofstream &ofile, double zval)
{
    ofile << scientific << zval;
    ofile << scientific << " " << Hydrogen_Atoms.X(1, 0);
    ofile << scientific << " " << HeI_Atoms.Xi(0) << endl;
    
    return;
}

void write_HI_X_array_to_file(ofstream &ofile, double zval, double *X, int nHIeq, int neq)
{
    ofile << scientific << zval;
    // Xe & HI Xi
    for(int l=0; l<Level_I.nHIeq+1; l++) ofile << scientific << " " << X[l];
    // rho
    ofile << scientific << " " << X[neq-1];
    // Xp, XeHeII
    ofile << scientific << " " << 1.0-Hydrogen_Atoms.X_tot();
    ofile << scientific << " " << parameters.fHe-HeI_Atoms.X_tot() << endl;
    
    return;
}

void write_HI_HeI_Xe_to_file(ofstream &ofile, double zval, double Xe)
{
    ofile << scientific << zval;
    ofile << scientific << " " << 1.0-Hydrogen_Atoms.X_tot();
    ofile << scientific << " " << parameters.fHe-HeI_Atoms.X_tot();
    ofile << scientific << " " << 1.0-Hydrogen_Atoms.X_tot()+parameters.fHe-HeI_Atoms.X_tot();
    ofile << scientific << " " << Xe << endl;
    return;
}

void write_correction_functions_to_file(ofstream &ofile, double zval, Cosmos &cosmos)
{
    // This function calculates and then writes the correction functions from CosmoRec and 
    // Recfast++ into their own files LH
    // Following the same schematic as the HI populations, we want to get the redshifts in with 
    // their associated hydrogen and HeI corrections.  
    
    // Refrain from considering the correction functions if the redshift is past the Recfast++ limit
    if (zval < 50.) return;
    double f_H, f_He, f_tot;
    f_tot = (1.0-Hydrogen_Atoms.X_tot()+parameters.fHe-HeI_Atoms.X_tot())/(cosmos.Xp(zval)+parameters.fHe-cosmos.XHe1s(zval))-1;
    f_H  = (1.0-Hydrogen_Atoms.X_tot())/cosmos.Xp(zval)-1;
    f_He = (parameters.fHe-HeI_Atoms.X_tot())/(parameters.fHe-cosmos.XHe1s(zval))-1;
    // Write correction functions to the file along with the redshift
    ofile << scientific << zval << " " << f_tot << " " << f_H << " " << f_He;
    ofile << endl;
    
    return;
}

//--------------------------------------------------------------------------------------------
void write_array_to_file(ofstream &ofile, double zval, double *X, int neq)
{
    ofile << scientific << zval << " " << cosmos.TCMB(zval);
    // write all the variables
    for(int l=0; l<neq; l++) ofile << scientific << " " << X[l];
    ofile << endl;
    
    return;
}

//--------------------------------------------------------------------------------------------
void write_X_Recfast_z_to_file(ofstream &ofile, double zval, Cosmos &cosmos)
{
    ofile << scientific << zval << " " 
          << cosmos.Xe_Seager(zval) << " " 
          << cosmos.Te(zval) << " "
//        << cosmos.dXe_dz_Seager(zval) << " " << cosmos.X1s(zval) << " " 
//        << cosmos.dX1s_dz(zval) << " " << cosmos.Te_Tg(zval) << " " << cosmos.H(zval)
//        << cosmos.Xp(zval) << " " << parameters.fHe-cosmos.XHe1s(zval)
          << endl;
    
    return;
}

//--------------------------------------------------------------------------------------------
void write_Xe_Te_to_file(ofstream &ofile, double zval, double Xe, double Te)
{
    ofile << scientific << zval << " " << Xe << " " << Te << endl;
    return;
}

//--------------------------------------------------------------------------------------------
// write spectral distortion from EMLA (added July 2nd 2012)
//--------------------------------------------------------------------------------------------
void write_HI_distortion_to_file(ofstream &ofile, double zval)
{
    if(zval>=parameters.zstart) return;
    
    double d1=0.0, d2=0.0;
    Transition_Data_A T=Hydrogen_Atoms.Level(2, 1).Get_Data_A(1, 0);
    ODE_effective::evaluate_Ly_n_channel(zval, cosmos.TCMB(zval), 
                                         Hydrogen_Atoms.X(1, 0), Hydrogen_Atoms.X(2, 1), 
                                         cosmos.NH(zval), cosmos.H(zval), 
                                         T.A21, T.lambda21, T.Dnu, d1, d2);    

    ofile << scientific << zval;
    ofile << scientific << " " << d1;

    d1=0.0; d2=0.0;
    ODE_effective::evaluate_2s_two_photon_decay(zval, cosmos.TCMB(zval), 
                                                Hydrogen_Atoms.X(0), Hydrogen_Atoms.X(1), 
                                                Hydrogen_Atoms.Level(1).Get_Dnu_1s(), 
                                                const_HI_A2s_1s, d1, d2);
    ofile << scientific << " " << d1 << endl;
    
    return;
}

//--------------------------------------------------------------------------------------------
void write_files(double z, ofstream *Xfiles)
{
    if(write_Xe_Te_sol)
    {
        //==========================================================
        // only Xe and Te
        //==========================================================
        if(iteration_count==0) write_X_Recfast_z_to_file(Xfiles[0], z, cosmos);
        write_Xe_Te_to_file(Xfiles[1], z, Level_I.X[0], Level_I.X[Level_I.neq-1]*cosmos.TCMB(z));
    }
    //==============================================================
    // only rho
    //==============================================================
//  write_array_to_file(Xfiles[2], z, &Level_I.X[Level_I.neq-1], 1);
 
    //==============================================================
    // Xe & hydrogen populations
    //==============================================================
    if(write_populations) write_HI_X_array_to_file(Xfiles[3], z, Level_I.X, Level_I.nHIeq, Level_I.neq);

    //==============================================================
    // write spectral distortion (added July 2nd 2012)
    //==============================================================
    if(write_HI_distortion) write_HI_distortion_to_file(Xfiles[5], z);
    
    //==============================================================
    // Helium populations
    //==============================================================
//  write_array_to_file(Xfiles[4], z, &Level_I.X[Level_I.nHIeq+1], Level_I.neq-(Level_I.nHIeq+1)-1);
    //
//  write_HI_HeI_1s_array_to_file(Xfiles[5], z);    
    
    //==============================================================
    // write correction functions between Recfast++ and CosmoRec LH
    //==============================================================
    if (write_corrections) write_correction_functions_to_file(Xfiles[6], z, cosmos);
    
    //==============================================================
    // write Xe for HI and HeI separately
    //==============================================================
    if (write_split) write_HI_HeI_Xe_to_file(Xfiles[7],z,Level_I.X[0]);
    
    return;
}


