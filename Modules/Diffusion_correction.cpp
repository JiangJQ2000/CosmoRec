//==========================================================================================================
//
// Module to account for the diffusion corrections to escape probabilities
//
//==========================================================================================================
int DPesc_splines_are_setup=0;

struct DPesc_splines_Data
{
    int memindex;
    double zmin, zmax;
    double DPmin, DPmax;
};

DPesc_splines_Data DPesc_splines_Data_Lyn[100];

int index_2g_spline;
int index_R_spline;
int index_Ly_n_spline;
int index_all_spline;

//==========================================================================================================
void setup_DF_spline_data(int n, vector<double> &DPesc_vec_z, vector<double> &DPesc_vec)
{
    int nz=DPesc_vec_z.size();
    //
    double *za=new double[nz];
    double *ya=new double[nz];
    
    for(int i=0; i<nz; i++)
    { 
        za[i]=DPesc_vec_z[nz-1-i]; 
        ya[i]=DPesc_vec[nz-1-i]; 
    }
    
    if(DPesc_splines_are_setup==0) 
        DPesc_splines_Data_Lyn[n].memindex=calc_spline_coeffies_JC(nz, za, ya,
                                                                   "setup_DF_spline_data :: "
                                                                   +int_to_string(n));   
    //
    else update_spline_coeffies_JC(DPesc_splines_Data_Lyn[n].memindex, nz, za, ya);                

    DPesc_splines_Data_Lyn[n].zmin=max(400.0, DPesc_vec_z[nz-1]);
    DPesc_splines_Data_Lyn[n].zmax=min(2000.0, DPesc_vec_z[0]);
    DPesc_splines_Data_Lyn[n].DPmin=DPesc_splines_Data_Lyn[n].DPmax=0.0;
    
    delete [] za;
    delete [] ya;
    
    return;
}

//==========================================================================================================
void setup_DF_spline_data(vector<double> &DF_vec_z, 
                          vector<vector<double> > &DF_Ly_n_vec, 
                          vector<vector<double> > &DF_2g_vec, 
                          vector<vector<double> > &DF_R_vec)
{   
    index_2g_spline=2;
    //
    for(int m=0; m<(int)DF_2g_vec.size(); m++) setup_DF_spline_data(m+index_2g_spline, DF_vec_z, DF_2g_vec[m]);
    index_R_spline=index_2g_spline+DF_2g_vec.size();
    //
    for(int m=0; m<(int)DF_R_vec.size(); m++) setup_DF_spline_data(m+index_R_spline, DF_vec_z, DF_R_vec[m]);
    index_Ly_n_spline=index_R_spline+DF_R_vec.size();
    //
    for(int m=0; m<(int)DF_Ly_n_vec.size(); m++) setup_DF_spline_data(m+index_Ly_n_spline, DF_vec_z, DF_Ly_n_vec[m]);
    index_all_spline=index_Ly_n_spline+DF_Ly_n_vec.size();
    //
    DPesc_splines_are_setup=1;

    return;
}

//==========================================================================================================
double DF_2g_spline(int i, double z)
{
    if(z<=DPesc_splines_Data_Lyn[index_2g_spline+i].zmin) return DPesc_splines_Data_Lyn[index_2g_spline+i].DPmin;
    if(z>=DPesc_splines_Data_Lyn[index_2g_spline+i].zmax) return DPesc_splines_Data_Lyn[index_2g_spline+i].DPmax;
    return calc_spline_JC(z, DPesc_splines_Data_Lyn[index_2g_spline+i].memindex);
}

//==========================================================================================================
double DF_R_spline(int i, double z)
{
    if(z<=DPesc_splines_Data_Lyn[index_R_spline+i].zmin) return DPesc_splines_Data_Lyn[index_R_spline+i].DPmin;
    if(z>=DPesc_splines_Data_Lyn[index_R_spline+i].zmax) return DPesc_splines_Data_Lyn[index_R_spline+i].DPmax;
    return calc_spline_JC(z, DPesc_splines_Data_Lyn[index_R_spline+i].memindex);
}

//==========================================================================================================
double DF_Ly_n_spline(int i, double z)
{
    if(z<=DPesc_splines_Data_Lyn[index_Ly_n_spline+i].zmin) return DPesc_splines_Data_Lyn[index_Ly_n_spline+i].DPmin;
    if(z>=DPesc_splines_Data_Lyn[index_Ly_n_spline+i].zmax) return DPesc_splines_Data_Lyn[index_Ly_n_spline+i].DPmax;
    return calc_spline_JC(z, DPesc_splines_Data_Lyn[index_Ly_n_spline+i].memindex);
}

//==========================================================================================================
//
// Ly-n 1+1 photon correction
//
//==========================================================================================================
double compute_DR_Ly_n_func(int i, double z, double X1s, double Tg, Gas_of_Atoms &HIA, int ni)
{
    double DRij=DF_Ly_n_spline(i, z)*( HIA.X(ni, 1)/X1s/3.0-exp(-const_h_kb*HIA.Level(ni, 1).Get_Dnu_1s()/Tg) );
//    cout << z << " (ni, li) " << ni << " " << 1 << " Ly-n " << DRij << endl;
    return DRij;
}

void evaluate_HI_Diffusion_correction_Ly_n(double z, Gas_of_Atoms &HIA, double NH, double H_z, double Tg, double *g)
{
    if(index_Ly_n_spline==index_all_spline) return;
    
    double DRij;
    double X1s=HIA.X(1, 0);
    int index=0;
    int nmax=8;
    
    for(int ni=2; ni<=nmax && index+index_Ly_n_spline<index_all_spline; ni++)
    {
        DRij=compute_DR_Ly_n_func(index, z, X1s, Tg, HIA, ni);
        g[0]+=-DRij;                          // remove electron from lower level
        g[HIA.Get_Level_index(ni, 1)]+= DRij; // add it to -state
        index++;
    }
    
    return;
}

//==========================================================================================================
//
// two-photon correction
//
//==========================================================================================================
double compute_DR_2g_func(int i, double z, double X1s, double Tg, Gas_of_Atoms &HIA, int ni, int li)
{
    double DRij=DF_2g_spline(i, z)*( HIA.X(ni, li)/X1s/(2.0*li+1.0)-exp(-const_h_kb*HIA.Level(ni, li).Get_Dnu_1s()/Tg) );
//  cout << z << " (ni, li) " << ni << " " << li << " 2g " << DRij << endl;
    return DRij;
}

void evaluate_HI_Diffusion_correction_2_gamma(double z, Gas_of_Atoms &HIA, double NH, double H_z, double Tg, double *g)
{
    if(index_2g_spline==index_R_spline) return;

    double DRij;
    double X1s=HIA.X(1, 0);
    int index=0;
    int nmax=8;
    int lmax=2;

    for(int ni=3; ni<=nmax && index+index_2g_spline<index_R_spline; ni++)
        for(int li=0; li<=lmax; li+=2)
        {
            DRij=compute_DR_2g_func(index, z, X1s, Tg, HIA, ni, li);
            g[0]+=-DRij;                           // remove electron from lower level
            g[HIA.Get_Level_index(ni, li)]+= DRij; // add it to -state
            index++;
        }
    
    return;
}

//==========================================================================================================
//
// Raman correction
//
//==========================================================================================================
double compute_DR_R_func(int i, double z, double X1s, double Tg, Gas_of_Atoms &HIA, int ni, int li)
{
    double DRij=DF_R_spline(i, z)*( HIA.X(ni, li)/X1s/(2.0*li+1.0)-exp(-const_h_kb*HIA.Level(ni, li).Get_Dnu_1s()/Tg) );
//  cout << z << " (ni, li) " << ni << " " << li << " Raman " << DRij << endl;  
    return DRij;
}

void evaluate_HI_Diffusion_correction_R(double z, Gas_of_Atoms &HIA, double NH, double H_z, double Tg, double *g)
{
    if(index_R_spline==index_Ly_n_spline) return;
    
    double DRij;
    double X1s=HIA.X(1, 0);
    int index=0;
    int nmax=7;
    int lmax=2;
    
    for(int ni=2; ni<=nmax && index+index_R_spline<index_Ly_n_spline; ni++)
        for(int li=0; li<=lmax && li<ni; li+=2)
        {
            DRij=compute_DR_R_func(index, z, X1s, Tg, HIA, ni, li);
            g[0]+=-DRij;                           // remove electron from lower level
            g[HIA.Get_Level_index(ni, li)]+= DRij; // add it to -state
            index++;
        }
    
    return;
}

//==========================================================================================================
void evaluate_HI_Diffusion_correction(double z, Gas_of_Atoms &HIA, double NH, double H_z, double Tg, double *g)
{
    evaluate_HI_Diffusion_correction_2_gamma(z, HIA, NH, H_z, Tg, g);
    evaluate_HI_Diffusion_correction_R      (z, HIA, NH, H_z, Tg, g);
    evaluate_HI_Diffusion_correction_Ly_n   (z, HIA, NH, H_z, Tg, g);
    
    return;
}

//==========================================================================================================
//
// Module to account for the 2s-1s corrections
//
//==========================================================================================================
int induced_flag=0;
bool DI1_2s_correction_on=0;
bool DI1_2s_splines_are_setup=0;

struct DI_2s_splines_Data
{
    int memindex;
    double zmin, zmax;
    double DI_2smin, DI_2smax;
};

DI_2s_splines_Data DI1_2s_splines_Data;


//==========================================================================================================
void setup_DI1_2s_spline_data(vector<double> &DI1_2s_vec_z, vector<double> &DI1_2s_vec)
{
    int nz=DI1_2s_vec_z.size();
    //
    DI1_2s_splines_Data.zmin=DI1_2s_vec_z[nz-1];
    DI1_2s_splines_Data.zmax=DI1_2s_vec_z[0];
    DI1_2s_splines_Data.DI_2smin=DI1_2s_vec[nz-1];
    DI1_2s_splines_Data.DI_2smax=DI1_2s_vec[0];
    //
    double *za=new double[nz];
    double *ya=new double[nz];
    
    for(int i=0; i<nz; i++)
    { 
        za[i]=DI1_2s_vec_z[nz-1-i]; 
        ya[i]=DI1_2s_vec[nz-1-i]; 
    }
    
    if(DI1_2s_splines_are_setup==0)
    {
        DI1_2s_splines_Data.memindex=calc_spline_coeffies_JC(nz, za, ya, "DI1_2s");     
        DI1_2s_splines_are_setup=1;
    }
    else update_spline_coeffies_JC(DI1_2s_splines_Data.memindex, nz, za, ya);             
    
    delete [] za;
    delete [] ya;
    
    return;
}


//==========================================================================================================
double DI1_2s_spline(double z)
{
    if(z<=DI1_2s_splines_Data.zmin) return DI1_2s_splines_Data.DI_2smin;
    if(z>=DI1_2s_splines_Data.zmax) return DI1_2s_splines_Data.DI_2smax;
    return calc_spline_JC(z, DI1_2s_splines_Data.memindex);
}

//==========================================================================================================
// correction to 2s-1s two-photon channel
//==========================================================================================================
void evaluate_HI_DI1_DI2_correction(double z, double X1s, double X2s, double nu21, double Tg, double *g)
{
    double x21=const_h_kb*nu21/Tg;  
    double DI1_2s=DI1_2s_spline(z);   
    DI1_2s*=X2s-X1s*exp(-x21);  
    double DRij=const_HI_A2s_1s*DI1_2s;
    
    g[0]+=-DRij; // remove electron from 1s level
    g[1]+= DRij; // add it to 2s-state
    
//  cout << z << " 2s-1s correction " << DRij << " DI1= " << DI1_2s << " " << exp(-x21) << endl;
    
    return;
}

