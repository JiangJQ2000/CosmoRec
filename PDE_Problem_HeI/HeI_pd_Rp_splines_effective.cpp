//===========================================================================================================
// Authors Jeffrey Fung & Jens Chluba Feb/March 2010
//
// comment: HeI populations have to be pre-loaded; Also the effective rates have to be initialized
//          Also note that here calc_Xi and all other splines are using the resolved level index 
//          instead of the HeI-level index!
//===========================================================================================================

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>

#include "HeI_pd_Rp_splines_effective.h"
#include "get_effective_rates.HeI.h"
#include "Load.Populations.HI.h"
#include "Load.Populations.HeI.h"
#include "HeI_Atom.h"
#include "Cosmos.h"
#include "physical_consts.h"
#include "routines.h"

using namespace std;

//===========================================================================================================
static bool HeI_pd_Rp_splines_effective_spline_memory_was_allocated=0;
int HeI_pd_Rp_splines_effective_verbose=1;

//===========================================================================================================
void Set_HeI_pd_Rp_splines_effective_verbosity(int v)
{ HeI_pd_Rp_splines_effective_verbose=v; return; }

//===========================================================================================================
//
// compute pd and Rp/Rm
//
//===========================================================================================================
void calc_HeI_Rp_Rm_pd_all(double z, Gas_of_HeI_Atoms &HeIA, double Tg, double Te, double Ne, 
                          double Xp, double NH, vector<double> &Rp_Rm, vector<double> &pd)
{
    int ntot=Rp_Rm.size();
    
    vector<double> Ai(ntot);
    vector<double> Bi(ntot), Bitot(ntot);
    vector<vector<double> > RijVec(ntot, vector<double>(ntot));

    //=======================================================================================================
    get_rates_all_HeI(Tg, Ai, Bi, Bitot, RijVec);

    //=======================================================================================================
    // Ai[nr] :: effective recombination rate to level nr
    // Bi[nr] :: effective photonionization rate from level nr
    // RijVec[i][j] :: effective transition rate to between level i and j
    // Bitot[nr] :: total rate for exciting level nr
    //=======================================================================================================
    
    double Rp, Rm;  
    for(int nr=0; nr<ntot; nr++)
    {
        int index_HeI_i=get_HeI_index(nr);
        
        //===================================================================================================
        // compute Rp & Rm
        //===================================================================================================
        Rp=Ne*Xp*Bi[nr]*Ai[nr];
        Rm=Bitot[nr];

        for(int j=0; j<nr; j++) Rp+=RijVec[j][nr]*calc_HeI_Xi(z, j);
        for(int j=nr+1; j<ntot; j++) Rp+=RijVec[j][nr]*calc_HeI_Xi(z, j);
 
        Rp_Rm[nr]=Rp*NH/Rm;
        
        if(index_HeI_i==1) pd[nr]=Rm/(const_HeI_A2s_1s+Rm);
        else pd[nr]=Rm/(HeIA.Get_A(HeIA.Get_n(index_HeI_i), HeIA.Get_l(index_HeI_i), 
                                   HeIA.Get_S(index_HeI_i), HeIA.Get_J(index_HeI_i), 1, 0, 0, 0)+Rm);
    }

    return;
}

//===========================================================================================================
//
// For precalulated death probabilities and emission rates
//
//===========================================================================================================
int HeI_Rp_nz_eff=2000;

struct HeI_Rp_pd_eff_Data
{
    double *Rp;
    int memindex;
};

vector<HeI_Rp_pd_eff_Data > HeI_pd_eff_Data;
vector<HeI_Rp_pd_eff_Data > HeI_Rp_eff_Data;
vector<HeI_Rp_pd_eff_Data > HeI_Dnem_eff_Data;

//===========================================================================================================
void clear_HeI_pd_Rp_effective_splines()
{
    if(HeI_pd_Rp_splines_effective_spline_memory_was_allocated)
    {
        for(int k=0; k<(int)HeI_pd_eff_Data.size(); k++) 
            free_spline_JC(HeI_pd_eff_Data[k].memindex, "set_up_splines_for_HeI_pd_Rp_effective");
    }
    
    HeI_pd_Rp_splines_effective_spline_memory_was_allocated=0;

    return;
}

//===========================================================================================================
void clear_HeI_pd_Rp_effective_memory()
{
    clear_HeI_pd_Rp_effective_splines();
    
    HeI_pd_eff_Data.clear();
    HeI_Rp_eff_Data.clear();
    HeI_Dnem_eff_Data.clear();
    
    return;
}

//===========================================================================================================
void set_up_splines_for_HeI_pd_Rp_effective(double zend, double zstart,
                                            Cosmos &cosmos, Gas_of_HeI_Atoms &HeA)
{   
    //===========================================================================
    int ntot=get_number_of_resolved_levels_HeI();
    int nz=HeI_Rp_nz_eff;
    double *Rp_zarr=new double[nz];
        
    //===========================================================================
    // check if spline memory was already allocated for as many levels
    //===========================================================================
    if(HeI_pd_Rp_splines_effective_spline_memory_was_allocated)
    {
        if(ntot+1!=(int)HeI_pd_eff_Data.size()) clear_HeI_pd_Rp_effective_splines();
    }
    else 
    {
        //=======================================================================
        // prepare spline data
        //=======================================================================
        HeI_pd_eff_Data.clear();
        HeI_Rp_eff_Data.clear();
        HeI_Dnem_eff_Data.clear();
        HeI_pd_eff_Data.resize(ntot+1);
        HeI_Rp_eff_Data.resize(ntot+1);
        HeI_Dnem_eff_Data.resize(ntot+1);
        
        //=======================================================================
        // allocate memory for pd && Rp-values
        //=======================================================================
        for(int i=0; i<ntot+1; i++)
        { 
            HeI_pd_eff_Data[i].Rp=new double[nz]; 
            HeI_Rp_eff_Data[i].Rp=new double[nz]; 
            HeI_Dnem_eff_Data[i].Rp=new double[nz]; 
        }
    }   
    
    //===========================================================================
    // fill z-array (z0 < z1 < ... < zn)
    //===========================================================================
    init_xarr(zend/1.01, zstart*1.01, Rp_zarr, nz, 0, 0);
    
    if(HeI_pd_Rp_splines_effective_verbose>=1) cout << endl;

    //===========================================================================
    // fill ALL y-arrays
    //===========================================================================
    double z, Tg, Te, NH;
    vector<double> Rp_Rm(ntot);
    vector<double> pd(ntot);
    vector<double> Dnem_arr(ntot);
    Transition_Data_HeI_A T;
    //
    for(int i=0; i<nz; i++)
    {
        z=Rp_zarr[i];
        Tg=cosmos.TCMB(z);
        Te=Tg*calc_HI_rho(z);
        NH=cosmos.NH(z);

        calc_HeI_Rp_Rm_pd_all(z, HeA, Tg, Te, NH*calc_HI_Xe(z), cosmos.fHe()-calc_HeI_X1s(z), NH, Rp_Rm, pd);

        //=====================================================================
        // common data
        //=====================================================================
        T=HeA.Get_Trans_Data(1, 1, 0, 0, 0);
        double x_c=const_h_kb*T.Dnu/Tg;
        double exp_21=exp(-x_c);
        double N1s=NH*calc_HeI_X1s(z);
        double N2s=NH*calc_HeI_Xi(z, 0);
        double N2p=NH*calc_HeI_Xi(z, 1);

        //=====================================================================
        // compute the 2s Rp_Rm using the net rate equation
        //=====================================================================
        int res_i=get_res_index_HeI(1);

        //=====================================================================
        // case pd_2s==1
        //=====================================================================
        Rp_Rm[res_i]=N2s;                // this is the normal multi-level term
        Dnem_arr[res_i]=N2s/N1s-exp_21;
        
        //=====================================================================
        // compute the 2p Rp_Rm using the net rate equation
        //=====================================================================
        res_i=get_res_index_HeI(2);
        T=HeA.Get_Trans_Data(2, 1, 0, 0, 0);
        x_c=const_h_kb*T.Dnu/Tg;
        exp_21=exp(-x_c);
        //
        double npl_c=exp_21;
        double nL=N2p/3.0/N1s;      
        double tauS=T.A21*pow(T.lambda21, 3)/(8.0*PI*cosmos.H(z))*(N1s*3.0 - N2p);      
        double PS=(1.0-exp(-tauS))/tauS;
        //
        double Dnem=(nL-npl_c)*(1.0+(1.0/pd[res_i]-1.0)*PS);
        Rp_Rm[res_i]=(Dnem+npl_c)*3.0*N1s;
        Dnem_arr[res_i]=Dnem;
        
        //=====================================================================
        // compute the nP - 1s and nD - 1s Rp_Rm using the net rate equation
        //=====================================================================
        for(int res_i=2; res_i<ntot; res_i++)
        {
            int HeI_i=get_HeI_index(res_i);
            double Ni=NH*calc_HeI_Xi(z, res_i);
            T=HeA.Get_Trans_Data(HeI_i, 1, 0, 0, 0);
            double wr=1.0*T.gwp/HeA.Get_gw(HeI_i);
            //
            x_c=const_h_kb*T.Dnu/Tg;
            exp_21=exp(-x_c);
            npl_c=exp_21;
            //
            nL=wr*Ni/N1s;
            tauS=T.A21*pow(T.lambda21, 3)/(8.0*PI*cosmos.H(z))*(N1s/wr - Ni);
            PS=(1.0-exp(-tauS))/tauS;
            //
            Dnem=(nL-npl_c)*(1.0+(1.0/pd[res_i]-1.0)*PS);
            Rp_Rm[res_i]=(Dnem+npl_c)*N1s/wr;
            Dnem_arr[res_i]=Dnem;
        }
        
        //=====================================================================
        // copy all stuff to arrays for splines
        //=====================================================================
        for(int nr=0; nr<ntot; nr++)
        { 
            HeI_pd_eff_Data[nr].Rp[i]=log(pd[nr]); 
            HeI_Rp_eff_Data[nr].Rp[i]=log(Rp_Rm[nr]); 
            HeI_Dnem_eff_Data[nr].Rp[i]=Dnem_arr[nr]; 
        }
        
        //=====================================================================
        // set the the entry with res-index=ntot to zero 
        // (needed for spline interpolation)
        //=====================================================================
        HeI_pd_eff_Data[ntot].Rp[i]=HeI_Rp_eff_Data[ntot].Rp[i]=HeI_Dnem_eff_Data[ntot].Rp[i]=0.0;
        
        if(HeI_pd_Rp_splines_effective_verbose>=1) 
            cout << "\033[2K" << " z= " << z << endl << "\033[1A";
    }
    
    if(HeI_pd_Rp_splines_effective_verbose>=1) 
        cout << " set_up_splines_for_HeI_pd_Rp:: rates tabulated " << endl;
    
    //========================================================================
    // calc splines + allocalate memory (if necessary)
    //========================================================================
    if(HeI_pd_Rp_splines_effective_spline_memory_was_allocated)
    {
        for(int n=0; n<ntot+1; n++)
        {
            update_spline_coeffies_JC(HeI_pd_eff_Data[n].memindex, nz, Rp_zarr, HeI_pd_eff_Data[n].Rp);
            update_spline_coeffies_JC(HeI_Rp_eff_Data[n].memindex, nz, Rp_zarr, HeI_Rp_eff_Data[n].Rp);
            update_spline_coeffies_JC(HeI_Dnem_eff_Data[n].memindex, nz, Rp_zarr, HeI_Dnem_eff_Data[n].Rp);
        }
    }
    else 
    {
        for(int n=0; n<ntot+1; n++)
        {
            HeI_pd_eff_Data[n].memindex=calc_spline_coeffies_JC(nz, Rp_zarr, HeI_pd_eff_Data[n].Rp, 
                                                               "set_up_splines_for_HeI_pd_Rp :: pd "
                                                               + int_to_string(n));

            HeI_Rp_eff_Data[n].memindex=calc_spline_coeffies_JC(nz, Rp_zarr, HeI_Rp_eff_Data[n].Rp, 
                                                               "set_up_splines_for_HeI_pd_Rp :: Rp "
                                                               + int_to_string(n));

            HeI_Dnem_eff_Data[n].memindex=calc_spline_coeffies_JC(nz, Rp_zarr, HeI_Dnem_eff_Data[n].Rp, 
                                                                 "set_up_splines_for_HeI_pd_Rp :: Dnem "
                                                                 + int_to_string(n));
        }
        
        HeI_pd_Rp_splines_effective_spline_memory_was_allocated=1;
    }
    
    if(HeI_pd_Rp_splines_effective_verbose>=1) 
        cout << " set_up_splines_for_HeI_pd_Rp:: splines done." << endl; 
    
    //========================================================================
    // clean up
    //========================================================================
    delete[] Rp_zarr;
    
    return;
}

//===========================================================================================================
//
// pd-functions
//
//===========================================================================================================
double calc_HeI_pd_i_splines_effective(double z, int ires)
{ return exp(calc_spline_JC(z, HeI_pd_eff_Data[ires].memindex)); }

//===========================================================================================================
//
// Rp-functions
//
//===========================================================================================================
double calc_HeI_Rp_Rm_i_splines_effective(double z, int ires)
{ return exp(calc_spline_JC(z, HeI_Rp_eff_Data[ires].memindex)); }

//===========================================================================================================
//
// Dnem-functions
//
//===========================================================================================================
double calc_HeI_Dnem_i_splines_effective(double z, int ires)
{ return calc_spline_JC(z, HeI_Dnem_eff_Data[ires].memindex); }

