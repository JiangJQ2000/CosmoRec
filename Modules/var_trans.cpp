//===========================================================================================================
//
// Level I <--> III
//
//===========================================================================================================

//===========================================================================================================
// variable: ( si=Xi/fi )
//===========================================================================================================

//===========================================================================================================
void transform_y2X(const Data_Level_III &LIII, double z, Data_Level_I &LI)
{
    for(int k=0; k<LI.neq; k++) LI.X[k]=LIII.X[k];
    return;
}

//===========================================================================================================
void transform_g_tX2g_ty(const Data_Level_I &LI, double z, Data_Level_III &LIII)
{ 
  for(int k=0; k<LIII.neq; k++) LIII.g[k]=LI.g[k];
  return;
}


//===========================================================================================================
//
// some auxillary routines for all variables 
// below the Data structures are not global!!!
//
//===========================================================================================================
void copy_y_to_X(const Data_Level_III &LIII, double z, Data_Level_I &LI)
{
    for(int k=0; k<LI.neq; k++) LI.X[k]=LIII.X[k];
    return;
}

//===========================================================================================================
void initialize_y(const Data_Level_I &LI, Gas_of_Atoms &H_Atoms, Gas_of_HeI_Atoms &HeIA, 
                  double z, Data_Level_III &LIII, int wait=0)
{
    for(int k=0; k<LIII.neq; k++)
    {
        if(flag_He==0 && LIII.index_HeI<=k && k<LIII.index_HeI+LIII.nHeIeq) continue;
        
        LIII.X[k]=LI.X[k];
    }
    
    return;
}
