//===========================================================================================================
//
// Add PCA functions from Farhang et al., 2012 & 2013
// Author: Jens Chluba (Aug 2017)
//
//===========================================================================================================
bool PCA_functions_are_loaded=0;

string PCA_path=COSMORECDIR+"./Modules/PCA_data/Planck_2015_smoothed/";

//===========================================================================================================
vector<double> PCA_zarr;
vector<vector<double> > PCA_EM(5);

//===========================================================================================================
// load PCA functions
//===========================================================================================================
void load_PCA_function(string name, vector<double> &zarr, vector<double> &EM)
{
    double dum;
    ifstream ifile(name.c_str());
    zarr.clear(); EM.clear();
    
    while(ifile)
    {
        ifile >> dum;
        if(!ifile) break;
        zarr.push_back(dum);
        ifile >> dum;
        ifile >> dum; EM.push_back(dum);
    }
    
    return;
}

void load_PCA_functions()
{
    for(int k=1; k<5; k++) load_PCA_function(PCA_path+"EM_"+int_to_string(k)+".dat", PCA_zarr, PCA_EM[k]);
    PCA_functions_are_loaded=1;
    return;
}

void Clear_PCA()
{
    PCA_functions_are_loaded=0;
    PCA_zarr.clear();
    PCA_EM.clear();
    return;
}

//===========================================================================================================
// add PCA modes
//===========================================================================================================
void add_PCA_functions(double mu1, double mu2, double mu3, double mu4,
                       const int nz, double *z_arr, double *Xe_arr)
{
    if(mu1==0.0 && mu2==0.0 && mu3==0.0 && mu4==0.0) return;
    
    if(PCA_functions_are_loaded==0) load_PCA_functions();
    
    int istart=0, npol=2;
    double y[5], dy[5];
    
    for(int k=0; k<nz; k++)
    {
        if(z_arr[k]>PCA_zarr[0] && z_arr[k]<PCA_zarr.back())
        {
            for(int m=1; m<5; m++)
                polint_JC(&PCA_zarr[0], &PCA_EM[m][0], PCA_EM[m].size(),
                          z_arr[k], istart, npol, &y[m], &dy[m]);
            
            Xe_arr[k]+=Xe_arr[k]*(y[1]*mu1+y[2]*mu2+y[3]*mu3+y[4]*mu4);
        }
    }
    
    return;
}

//===========================================================================================================
//===========================================================================================================
