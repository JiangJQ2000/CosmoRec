//==================================================================================================
// setting for two-gamma and Raman corrections
//==================================================================================================
int ntg_max=0;
int nR_max=0;
bool hydrogen_and_helium_atom_are_set=0;

//==================================================================================================
// for data output
//==================================================================================================
string path, addname;
string detailedname;
ofstream *Xfiles=new ofstream[9];
vector<vector<double> > output_CosmoRec;
vector<vector<double> > pass_on_the_Solution_CosmoRec;
vector<vector<double> > pass_on_the_Solution_CosmoRec_HeI; //JF 

//==================================================================================================
// output files
//==================================================================================================
bool write_populations=0;   // output solution for the populations of the resolved states
bool write_corrections=0;   // output correction functions from Recfast++ and CosmoRec (LH)
bool write_Xe_Te_sol=1;     // output Xe & Te
bool write_HI_distortion=0; // output approximate HI distortion (added July 2nd 2012)
bool write_split=0;         // output that splits the Xe for H and He so they can be individually analysed
bool write_RF_visi=0;       // output thomson visibility function in RF mode (LH)

//==================================================================================================
// path to the effective rate database
//==================================================================================================
string Rec_database_path=COSMORECDIR+COSMOREC_BASE_PATH;

//==================================================================================================
// below z=200 extend the solution with Recfast system
//==================================================================================================
int flag_compute_Recfast_part=0;
double ze_Recfast=0.0;
int nz_Recfast_part=100;
double Xi_HeI_switch=1.0e-7;

//==================================================================================================
// Variables for module that accounts for the diffusion corrections to escape probabilities
//==================================================================================================
bool Diffusion_correction=0;
bool Diffusion_correction_is_on=0;
int Diffusion_flag;
//
int Diff_iteration_max=2, Diff_iteration_min=0;
int iteration_count=Diff_iteration_min;
// 
double Diff_corr_zmin=500.0;
double Diff_corr_zmax=2000.0;

//==================================================================================================
// Helium diffusion correction (added Feb-May 2011 by Jens Chluba & Jeffrey Fung)
//==================================================================================================
bool Diffusion_correction_HeI=0;
bool Diffusion_correction_HeI_is_on=0;
double HeI_was_switched_off_at_z=0.0;

//==================================================================================================
// Variables for dark matter annihilation module
//==================================================================================================
bool DM_annihilation=0;
double fDM_CosmoRec=0.0;

//==================================================================================================
// Variables for magnetic field heating module
//==================================================================================================
bool magnetic_fields=0;
double B0_CosmoRec=0.0, nB_CosmoRec=0.0, kD_CosmoRec=0.0;

bool decaying_turb_on =0;  // default mode set in "./Modules/Magnetic_fields.cpp"
bool ambipolar_heat_on=0;  // default mode set in "./Modules/Magnetic_fields.cpp"
bool KK_or_PF_Lorentz =0;  // choose parametrization of KK (==0) or PF (==1)

// set functions defined in './Modules/Magnetic_fields.cpp'
void Switch_PMF_off();
void set_PMF_parameters_KK_PF(double B0, double nB, int PMF_mode, double h0);

//==================================================================================================
// CosmoRec specific variables for variation of fundamental constant
//==================================================================================================
// Scaling variables for alpha, reduced mass and omission modes (Luke Hart)
struct Fund_Const_scaling
{
    size_t runmode;
    double a_scale;
    double m_scale;
    double p_scale;
    // Scaling coefficients that we multiply rates/cross sections
    double f_a, f_b, f_t;
    
    Fund_Const_scaling() {
        runmode=0;
        a_scale=m_scale=1.0;
        p_scale=0.0;
        f_a=f_b=f_t=1.0;
    }
};

Fund_Const_scaling FC_CosmoRec;

// Switch to determine whether we are rescaling Helium recombination or not
int FC_He_switch = 1;

// set functions defined in './Modules/Variation_constants.cpp'
void Switch_VFC_off();
void Set_VFC_pars(double aS, double mS, double pS, int RM);


//==================================================================================================
// memory for CosmoRec-ODE solver
//==================================================================================================
ODE_solver_Rec::ODE_solver_Solution Sz;
ODE_solver_Rec::ODE_solver_accuracies tols;
ODE_solver_Rec::ODE_Solver_data ODE_Solver_info_CR;

//==================================================================================================
// for passing on the data about cosmology
//==================================================================================================
struct Parameters_form_initfile
{
    int verbosity;
    int nz;
    double zstart, zend;

    string path, addname;

    Cosmos::parameters cosm_pars;
    
    Parameters_form_initfile()
    {
        verbosity=show_CosmoRec_mess;
        nz=2000;
        zstart=3000; zend=0.0;
    }
    
    void show()
    {
        cout << " Parameters_form_initfile :: show()" << endl;
        cout << " Verbosity level " << verbosity << endl;
        cout << " nz = " << nz << " zstart= " << zstart << " zend= " << zend << endl;
        cosm_pars.show();
    }
};

Parameters_form_initfile parameters;

//==================================================================================================
// Level I: Xi=Ni/NHtot 
//==================================================================================================
struct Data_Level_I
{
    int neq;
    int nHIeq;
    int nHeIeq;
    int index_HI;
    int index_HeI;
    double *X, *g;
};

Data_Level_I Level_I, Level_I_temp, Level_I_temp_JAC;

//==================================================================================================
double *zarr;               // redshift array

//==================================================================================================
// for Cosmology object
//==================================================================================================
Cosmos cosmos;

//==================================================================================================
//==================================================================================================
