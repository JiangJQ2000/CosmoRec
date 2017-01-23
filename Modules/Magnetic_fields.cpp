//========================================================================================
// Author: Jens Chluba (Oct. 2014)
// Johns Hopkins University
//========================================================================================

//========================================================================================
// This simple module allows to add the effect of heating by magnetic fields.
// Expressions from Sethi & Subramanian, 2005 and Kunze & Komatsu, 2013 are used.
//========================================================================================

//========================================================================================
// convert parameters between KK and PF parametrizations
//========================================================================================
void set_PMF_parameters_KK_PF(double B0, double nB, double h0)
{
    if(B0>0.0)
    {
        //======================================================================
        // KK parametrization and interpretation; B0 and nB unchanged; (default)
        //======================================================================
        double kD=286.91/B0;                          // in Mpc^-1; cos(theta)=1
        
        //======================================================================
        // set dissipation scale and convert B_lam --> B0 when using PF vars
        //======================================================================
        if(KK_or_PF_Lorentz)
        {
            double lambda_Mpc=1.0, k_lambda=TWOPI/lambda_Mpc;
            double B_lambda=B0, nBp3=nB+3.0, nBp5=nB+5.0;
            
            kD=pow(2.9e+4*h0, 1.0/nBp5)/pow(B_lambda, 2.0/nBp5)*pow(k_lambda, nBp3/nBp5);
            //kD=pow(5.5e+4*h0, 1.0/nBp5)/pow(B_lambda, 2.0/nBp5)*pow(k_lambda, nBp3/nBp5);
            B0=sqrt( B_lambda*B_lambda*pow(kD*lambda_Mpc, nBp3)/Gamma_JC(nBp5/2.0) );
        }
        
        B0_CosmoRec=RECFAST_atomic_constants::RF_B0_heating=B0;
        nB_CosmoRec=RECFAST_atomic_constants::RF_nB_heating=nB;
        kD_CosmoRec=RECFAST_atomic_constants::RF_kD_heating=kD;
    }
    else
    {
        B0_CosmoRec=RECFAST_atomic_constants::RF_B0_heating=0.0;
        nB_CosmoRec=RECFAST_atomic_constants::RF_nB_heating=0.0;
        kD_CosmoRec=RECFAST_atomic_constants::RF_kD_heating=0.0;
        RECFAST_atomic_constants::RF_eval_ion_Tr = 0; // do old recfast evaluation
        //RECFAST_atomic_constants::RF_eval_ion_Tr = 1; // do new recfast evaluation
    }
    
    return;
}

//========================================================================================
// magnetic field decay rate following Kunze & Komatsu, 2013.
//========================================================================================
double Gamma_decay(double z, double Hz, double B0, double nB, double kd)
{
    double zi=1088.0;
    
    // B0 is in nG
    double m=2.0*(nB+3.0)/(nB+5.0);
    double td_ti=14.8/B0/kd;
    //double rhoB=9.545e-8*B0*B0*cosmos.rho_g_g_cm3(z)*const_cl*const_cl;
    double rhoB=B0*B0*1.0e-18/(2.0*FOURPI)*pow(1.0+z, 4);
    
    double Gamma=0.0;
    double a=log(1.0+td_ti);
    double sig=5.0e+3;        // width if the Gaussian transition
    int trans=1;              // choose type of transition.
    
    if(z<zi) Gamma=rhoB*Hz*1.5*m*pow(a, m)/pow(a+1.5*log((1.0+zi)/(1.0+z)), m+1.0);
    // make smooth transition from no heating to heating case
    else if(z>=zi && trans==0)
    {
        double G_zi=rhoB*Hz*1.5*m/a;
        //Gamma=0.0;
        //Gamma=G_zi*pow((1.0+zi)/(1.0+z), 4);
        //Gamma=G_zi*exp(-pow(z-zi, 2)/sig);
        Gamma=G_zi*exp(-pow(z-zi, 2)/sig)*pow((1.0+zi)/(1.0+z), 4);
    }
    else if(z>=zi && trans==1)
    {
        // smooth 1. derivative transition. Gaussian at z>=z0, polynomial at zi<=z<=z0
        double G_zi=rhoB*Hz*1.5*m/a*pow((1.0+zi)/(1.0+z), 4);
        double dG_dzi=G_zi*(4.0+1.5*(1.0+m)/a)/(1.0+zi);
        double k=7.0;
        double z0=zi*1.001;
        double g=dG_dzi/( k*pow(zi-z0, k-1) );
        double g0=G_zi-g*pow(zi-z0, k);

        if(z<=z0) Gamma=g0+g*pow(z-z0, k);
        else Gamma=g0*exp(-pow(z-z0, 2)/sig);
    }

    return Gamma;
}

//========================================================================================
// magnetic field heating rate for ambipolar diffusion Kunze & Komatsu, 2013.
//========================================================================================
double Gamma_ambi(double z, double Tm, double Xp, double B0, double nB, double kd)
{
    // fit to the integral KK
    double L2=0.8313*pow(nB+3.0, 1.105)*(1.0-1.020e-2*(nB+3.0));
    
    // fit to the integral of Finelli
    if(KK_or_PF_Lorentz)
        L2=0.6615*pow(nB+3.0, 0.8874)*(1.0-0.1367*(nB+3.0)+7.574e-3*pow(nB+3.0, 2));
    
    // similar to Sethi & Subramanian, 2005
    //double L2=1.0;

    // B0 is in nG
    double rhoB=B0*B0*1.0e-18/(2.0*FOURPI)*pow(1.0+z, 4);
    double sig=0.649*pow(Tm, 0.375)*1.0e-9;
    double g=sig/(2.0*const_mH_gr);
    double rhobar=const_mH_gr*cosmos.Nb(z);
    double Gamma=(1.0-Xp)/Xp/g*pow((1.0+z)*kd/const_Mpc*rhoB/rhobar, 2)*L2;
    
    // fit that reproduces KK. It seems this is ~1/(8pi)^2 times smaller than the correct rate
    // conversion J/m^3 --> 10 ergs / cm^3
    //double Gamma=10.0*(1.0-Xp)/Xp*pow(1.0+z, 4)/pow(Tm, 0.375)*5.7e-40*L2*pow(B0/3.0, 2);
    
    //cout << Gamma/L2/((1.0-Xp)/Xp) << " "
    //     << 10.0*pow(1.0+z, 4)/pow(Tm, 0.375)*5.7e-40*pow(B0/3.0, 2)*pow(2.0*FOURPI, 2) << endl;
    
    return Gamma;
}

//========================================================================================
// evalutation of ODE terms
//========================================================================================
void evaluate_magnetic_field_heating(double z, double fHe, double Tg, double Tm,
                                     double XHII, double XHeII,
                                     double B0, double nB, double kD,
                                     double &dX1s_dt,
                                     double &drho_dt)
{
    double Xe=XHII+XHeII, NH=cosmos.NH(z);
//    double fheat=2.0/3.0/const_kB/NH;  // Xe term canceled
    double fheat=2.0/3.0/const_kB/(1.0+fHe+Xe)/NH;
    
    double dE_dt=0.0;
    if(decaying_turb_on)  dE_dt+=Gamma_decay(z, cosmos.H(z), B0, nB, kD);
    if(ambipolar_heat_on) dE_dt+=Gamma_ambi (z, Tm, XHII, B0, nB, kD);
    
    //====================================================================================
    // output rate
    //====================================================================================
/*    double za[1000];
    init_xarr(1200.0, 1000.0, za, 1000, 0, 0);
    
    ofstream ofile("./outputs.Mag_test/test.dat");
    
    for(int k=0; k<1000; k++)
        ofile << za[k] << " " << Gamma_decay(za[k], cosmos.H(za[k]), B0, nB, kD) << endl;
    ofile.close();
    exit(0);
*/    
    //====================================================================================
    // add collisional ionization for H & He (kci in m^3 / sec from Bell et al. 1983)
    // Usually this does not make a difference due to the exponential suppression.
    //====================================================================================
    double kciH =5.466e-9*1.07*sqrt(Tm/1.0e+4)*exp(-157821.462/Tm);
    //double kciHe=5.466e-9*0.37*sqrt(Tm/1.0e+4)*exp(-285471.174/Tm);
    
    dX1s_dt-=kciH*Xe*(1.0-XHII)*NH;
    
    // cooling term
    dE_dt-=13.6*1.602176487e-19*kciH*Xe*(1.0-XHII)*NH*NH;

    drho_dt+=fheat*dE_dt/Tg;  // rho = Te/Tg --> extra factor of 1/Tg

    return;
}

//========================================================================================
//========================================================================================
