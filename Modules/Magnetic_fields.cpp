//========================================================================================
// Author: Jens Chluba (Oct. 2014)
// Johns Hopkins University
//========================================================================================

//========================================================================================
// decide in which mode to run
//========================================================================================
bool decaying_turb_on =1;
bool ambipolar_heat_on=0;

//========================================================================================
// This simple module allows to add the effect of heating by magnetic fields.
// Expressions from Sethi & Subramanian, 2005 and Kunze & Komatsu, 2013 are used.
//========================================================================================

//========================================================================================
// magnetic field decay rate following Kunze & Komatsu, 2013.
//========================================================================================
double Gamma_decay(double z, double Hz, double B0, double nB)
{
    double zi=1088.0;
    
    // B0 is in nG
    double m=2.0*(nB+3.0)/(nB+5.0);
    double kd=286.91/B0;                                        // in Mpc^-1; cos(theta)=1
    double td_ti=14.8/B0/kd;
    //double rhoB=9.545e-8*B0*B0*cosmos.rho_g_g_cm3(z)*const_cl*const_cl;
    double rhoB=B0*B0*1.0e-18/(2.0*FOURPI)*pow(1.0+z, 4);
    
    // make smooth transition from no heating to heating case
    if(z>=zi){ rhoB*=exp(-pow(z-zi, 2)/5.0e+3); z=zi; }
    
    double Gamma=rhoB*Hz*1.5*m*pow(log(1.0+td_ti), m)
                   /pow(log(1.0+td_ti)+1.5*log((1.0+zi)/(1.0+z)), m+1.0);
 
    return Gamma;
}

//========================================================================================
// magnetic field heating rate for ambipolar diffusion Kunze & Komatsu, 2013.
//========================================================================================
double Gamma_ambi(double z, double Tm, double Xe, double B0, double nB)
{
    // fit to the integral KK
    double L2=0.8313*pow(nB+3.0, 1.105)*(1.0-1.020e-2*(nB+3.0));
    
    //fit to the integral for Finelli
    //double L2=0.6615*pow(nB+3.0, 0.8874)*(1.0-0.1367*(nB+3.0)+7.574e-3*pow(nB+3.0, 2));
    
    // Sethi & Subramanian, 2005
    //double L2=1.0;

    // B0 is in nG
    double kd=286.91/B0;                                  // in Mpc^-1; cos(theta)=1
    double rhoB=B0*B0*1.0e-18/(8.0*RF_PI)*pow(1.0+z, 4);
    double sig=0.649*pow(Tm, 0.375)*1.0e-9;
    double g=sig/(2.0*const_mH_gr);
    double rhobar=const_mH_gr*cosmos.Nb(z);
    double Gamma=(1.0-Xe)/Xe/g*pow((1.0+z)*kd/const_Mpc*rhoB/rhobar, 2)*L2;
    
    // fit that reproduces KK. It seems this is ~1/(8pi)^2 times smaller than the correct rate
    // J/m^3 --> 10 ergs / cm^3
    //double Gamma=10.0*(1.0-Xe)/Xe*pow(1.0+z, 4)/pow(Tm, 0.375)*5.7e-40*L2*pow(B0/3.0, 2);
    
    //cout << Gamma/L2/((1.0-Xe)/Xe) << " "
    //     << 10.0*pow(1.0+z, 4)/pow(Tm, 0.375)*5.7e-40*pow(B0/3.0, 2)*pow(2.0*FOURPI, 2) << endl;
    
    return Gamma;
}

//========================================================================================
// evalutation of ODE terms
//========================================================================================
void evaluate_magnetic_field_heating(double z, double fHe, double Tg, double Tm,
                                     double XHII, double XHeII,
                                     double B0, double nB,
                                     double &dX1s_dt,
                                     double &drho_dt)
{
    double Xe=XHII+XHeII, NH=cosmos.NH(z);
    double fheat=2.0/3.0/const_kB/(1.0+fHe+Xe)/NH;
//    double fheat=2.0/3.0/const_kB/NH;  // Xe term canceled
    
    double dE_dt=0.0;
    if(decaying_turb_on)  dE_dt+=Gamma_decay(z, cosmos.H(z), B0, nB);
    if(ambipolar_heat_on) dE_dt+=Gamma_ambi (z, Tm, Xe, B0, nB);
    
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
