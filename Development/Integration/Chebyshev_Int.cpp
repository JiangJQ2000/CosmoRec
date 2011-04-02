//-----------------------------
// Author Jens Chluba Jan 2010
//-----------------------------
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <complex>
#include <vector>

#include "Definitions.h"
#include "Chebyshev_Int.h"
#include "routines.h"

using namespace std;

//====================================================================================
// functions for mapping of f(x) to interval mu=(-1, 1)
//====================================================================================
static double Chebyshev_lambda;
static double Chebyshev_a;
static double Chebyshev_b;
static double Chebyshev_c;

double Chebyshev_x(double mu){ return Chebyshev_lambda*(Chebyshev_b+mu)/(Chebyshev_a-mu); }

double (*Chebyshev_Int_fptr)(double);
double Chebyshev_Int_f1(double mu)
{ 
    double dx_dmu=Chebyshev_c/pow(Chebyshev_a-mu, 2);
    return dx_dmu*Chebyshev_Int_fptr(Chebyshev_x(mu)); 
}

double Chebyshev_test_fcn(double x)
{ 
    return 1.0/(1.0+x*x);
    //return 1.0/(exp(x)-1.0)/x; 
    //return x*x*x*x/(exp(x)-1.0); 
    //return x*x*x/(exp(x)-1.0); 
    //return x*x/(exp(x)-1.0); 
    //return x*x;
}

double Chebyshev_test_fcn_log(double lgx)
{ return exp(lgx)*Chebyshev_test_fcn(exp(lgx)); }

//====================================================================================
// computations of weights & abscissae
//====================================================================================
static bool Chebyshev_table_is_set=false;
vector<vector<double> > Chebyshev_w;
vector<vector<double> > Chebyshev_xi;

double sum_Chebyshev_weights(int n, double theta_k)
{
    double r=0.0;
    for(int m=(n+1)/2; m>=1; m--) r+=sin((2*m-1)*theta_k)/(2*m-1);
    return r*4*sin(theta_k)/(n+1);
}

void compute_nested_Chebyshev_table(int pmax)
{
    cout << " Computing nested Table for Chebychev integration rules up to p= " << pmax << endl; 
    //
    int mess=0;
    Chebyshev_w.clear();
    Chebyshev_xi.clear();
    
    int nmax;
    vector<double> dum;
    double theta_k, tt;
    
    // empty element for p=0
    Chebyshev_w.push_back(dum);
    Chebyshev_xi.push_back(dum);
    
    cout.precision(16);
    for(int p=1; p<=pmax; p++)
    {
        nmax=(int)pow(2.0, p)-1;

        dum.clear();
        if(mess>=1) cout << "\n nmax " << nmax << endl;
        if(mess>=1) cout << " wi " << nmax << endl;
        for(int k=1; k<nmax/2+2; k++) 
        {
            theta_k=PI*k/(nmax+1);
            tt=sum_Chebyshev_weights(nmax, theta_k);
            dum.push_back(tt);
        }
        
        if(mess>=1) for(int i=0; i<(int)dum.size(); i++) cout << scientific << dum[i] << ", ";      
        Chebyshev_w.push_back(dum);

        dum.clear();
        if(mess>=1) cout << "\n xi " << endl;
        for(int k=1, l=1; k<nmax/2+1; k++) 
        {
            theta_k=PI*k/(nmax+1);
            tt=cos(theta_k);

            if(p>1 && l==k)
            { 
                dum.push_back(tt);  
                l+=2; 
            }       
        }
        if(p==1) dum.push_back(0.0);

        if(mess>=1) for(int i=0; i<(int)dum.size(); i++) cout << scientific << dum[i] << ", ";      
        Chebyshev_xi.push_back(dum);
        
        if(mess>=1) cout << endl;
        if(mess>=1) wait_f_r();
    }
    
    if(mess>=1) cout << endl;   
    if(mess>=1) wait_f_r();
    
    Chebyshev_table_is_set=true;
    
    return;
}
 
//====================================================================================
// Functions to compute the integral
//====================================================================================
vector<double> Chebyshev_fvals;

int compute_integral_function_Chebyshev(int p, double (*f)(double), int *neval, double *r)
{
    int debug=0;
    *neval=(int)(2.0*pow(2.0, p-1)-1.0);
    *r=0.0;
    
    //======================================================
    // only add new points
    //======================================================
    double dum, mu;
    int index_f=0, nxi=Chebyshev_xi[p].size();
    vector<double> loc_fvals;
    // copy first point
    loc_fvals.push_back(Chebyshev_fvals[index_f]);  // compute f(mu=0) outside!
    //
    index_f++;
    for(int k=nxi-1; k>0; k--)
    {
        mu=Chebyshev_xi[p][k];
        //---------------------------------
        // calc new point
        dum=f(mu)+f(-mu);
        loc_fvals.push_back(dum);
        //---------------------------------
        // copy old point
        loc_fvals.push_back(Chebyshev_fvals[index_f]);
        index_f++;
    }
    //======================================================
    // finish last point
    //======================================================
    mu=Chebyshev_xi[p][0];
    dum=f(mu)+f(-mu);
    //
    loc_fvals.push_back(dum);
    
    //======================================================
    // save things for next call
    //======================================================
    Chebyshev_fvals=loc_fvals;
    
    int nfvals=(int)Chebyshev_fvals.size();
    for(int k=nfvals-1; k>=0; k--) *r+=Chebyshev_fvals[nfvals-1-k]*Chebyshev_w[p][k];
    
    if(debug>=1) cout << "\n order " << p << " neval= " << *neval << " " << *r << endl; 
    if(debug>=2) wait_f_r();
    
    return 0;
}

//====================================================================================
// Epsilon algorithm
//====================================================================================
void update_remaining_epsilon_table(int m, int n, vector<vector<double> > &epstable)
{
    if(m>0)
    {
        double Del=epstable[m+1][n-1]-epstable[m][n-1];
        cout << " Check " << m << " " << n << " " << Del << " " << epstable[m][n-2]+1.0/Del << endl;
    
        epstable[m].push_back(epstable[m+1][n-2]+1.0/Del);
    
        update_remaining_epsilon_table(m-1, n+1, epstable);
    }
        
    return; 
}

void show_epsilon_table(vector<vector<double> > &epstable)
{
    int mmax=epstable.size();
    
    cout.precision(10);
    for(int m=0; m<mmax; m++)
    {
        cout << m << " ";
        for(int n=0; n<(int)epstable[m].size(); n++) cout << epstable[m][n] << " ";
        cout << endl;
    }

    wait_f_r();
    
    return; 
}

//====================================================================================
// Driver that loops over the different rules
//====================================================================================
int Integrate_using_Chebyshev(double a, double x0, double b, double epsrel, double epsabs, double (*fptr)(double), int *neval, double *r, string mess)
{
    if(Chebyshev_table_is_set==false) compute_nested_Chebyshev_table(12);
    
    int debug=0;
    *r=0.0; 
    if(a>=b) return 2;
    
    if(debug>=2) cout << " Entering chebyshev integrator xl= " << a << " xu= " << b << " x0= " << x0 << endl;
    if(x0<=a || x0>=(a+b)/2.0*(1.0-1.0e-2) ) x0=(3.0*a+b)/4.0; // in this case Lambda ~ 1
    //
    Chebyshev_a=(b-a)/(a+b-2.0*x0);
    Chebyshev_b=x0*(b-a)/((a+b)*x0-2.0*a*b);
    Chebyshev_lambda=((a+b)*x0-2.0*a*b)/(a+b-2.0*x0); // == a/b
    //
    Chebyshev_c=Chebyshev_lambda*(Chebyshev_a+Chebyshev_b);
    //
    if(debug>=3) cout << " a = " << Chebyshev_a << " b= " << Chebyshev_b << " c= " << Chebyshev_c << " lam= " << Chebyshev_lambda << " x0= " << x0 << endl;
    //                            
    Chebyshev_Int_fptr=fptr;
    if(Chebyshev_lambda==0.0) cout << " Integrate_using_Chebyshev:: check x0 " << x0 << " a= " << a << " b= " << b << endl; 
    
    int ifail;
    double fi=Chebyshev_Int_f1(0.0);
    Chebyshev_fvals.clear();
    Chebyshev_fvals.push_back(fi);
    //--------------------------------
    // Integral for 1-point formula
    //--------------------------------
    *r=2*fi;
    *neval=1;
    if(debug>=2) cout << " f(mu=0.0)= " << fi << endl << endl;
    
    double r1, Dr;
    //
    for(int it=2; it<=12; it++)
    {
        ifail=compute_integral_function_Chebyshev(it, Chebyshev_Int_f1, neval, &r1);
        Dr= *r-r1;
        *r=r1;
        
        //======================================================
        // check error
        //======================================================
        if(fabs(Dr)<=max(epsabs, fabs( *r)*epsrel))
        {
            if(debug>=1) cout << " result= " << *r << " neval= " << *neval << " check output-file " << endl;
            if(debug>=2) wait_f_r();
            return 0;
        }
    }
    
    cout << " Integral did not converge using Chebyshev rules " << a << " " << x0 << " " << b << " " << *r << " " << *neval << " mess= " << mess << endl; 
    
    return 1;           
}

