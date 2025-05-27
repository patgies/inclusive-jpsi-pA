#include "def.hpp"
#include "amplitudelib.hpp"
#include "fourier.h"
#include <cmath>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_bessel.h>
#include <fstream>
#include <iostream>

using namespace std;


double SQR(double x) 
{ 
    return x*x;  
}


// Gamma functions collinear limit in appendix B.2.2 of arXiv:1309.7337v2

double Gamma(int kappa, double m, double p, double kx, double ky)
{ 
    double ls = SQR(kx- p/2.) + SQR(ky);  //considering py=0!
    double x = ls  + m*m;                 //energy denominator
    
    if (kappa == 1)
    {
      return (2./(3.*m*m*m))*SQR(1. - (m*m/x)); 
    }
    if (kappa == 2)
    {
      return (ls/(3.*m*m*m*x*x)) * SQR(1. + ((2.*m*m)/x));
    }
    if (kappa == 3)
    {
      return (2.*ls/(3.*m*m*m*x*x)) * (SQR(1. - ((m*m)/x)) + (m*m* SQR(p*(kx-p/2.))/ (2.*ls*x*x)));
    }
    if (kappa == 4)
    {
      return (2.*ls/(15.*m*m*m*x*x)) * (1. - (2.*m*m)/x + ((m*m/2.)* ((14.*m*m + (3.*SQR(p*(kx-p/2.)))/ls)/ (x*x))));      
    }
        
    if (kappa == 5)
    {
      return (ls/(3.*m*m*m*x*x)) * (1. - (4.*m*m)/(3.*x) + ((2.*m*m/3.)* ((4.*m*m + (SQR(p*(kx-p/2.)))/ls)/ (x*x))));
    }
    if (kappa == 6)
    {
      return ls/(m*x*x);
    }
    if (kappa == 7)
    {
      return -(1./(12.*m*m*m*x*x)) * (p*p - (4.*SQR(p*(kx-p/2.)))/x + (4.*ls*ls* (4.*m*m + (SQR(p*(kx-p/2.)))/ls)/ (x*x)));
    }
    return 0;
}


double integrand(double* vec, size_t dim, void* p)
{
    parameters *par = ((parameters*)p);   

    double kx =  vec[0];
    double ky =  vec[1];
    

    double G = Gamma(par->kappa, par -> mc , par -> pt, kx, ky);
    
    //Dipole amplitudes
    double F1 = par -> dipole -> S_k( sqrt(SQR(kx) + SQR(ky)), par -> xbj, FUNDAMENTAL, 1.0);
    double F2 = par -> dipole -> S_k (sqrt(SQR((par -> pt)- kx) + SQR(-ky)), par -> xbj, FUNDAMENTAL, 1.0);
    
        
    return G*F1*F2;
   
}



