#include "def.hpp"
#include "amplitudelib.hpp"
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

double besselFunction(double nu, double x)
{   
    gsl_sf_result result;

    int gsl_status1 = gsl_sf_bessel_Knu_e(nu, x, &result);
    
    if (gsl_status1 != GSL_SUCCESS) 
    {
        return 0.0; 
    }
  
    return result.val;

}

double Gamma(int kappa, double m, double p,  double lr, double lrp, double rx, double rpx, double ry,double rpy)
{ 
    double den1 = (-2.*m/(2.*lr))*besselFunction(1.,lr*m);
    double den2 = (-2.*m/(2.*lrp))*besselFunction(1.,lrp*m);
    double den3 = (1. / m )*(  besselFunction(1.,lr*m)/(2.*lr) + (m/4.)*(-besselFunction(0.,m*lr)-besselFunction(2.,m*lr)) );
    double den4 = (1. / m )*(  besselFunction(1.,lrp*m)/(2.*lrp) + (m/4.)*(-besselFunction(0.,m*lrp)-besselFunction(2.,m*lrp)) );
    
    if (kappa == 1)
    {
     return (2./3.) * m * besselFunction(0.,lr*m) * besselFunction(0.,lrp*m); 
    }
    if (kappa == 2)
    {
        return (((rx*rpx)+(ry*rpy))/m)*den1*den2;
        
    }
    if (kappa == 3)
    {
        return (((rx*rpx)+(ry*rpy))/(3.*m*m*m)) *( ( den1*den2) + ( 2.*m*m * den3*den2) + ( 2.*m*m * den1*den4 ) + ( 4.*m*m*m*m * den3*den4));
    }
    if (kappa == 4)
    {
       return ((2.*((rx*rpx)+(ry*rpy)))/(3.*m*m*m))*( (den1*den2) - (m*m*( den1*den4 + den2*den3)) + (((m*m*m*m) + (m*m*p*rx*p*rpx/(2.*((rx*rpx)+(ry*rpy)))) )*(den3*den4))  );
       
    }
    
    if (kappa == 5)
    {
       return ((2.*((rx*rpx)+(ry*rpy)))/(15.*m*m*m))*( (den1*den2) - (m*m*( den1*den4 + den2*den3)) + ((7.*m*m*m*m) + (3.*m*m*p*rx*p*rpx/(2.*((rx*rpx)+(ry*rpy))) )*(den3*den4)));
    }
    
    return 0;
}


int integrand_cuba(const int *ndim, const double x[],
  const int *ncomp, double f[], void *userdata)
{
    parameters *par = ((parameters*)userdata);   //typecast

    if (*ndim != 6) {
        cerr << "Unknown dimension! " << endl;
        exit(1);
    };

    const double MAX=99;
    const double MIN=-MAX;
    double[7] vec;
    for (int i=0; i<6; i++)
    {
        vec[i] = MIN + x[i]*(MAX-MIN);
    }

    double jacobian = std::pow(MAX-MIN,6);

    return integrand(vec, 6, userdata) * jacobian;
   
}
double integrand(double* vec, size_t dim, void* p)
{
    parameters *par = ((parameters*)p);   //typecast

    double rx = vec[0];
    double rpx = vec[1];
    double Deltax = vec[2];
    double ry = vec[3];
    double rpy = vec[4];
    double Deltay = vec[5];

    double lr =   std::sqrt ( SQR(rx) + SQR(ry));   
    double lrp =  std::sqrt ( SQR(rpx) + SQR(rpy));
    double l1 =  std::sqrt ( SQR(0.5*(rx + rpx) - Deltax) + SQR(0.5*(ry + rpy) - Deltay));   
    double l2 =  std::sqrt ( SQR(0.5*(rx + rpx) + Deltax) + SQR(0.5*(ry + rpy) + Deltay));
    double l3 = std::sqrt ( SQR(0.5*(rx - rpx) - Deltax) + SQR(0.5*(ry - rpy) - Deltay));   
    double l4 =  std::sqrt ( SQR(0.5*(rx - rpx) + Deltax) + SQR(0.5*(ry - rpy) + Deltay));

    double Dr = par->dipole->S(lr, par->xbj);          
    double Drp = par->dipole->S(lrp, par->xbj);
    double D1 = par->dipole->S(l1, par->xbj);
    double D2 = par->dipole->S(l2, par->xbj);
    double D3 = par->dipole->S(l3, par->xbj);
    double D4 = par->dipole->S(l4, par->xbj);


    

    double G = Gamma(par->kappa, par -> mc , par -> pt, lr , lrp, rx,rpx, ry,rpy);

  
    
    double Q =  (std::log((D1 * D2) / (D3 * D4)) / (std::log((Dr * Drp) / (D3 * D4)))) * (-(Dr * Drp) + (D3 * D4)); 
 
    if (std::isnan(Q) or isinf(Q)) {
         //cout << "NaN with D1=" << D1 << " D2=" << D2 <<" D3=" << D3 << " D4="<< D4  << "Dr: " << Dr << "Drp: " << Drp << endl;
        return 0;
    }
    
    
    
    return std::cos(par->pt * Deltax)*Q*G;
    //return Q*G;
    
    
}



