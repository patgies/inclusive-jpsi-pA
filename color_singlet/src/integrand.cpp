#include "def.hpp"
#include <cmath>
#include "amplitudelib.hpp"
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_bessel.h>
#include <fstream>
#include <iostream>


// Integrand of the cross section in eq. (3.37) in arXiv:1309.7337v2
//
//

double SQR(double x) 
{ 
    return x*x;  
}

//Bessel functions in appendix (B.13)

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

// Functions for specific quantum state in appendix (B.1.2) 
     
double W(int kappa, double m, double p,  double lr, double lrp, double rx, double rpx, double ry, double rpy)    // p=px!
{ 
    double den1 = (-m*2./(2.*lr))*besselFunction(1,lr*m);
    double den2 = (-m*2./(2.*lrp))*besselFunction(1,lrp*m);
    double den3 = (1. / m )*(  besselFunction(1,lr*m)/(2.*lr) + (m/4.)*(-besselFunction(0,m*lr)-besselFunction(2,m*lr)) );
    double den4 = (1. / m )*(  besselFunction(1,lrp*m)/(2.*lrp) + (m/4.)*(-besselFunction(0,m*lrp)-besselFunction(2,m*lrp)) );
    double fact = ((rx*rpx)+(ry*rpy));
    
    if (kappa == 1)  //3S1
    {
     return (2./3.) * m * besselFunction(0,lr*m) * besselFunction(0,lrp*m); 
    }
    if (kappa == 2)   //3P0
    {
        return (fact/(3.*m*m*m)) *( ( den1*den2) + ( 2.*m*m * den3*den2) + ( 2.*m*m * den1*den4 ) + ( 4.*m*m*m*m * den3*den4));
    }
    if (kappa == 3)   //3P1
    {
       return ((2.*fact)/(3.*m*m*m))*( (den1*den2) - (m*m*( den1*den4 + den2*den3)) + (((m*m*m*m) + (m*m*p*rx*p*rpx/(2.*fact)) )/(den3*den4))  );
    }
    if (kappa == 4)   //3P2
    {
       return ((2.*fact)/(15.*m*m*m))*( den1*den2 - m*m*( den1*den4 + den2*den3) + ((7.*m*m*m*m) + (3.*m*m*p*rx*p*rpx/(2.*fact)) )/(den3*den4));
    }
    
    if (kappa == 5)  //1S0
    {
       return (fact/m) * den1*den2;
    }
    
    return 0;
}


int integrand_cuba(const int *ndim, const double x[],
  const int *ncomp, double f[], void *userdata)
{
    

    if (*ndim != 6) {
        cerr << "Unknown dimension! " << endl;
        exit(1);
    };

   // Polar coordinates
   // Params: r, theta_r, rp, theta_rp, delta, theta_delta
   const double MAXR=20;
   double vec[6];
   vec[0] = MAXR*x[0];
   vec[2] = MAXR*x[2];
    vec[4] = MAXR*x[4];

    vec[1]=2.*M_PI*x[1];
    vec[3]=2.*M_PI*x[3];
    vec[5]=2.*M_PI*x[5];

    double jacobian = vec[0]*vec[2]*vec[4]*std::pow(MAXR,3)*std::pow(2.*M_PI,3);

    f[0]=integrand(vec, 6, userdata)* jacobian;
    //delete[] vec;
    return 0;
   
}
double integrand(double* vec, size_t dim, void* p)
{
    parameters *par = ((parameters*)p);   //typecast

    //Integration variables in coordinate space
    */
    double rx = vec[0]*std::cos(vec[1]);
    double rpx = vec[2]*std::cos(vec[3]);
    double Deltax = vec[4]*std::cos(vec[5]);
    double ry = vec[0]*std::sin(vec[1]);
    double rpy = vec[2]*std::sin(vec[3]);
    double Deltay = vec[4]*std::sin(vec[5]);

    //Dipole lengths
    
    double pr =  std::sqrt ( SQR(par ->pt)+ SQR(rx));
    double prp = std::sqrt ( SQR(par ->pt) + SQR(rpx));
    double lr =   std::sqrt ( SQR(rx) + SQR(ry));   
    double lrp =  std::sqrt ( SQR(rpx) + SQR(rpy));
    double l1 = (0.5) * std::sqrt ( SQR(rx + rpx - 2.0*Deltax) + SQR(ry + rpy - 2.0*Deltay));   
    double l2 = (0.5) * std::sqrt ( SQR(rx + rpx + 2.0*Deltax) + SQR(ry + rpy + 2.0*Deltay));
    double l3 = (0.5) * std::sqrt ( SQR(rx - rpx - 2.0*Deltax) + SQR(ry - rpy - 2.0*Deltay));   
    double l4 = (0.5) * std::sqrt ( SQR(rx - rpx + 2.0*Deltax) + SQR(ry - rpy + 2.0*Deltay));

    
    //Dipoles
    
    double Dr = par->dipole->S(lr, par->xbj);          
    double Drp = par->dipole->S(lrp, par->xbj);
    double D1 = par->dipole->S(l1, par->xbj);
    double D2 = par->dipole->S(l2, par->xbj);
    double D3 = par->dipole->S(l3, par->xbj);
    double D4 = par->dipole->S(l4, par->xbj);


    //Function specific quantum state
    
    double state = W(par->kappa, par -> mc , par -> pt, lr , lrp, rx,rpx, ry,rpy);

    //std::ofstream outputFile("Nan.txt", std::ios::app);
 
 
    //Quadrupole
    
    double fm = std::log( Dr*Drp / (D3*D4));    
    double fl = std::log( D1*D2 / (D3*D4));
    double delta = std::sqrt( fm*fm + (4./9.)*  fl * std::log( D1*D2 / (Dr*Drp)) );
    double e1 = std::exp((9./16.) * delta);
    double e2 = std::exp( (-9./16.) * delta);
    double e3 = std::exp( ((-9./16.) * fm) + ((1./8.) * fl) );            //CF=4/3
 
    
    double Q = (Dr * Drp) *  ( ( (delta + fm)/(2.*delta) - (fl/delta))*e1 + ((delta-fm)/(2.*delta) + fl/delta)*e2 ) * e3; 
    
   
    
    if (std::isnan(Q) or isinf(Q)) {
        //cerr << "NaN with D1=" << D1 << " D2=" << D2 <<" D3=" << D3 << " D4="<< D4  << endl;
        //outputFile << "NaN detected in Q calculation." << D3<<""<< D4<< std::endl;
        return 0;
    }
    //outputFile.close();
    
    
    return std::cos(par->pt * Deltax)* (Q - Dr*Drp) * state;

    
}



