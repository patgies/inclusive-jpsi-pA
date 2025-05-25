#ifndef def_hpp 
#define def_hpp
#include "amplitudelib.hpp"

//Parameters used in the quadrupole

struct parameters          
{
    AmplitudeLib *dipole;    
    double xbj; 
    double Nc; 
    double mc;
    int kappa;
    double pt;
};

int integrand_cuba(const int *ndim, const double x[],
  const int *ncomp, double f[], void *userdata);


double integrand(double* vec, size_t dim, void* p);


#endif
