#ifndef def_hpp 
#define def_hpp
#include "amplitudelib.hpp"



struct parameters          
{
    AmplitudeLib *dipole;    
    double xbj;
    double Nc; 
    double mc;
    double pt;
};


double integrand(double* vec, size_t dim, void* p);


#endif
