#include "amplitudelib.hpp"
#include "def.hpp"
#include <string>      
#include <iostream>    
#include <iomanip>   
#include <cmath>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>    
#include <gsl/gsl_errno.h>
#include "tools.hpp"
#include "cuba.h"


using namespace std;



int main(int argc, char* argv[])  
{
     
    string datafile="./data/proton/mve.dat";   
    
    // Read data
    AmplitudeLib inst(datafile);  
    inst.SetOutOfRangeErrors(false);
   
    gsl_set_error_handler_off();

    double pt = StrToReal(argv[1]);

    double mjpsi = 3.097; // mass vector mesonGeV
    double ss = 8160.0; // sqrt(s) energy in GeV 
    double y = 3.25;     //rapidity
    double mt = std::sqrt(pt*pt + mjpsi*mjpsi); //transverse mass
    double xbj = (mt /ss)*std::exp(-y);     //Bjorken-x target
    

    parameters param;   
    param.xbj= xbj;
    param.dipole = &inst;  
    param.Nc = 3.; 
    param.mc = 1.5;   //mass c quark in GeV
    param.kappa = 1;  //meson in specific quantum state
    param.pt = pt;

    inst.InitializeInterpolation(param.xbj);

    // Setup random number generator
    const gsl_rng_type *T;  
    gsl_rng *r;    
    gsl_rng_env_setup ();  
    T = gsl_rng_default;  
    r = gsl_rng_alloc (T); 


    // Function to be integrated for GSL
    gsl_monte_function F;
    F.f = &integrand;
    F.dim = 6;
    F.params = &param;

    
    double max = 99;
    double low[] = {-max,-max,-max,-max,-max,-max};   
    double up[] = {max,max, max, max, max, max};

    long long int calls = StrToReal(argv[2]);

    double res, err;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (F.dim);


    const int VERBOSE=0;
    int nnew=calls/2, nmin=200; // nnew=10e3
    double flatness=21; //25;
    int neval, fail, nregions;
    double *prob = new double[F.dim];
    Suave(F.dim, 1, integrand_cuba, &param, 1, 1e-4, 0, VERBOSE, 0, calls/2, calls*5, 
                nnew, nmin, flatness, "", NULL, &nregions, &neval, &fail, &res, &err, prob  );
    

    cout << param.pt << " " << res << " " << err << endl; // " relative uncertainty " << err/res << " chi^2/dof " <<  gsl_monte_vegas_chisq (s) << " kappa: " << param.kappa << " E: " << ss  << endl;


    gsl_monte_vegas_free (s);


    return 0;
}
    
