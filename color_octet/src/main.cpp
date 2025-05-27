
#include "amplitudelib.hpp"
#include "def.hpp"
#include "tools.hpp"
#include "fourier.h"
#include <string>      
#include <iostream>  
#include <iomanip>   
#include <fstream>
#include <cmath>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>    
#include <gsl/gsl_errno.h>
#include <sstream>



using namespace std;




int main(int argc, char* argv[])  
{
    
   string datafile="./data/proton/mve.dat"; 

    
    // Read data
    AmplitudeLib inst(datafile);  
    inst.SetOutOfRangeErrors(false);
   
    // Compute the J0 transfer
     gsl_set_error_handler_off();
     init_workspace_fourier(1000);
     set_fourier_precision(1.0e-12,1.0e-12);
   
  

    double pt = StrToReal(argv[1]);
    
    double mjpsi = 3.097; //GeV
    double ss = 8160.0; //GeV 
    double y = 2.25;
    double mt = std::sqrt(pt*pt + mjpsi*mjpsi); 
    double xbj =(mt /ss)*std::exp(-y);
    double Q = mjpsi;
    
    

    parameters param;   
    param.xbj= xbj;
    param.dipole = &inst;  
    param.Nc = 3.; 
    param.mc = 1.5;   //mass c quark in GeV
    param.kappa = 6;
    param.pt = pt;

    inst.InitializeInterpolation(param.xbj);
    
    inst.SetFTMethod(ACC_SERIES);
    
    

    // Setup random number generator
    const gsl_rng_type *T;  
    gsl_rng *r;    
    gsl_rng_env_setup ();  
    T = gsl_rng_default;  
    r = gsl_rng_alloc (T); 


    // Function to be integrated for GSL
    gsl_monte_function F;
    F.f = &integrand;
    F.dim = 2;
    F.params = &param;

    
    double max = 99;
    double low[] = {-max,-max};   
    double up[] = {max,max};

    size_t calls =1e4; // Number of function calls per iteration
    
    double res, err;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (F.dim);

    // Quick warmup
    gsl_monte_vegas_integrate (&F, low, up, F.dim, calls/10, r, s, &res, &err);
	//cout << "warmup done" << endl;
   
   int iter=0;
   
    do
      {
        gsl_monte_vegas_integrate (&F, low, up, F.dim, calls, r, s, &res, &err);
        //cout << "# result = " << res << " relative uncertainty " << err/res << " chi^2/dof " <<  gsl_monte_vegas_chisq (s) << endl;
        iter=iter+1;

      }
    while (abs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5 and iter < 5); 
    //or iter < 4);
    
    

    cout << param.pt << " " << res <<  " re. uncert. " << err/res << " chi^2/dof " <<  gsl_monte_vegas_chisq (s) << " k: " << param.kappa << " E: " << ss << " y: " << y << endl;


    gsl_monte_vegas_free (s);


    return 0;
}
    
