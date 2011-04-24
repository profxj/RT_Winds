#include "locate_array.hh"
#include "physical_constants.hh"
#include "cdf.hh"
#include <gsl/gsl_integration.h>


class VOIGT
{
  double *profile;
  double a_param, voigt_u;
  LOCATE_ARRAY x;
  int n_x;
  CDF *emit;
  static double CallIntegrand(double,void*); // Michele code only
  double Integrand(double,void*); // Michele code only

public:
  
  void New(int n, double,  double);
  void Compute_Profile();
  void Print();
  double Profile(double);
  double Scatter_Velocity(double);

};

// Miki stuff
struct voigt_params{ 
  double a; 
  double u;
};

//Integral part of the Voig function in the approximation by zaghloul et al. 2007
double miki_voigt (double , void *); 
