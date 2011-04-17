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

struct CCallbackHolder
{
  VOIGT* cls;
  void* data;
};

struct voig_params{ 
  double a; 
  double u;
};
