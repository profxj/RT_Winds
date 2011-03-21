#include "locate_array.hh"
#include "physical_constants.hh"
#include "cdf.hh"

class VOIGT
{
  double *profile;
  double a_param;
  LOCATE_ARRAY x;
  int n_x;
  CDF *emit;

public:
  
  void New(int n, double,  double);
  void Compute_Profile();
  void Print();
  double Profile(double);
  double Scatter_Velocity(double);

};
