#include "locate_array.hh"
#include "physical_constants.hh"
#include "cdf.hh"

class EMIT_LINE
{
  double *cumul;
  double EW_x,sigma_x;
  LOCATE_ARRAY wave;
  int n_x;

public:
  
  void New(int n, double,  double);
  void Compute_Cumul();
  //  void Print();
  double Profile(double);
  double Scatter_Velocity(double);

};
