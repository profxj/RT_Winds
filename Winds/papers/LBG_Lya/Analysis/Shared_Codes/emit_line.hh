#include "locate_array.hh"
#include "physical_constants.hh"
#include "cdf.hh"

class EMIT_LINE
{
  double *wave;
  double EW_x,sigma_x;
  LOCATE_ARRAY cumul;
  int n_x;

public:
  
  void New(int n, double,  double, double, double);
  void Compute_Cumul();
  //  void Print();
  double Get_Wave(double);

};
