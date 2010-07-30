#ifndef _ANVOIGT_H
#define _ANVOIGT_H

#include "locate_array.hh"
#include "physical_constants.hh"
#include "cdf.hh"
#include <gsl/gsl_rng.h>


class ANVOIGT
{
  double u0, eu0;
  gsl_rng  *rangen;    // random number generator

public:
  
  ANVOIGT();
  void   Set_U0(double);
  double Profile(double,double);
  double Sample_U(double,double);

};

#endif
