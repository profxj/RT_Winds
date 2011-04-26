#include <time.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "anvoigt.hh"
#define  SQRT_PI 1.77245


ANVOIGT::ANVOIGT()
{
  u0  = 0;
  eu0 = 1;

  // setup and seed random number generator
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  gsl_rng_default_seed = (unsigned int)time(NULL);
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);
}


void ANVOIGT::Set_U0(double upass)
{
  u0     = upass;
  eu0    = exp(-u0*u0);
}


double ANVOIGT::Profile(double x, double a)
{
  double xsq = x*x;
  double c   = (xsq - 0.855)/(xsq + 3.42);

  double q;
  if (c < 0) q = 0;
  else 
  {
    double pic = 5.674*c*c*c*c -9.207*c*c*c + 4.421*c*c + 0.1117*c;
    q = (1 + 21/xsq)*a/PI/(xsq + 1)*pic; 
  }
 
  double H = q*SQRT_PI + exp(-xsq);
  return H;
}



double ANVOIGT::Sample_U(double x, double a)
{
  double u,th;

  double sgn = 1;
  if (x < 0) {sgn = -1; x = -1*x;}

  double theta0 = atan((u0 - x)/a);
  double denom  = (1 - eu0)*theta0 + (1 + eu0)*PI/2;
  double p0 = (theta0 + PI/2)/denom;
  
  int stop = 0;
  while (!stop) 
  {
    double r1 = gsl_rng_uniform(rangen);
    double r2 = gsl_rng_uniform(rangen);
    double r3 = gsl_rng_uniform(rangen);

    if (r1 < p0) 
    {
      th = -1*PI/2.0 + r2*(theta0 + PI/2);
      u = a*tan(th) + x; 
      if (r3 < exp(-u*u)) stop = 1;
    }
    else 
    {
      th = theta0 + r2*(PI/2 - theta0);
      u = a*tan(th) + x; 
      if (r3 < exp(-u*u)/eu0) stop = 1;
    }
  }

  return sgn*u;
}
