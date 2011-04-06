#include <time.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "voigt.hh"


void VOIGT::New(int n, double xmax, double a)
{
  n_x = n;
  a_param = a;
  x.Uniform_New(-1*xmax,2*xmax/n_x,n_x);
  profile = new double[n_x];
  emit = new CDF[n_x];
  for (int i=0;i<n_x;i++) emit[i].New(n_x);
  Compute_Profile();
}

void VOIGT::Compute_Profile()
{
  int i;
  double y,v,e,expy,expv;
  double xsq, c, q, pic,SQRT_PI;
  
  double ymax = 10*x.Max();
  double dy   = a_param*0.01;
  double ap   = a_param*a_param;

  // zero out
  for (i=0;i<n_x;i++) profile[i] = 0;

  // // do integral over yfrom -infinity to +infinity
  // for (y = -1*ymax;y<ymax;y+=dy)
  // {
  //   expy = exp(-1*y*y);
  //   for (i=0;i<n_x;i++)
  //   {
  //     v = x.Center(i);
  //     profile[i] += expy/((v-y)*(v-y) + ap)*dy;
  //   }
  // }

  // // renormalize
  // for (i=0;i<n_x;i++) profile[i] = a_param/PI*profile[i]/sqrt(PI);


  // Numerical Approx
  SQRT_PI = sqrt(PI);
  for (i=0;i<n_x;i++)
    {
      v = x.Center(i);
      xsq = v*v;
      c   = (xsq - 0.855)/(xsq + 3.42);

      if (c < 0) q = 0;
      else 
	{
	  pic = 5.674*c*c*c*c -9.207*c*c*c + 4.421*c*c + 0.1117*c;
	  q = (1 + 21/xsq)*a_param/PI/(xsq + 1)*pic; 
	}
      
      profile[i] = q*SQRT_PI + exp(-xsq);
      // printf("# x, H(a,x) =  %.3e %.3e \n",v,profile[i]);
    }
  
  // renormalize
  for (i=0;i<n_x;i++) profile[i] = profile[i]/sqrt(PI);

  int j;
  for (i=0;i<n_x;i++) 
  {
    y = x.Center(i);
    for (j=0; j<n_x; j++) 
    {
      v = x.Center(j);
      expv = exp(-1*v*v);
      e = expv/((v-y)*(v-y) + ap);
      emit[i].Put(j,e);
    }
    emit[i].Normalize(); 
  }
}


double VOIGT::Scatter_Velocity(double thisx)
{  
  int i = x.Locate(thisx);
  if (i < 0) return 0;
  int z = emit[i].Sample();
  return x.Sample(z);
}

 
double VOIGT::Profile(double y)
{
  int i = x.Locate(y);
  if (i < 0) return 0;
  if (i >= x.Size()-1) return profile[i];

  // linear interpolation
  double dy = (y - x.Center(i))/(x.Center(i+1) - x.Center(i));
  double val = profile[i] + dy*(profile[i+1] - profile[i]);
  return val;
}

void VOIGT::Print()
{
  for (int i=0;i<n_x;i++)
    printf("%12.5e %12.5e\n",x.Center(i),profile[i]);
}
