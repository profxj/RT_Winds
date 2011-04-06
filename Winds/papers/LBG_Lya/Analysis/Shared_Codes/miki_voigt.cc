#include <time.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_cdf.h>
#include "miki_voigt.hh"

// Now using an algorithm from Michele (Zaghloul et al. 2007)


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
  double Hau,error,ffactor;
  
  double ymax = 10*x.Max();
  double dy   = a_param*0.01;
  double ap   = a_param*a_param;

  //initilize integral
  gsl_integration_workspace * w=gsl_integration_workspace_alloc (1000);
       
  //intialize function
  gsl_function F;
  //>>>MF: This line won't work since you are
  //passing double (VOIGT::*)(double, void*) to double (*)(double, void*)
  //This is a problem with C/C++ interaction.
  //Here http://www.cplusplus.com/forum/general/10435/ there is a way to fix it.
  F.function = &VOIGT::CallIntegrand;

  // zero out
  for (i=0;i<n_x;i++) profile[i] = 0;

  // Numerical Approx
  for (i=0;i<n_x;i++)
    {
      voigt_u = x.Center(i);
      
      if(a_param < 200){
	//compute the inegral term for H(a,u)
	//F.params=&p;
	//gsl_integration_qags(&F,0.,voigt_u,0,1e-8,1000,w,&Hau,&error); 
	gsl_integration_qag(&F,0.,voigt_u,0,1e-8,1000,6,w,&Hau,&error); 
	//compute the factor in the front (2 cases)
	if(a_param < 26.6){
	  //use full function
	  ffactor=exp(-voigt_u*voigt_u)*cos(2*a_param*voigt_u)*exp(a_param*a_param)*erfc(a_param);
	  Hau=ffactor+1.12837917*Hau;
	} else {
	  //use series expansion
	  ffactor=exp(-voigt_u*voigt_u)*cos(2*a_param*voigt_u)/(1.7724538*a_param)*
	    (1.-(0.5/pow(a_param,2.))+(0.75/pow(a_param,4.))-(1.875/pow(a_param,6.))+
	     (6.5625/pow(a_param,8.))-(29.53125/pow(a_param,10.))+(162.421875/pow(a_param,12.)));
	  Hau=ffactor+1.12837917*Hau;
	}
      } else {
	//for large a, use approximate equation
	Hau=a_param/(a_param*a_param+voigt_u*voigt_u)/1.7724538;
      }
      profile[i] = Hau;
      // printf("# x, H(a,x) =  %.3e %.3e \n",v,profile[i]);
    }
  
  // printf("# x, H(a,x) =  %.3e %.3e \n",v,profile[i]);
  
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


static double VOIGT::CallIntegrand(double y,void * v)
{
  CCallbackHolder * h = static_cast<CCallbackHolder*>(v);
  return h->cls->Integrand(h->data);
}

double VOIGT::Integrand(double y,void * params)
{
  double f=exp(-(voigt_u*voigt_u-y*y))*sin(2*a_param*(voigt_u-y));
  return f;
}
