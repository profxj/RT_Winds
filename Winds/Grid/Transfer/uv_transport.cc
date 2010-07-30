#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "locate_array.hh"
#include "model.hh"
#include "physical_constants.hh"


extern gsl_rng      *rangen;    // random number generator
extern int verbose;             // output parameter
extern MODEL model;

struct PHOTON {
  int    ic;
  double r[3], D[3];
  double lam, xloc, lloc;
  double E_p;
  int    escaped;
};


//--------------------------------------------
// Calculate a spectrum using Monte Carlo
//--------------------------------------------
void Run_UV_Transport(double n_photons, double L_UV, int n_iter)
{
  void Run_UV_MC(double, double);
  
  for (int i=0;i<n_iter;i++)
  {
    model.Wipe_J_UV();
    model.Compute_Ionization_State();
    model.Compute_Photo_Opacity();
    Run_UV_MC(n_photons, L_UV);
    model.Reduce_J_UV();
    model.Print_J_UV();
  }

}


void Run_UV_MC(double n_photons, double L_UV)
{

  double stepfrac = 0.3;

  // local variables
  int    i, scatter;
  double tau_r, opac, vdotD;
  double l_step, step, max_step;
  double mu, phi, sin_theta;
  PHOTON p;

  // set the start timer 
  time_t start_tp,end_tp;
  time(&start_tp);

  // send the photons
  for (i=0;i<n_photons;i++)
  {
    double xcen = model.xmax/2.0;
    // Emit the photon
    mu  = 1 - 2.0*gsl_rng_uniform(rangen);
    phi = 2.0*PI*gsl_rng_uniform(rangen);
    sin_theta = sqrt(1 - mu*mu);
    double r_emit = 0; // xcen
    p.r[0] = xcen + r_emit*sin_theta*cos(phi);
    p.r[1] = xcen + r_emit*sin_theta*sin(phi);
    p.r[2] = xcen + r_emit*mu; 
    mu  = 1 - 2.0*gsl_rng_uniform(rangen);
    phi = 2.0*PI*gsl_rng_uniform(rangen);
    sin_theta = sqrt(1 - mu*mu);
    p.D[0] = sin_theta*cos(phi);
    p.D[1] = sin_theta*sin(phi);
    p.D[2] = mu;
    p.E_p = L_UV/(1.0*n_photons); 
    p.lam  = 500;
    p.escaped = 0;
    // -- done emit
    
    // propogate until escaped
    while (1)
    {
      // locate zone and escape if so
      p.ic = model.Locate_Zone(p.r);
      if (p.ic < 0) {p.escaped = 1; break; }

      // photon wavelength in local frame
      p.lloc = p.lam*(1 - model.vdotD(p.ic, p.D)/C_LIGHT);
      
      // default step size
      max_step = model.dx*stepfrac;
      step = max_step;
      scatter = -1;

      // calculate random step size to each possible line scatter
      tau_r =  -1.0*log(1.0 - gsl_rng_uniform(rangen));
      opac  = model.photo_opacity(p.ic,p.lloc);
      // debug
      opac = 0;
      l_step = tau_r/opac;
      if (opac == 0) l_step = VERY_LARGE_NUMBER;
      if (l_step < max_step) {step = l_step; scatter = 1; }

      // take the step
      p.r[0] += p.D[0]*step;
      p.r[1] += p.D[1]*step;
      p.r[2] += p.D[2]*step;

      // add in energy to model
      model.J_UV[p.ic] += p.E_p*step; 
      
      // Do line scatter if necessary
      if (scatter == 1) 
      {
	printf("SCAT\n");
	double z = gsl_rng_uniform(rangen);
	if (z < 0.38)
	{
	  // scatter isotropically
	  double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
	  double phi = 2.0*PI*gsl_rng_uniform(rangen);
	  double sin_theta = sqrt(1 - mu*mu);
	  p.D[0] = sin_theta*cos(phi);
	  p.D[1] = sin_theta*sin(phi);
	  p.D[2] = mu;      
	}
	else break;
	
      }	
      
    }
  }
  // calculate the elapsed time 
  time(&end_tp);
  float time_wasted=difftime(end_tp,start_tp)/60.0;
  if (verbose)
    printf("#\n# UV Transport DONE took %.3f minutes (%.2f hours)\n",
	   time_wasted,time_wasted/60.0);
  
}




