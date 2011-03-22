#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "locate_array.hh"
#include "anvoigt.hh"
#include "spectrum.hh"
#include "model_lbg.hh"
#include "lines.hh"
#include "physical_constants.hh"

extern gsl_rng      *rangen;    // random number generator
extern int verbose;             // output parameter
extern MODEL model;
extern SPECTRUM spectrum;
extern double D_obs[3], cos_obs[4];
extern LINES lines;
ANVOIGT voigt;


struct PHOTON {
  int    ic;
  double r[3], D[3];
  double lam, xloc, lloc;
  double E_p;
  int    escaped;
};

double stepfrac = 0.7;

//--------------------------------------------
// Calculate a spectrum using Monte Carlo
//--------------------------------------------
void Run_Monte_Carlo(double n_photons)
{
  // local variables
  int    i, scatter, l_scat;
  int l, flg_resonance[50];
  double tau_r, tau_x, rad, r_sq, cover;
  double l_step, step, max_step;
  double Pe, lobs, r_rot[3];
  PHOTON p;

  // functions to call
  void MPI_Average_Array(double *, int);
  double P_esc(double, double*, int);
  double Get_vdotD(double,double*, double*);  // Vector only
  void Emit(PHOTON*, double);
  void Line_Scatter(PHOTON&, double, double*);
  void Get_Rotated_Coords(double *, double *);

  // set the start timer 
  time_t start_tp,end_tp;
  time(&start_tp);

  // send the photons
  for (i=0;i<n_photons;i++)
  {
    Emit(&p, model.r_emit);
    // Energy per photon
    p.E_p = model.L_tot/(1.0*n_photons); 

    // add in straight light
    Pe   = P_esc(p.lloc,p.r, -1);  // This one need not be unity
    // printf(" \n", i);
    lobs = p.lloc*(1 - Get_vdotD(model.r_emit, p.r,D_obs)/C_LIGHT);
    Get_Rotated_Coords(p.r,r_rot);
    spectrum.Count(1,lobs,abs(r_rot[0]),abs(r_rot[1]),p.E_p*Pe);

    // if ((i % 5000) == 0)  printf("i = %d \n", i);

    // Did we escape right out?
    if (Pe > 0.99999) continue;
    // printf("# Got here \n");

    // Reset the flag
    for (l=0;l<50;l++) flg_resonance[l] = 0;

    // propogate until escaped (unless nothing stands in the way)
    while (1)
    {
      // locate zone and escape if so
      //p.ic = model.Locate_Zone(p.r);
      //if (p.ic < 0) {p.escaped = 1; break; }

      // photon wavelength in local frame
      rad = sqrt(p.r[0]*p.r[0] + p.r[1]*p.r[1] + p.r[2]*p.r[2]);
      p.lloc = p.lam/(1 - Get_vdotD(rad, p.r, p.D)/C_LIGHT);

      // Variable step size (kpc)
      if (rad < 2.0) step = 1e-4;
      else {
	if (rad > 10.0) step = 0.1; else step = 0.01;
      }
      max_step = step;

      scatter = -1;
      l_scat = -1;

      // calculate random step size to each possible line scatter
      if (rad > model.r_inner) { 
	//l_step = VERY_LARGE_NUMBER;
      	//int l_scat = 0;
	for (int j=0;j<lines.n();j++)
	  {
	    p.xloc = (p.lloc/lines.lambda(j) - 1)*C_LIGHT/model.v_interact();
	    // In resonance
	    if ((p.xloc*p.xloc) <  1 && flg_resonance[j] == 0) {
	      cover = model.Covering(rad);
	      flg_resonance[j] = 1;
	      if (gsl_rng_uniform(rangen) < cover) {step=1e-4; scatter=1; l_scat=j;}
	    }
	  }
      }
      
      // take the step
      p.r[0] += p.D[0]*step;
      p.r[1] += p.D[1]*step;
      p.r[2] += p.D[2]*step;

      // Did we escape?
      r_sq = p.r[0]*p.r[0] + p.r[1]*p.r[1] + p.r[2]*p.r[2];
      if (r_sq > model.r_outer*model.r_outer) {break;}
      
      // Do line scatter if necessary
      if (scatter == 1) 
      {
	// scatter the photon
	double vscat[3];
	Line_Scatter(p, lines.lambda(l_scat), vscat);

	for (l=0;l<50;l++) flg_resonance[l] = 0;

	// count scattered light
	lobs = p.lloc*(1 - Get_vdotD(rad, p.r,D_obs)/C_LIGHT);
	Pe   = P_esc(p.lloc,p.r,l_scat)*lines.P_scat(l_scat);  // Check to come into resonance with redder lines

	Get_Rotated_Coords(p.r,r_rot);
	spectrum.Count(1,lobs,abs(r_rot[0]),abs(r_rot[1]),p.E_p*Pe); 

	// count branching lines
	double dvdp = vscat[0]*D_obs[0] + vscat[1]*D_obs[1] + vscat[2]*D_obs[2];
	for (int j=0;j<lines.n_branch(l_scat);j++) 
	{
	  lobs = lines.blam(l_scat,j)*(1 - Get_vdotD(rad,p.r,D_obs)/C_LIGHT);
	  lobs = lobs*(1 - dvdp/C_LIGHT);
	  Pe   = lines.bprob(l_scat,j);
	  spectrum.Count(1,lobs,r_rot[0],r_rot[1],p.E_p*Pe); 
	}
	
	// see if photon is branched away
	double r1 = gsl_rng_uniform(rangen);
	if (r1 > lines.P_scat(l_scat)) break;
      }

      
    }	

    // debug, for direct count
    //spectrum.Count(1,p.lam,p.r[0],p.r[1],p.E_p); 

    // counter
    //  if (verbose)
    // printf5B("%d\n",i);

    //if (i % (int)(n_photons/10) == 0) printf("count %d\n",i);
  }
   
  // reduce and print spectrum
  spectrum.MPI_Average_All();
  spectrum.Normalize();
  if (verbose) spectrum.Print_Spectrum();
    
  // calculate the elapsed time 
  time(&end_tp);
  float time_wasted=difftime(end_tp,start_tp)/60.0;
  if (verbose)
    printf("#\n# CALCULATION took %.3f minutes (%.2f hours)\n",
	   time_wasted,time_wasted/60.0);
  
}

//   printf("%12d; %12.4e %12.4e %12.4e | %12.4e %12.4e %12.4e | %12.4e %12.4e\n",
//     p.ic, p.r[0],p.r[1],p.r[2],p.D[0],p.D[1],p.D[2],p.lam, p.lloc);


double Get_vdotD(double rad, double *r, double *D)
{
  double vel = model.Velocity(rad);
  double vd = vel*(D[0]*r[0] + D[1]*r[1] + D[2]*r[2])/rad;  
  return vd;
}


void Emit(PHOTON *p, double rphot)
{
  // Get initial positions and direction
  double lam_emit;

  double phi_core = 2*PI*gsl_rng_uniform(rangen); 
  double cosp_core  = cos(phi_core);
  double sinp_core  = sin(phi_core);
  double cost_core  = 1 - 2.0*gsl_rng_uniform(rangen);
  double sint_core  = sqrt(1-cost_core*cost_core);
  // real coordinates
  double remit = rphot * gsl_rng_uniform(rangen);
  p->r[0] = remit*sint_core*cosp_core;
  p->r[1] = remit*sint_core*sinp_core;
  p->r[2] = remit*cost_core;
  // Initial isotropic emission
  double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
  double phi = 2.0*PI*gsl_rng_uniform(rangen);
  double sin_theta = sqrt(1 - mu*mu);
  p->D[0] = sin_theta*cos(phi);
  p->D[1] = sin_theta*sin(phi);
  p->D[2] = mu;

  // local and observer frame wavelength
  lam_emit = model.l_start + (model.l_stop-model.l_start)*gsl_rng_uniform(rangen);
  p->lloc  = lam_emit;
  double vdotD = Get_vdotD(remit, p->r, p->D);  // This is always zero
  p->lam  = p->lloc*(1 - vdotD/C_LIGHT);

  p->escaped = 0;

}


void Line_Scatter(PHOTON &p, double lam, double *uvec)
{
  double u0,u1,u2, R10, R11;

  p.xloc = (p.lloc/lam - 1)*C_LIGHT/model.v_interact();

  // Offset in radius (only proper along radial, but hopefully close enough)
  double rad = sqrt(p.r[0]*p.r[0] + p.r[1]*p.r[1] + p.r[2]*p.r[2]);
  double dvdr = model.DvDr(rad);
  double step = (gsl_rng_uniform(rangen)-0.5) * model.v_doppler() / dvdr;  // kpc
  if (dvdr > 0) {
    p.r[0] += p.D[0]*step;
    p.r[1] += p.D[1]*step;
    p.r[2] += p.D[2]*step;
  }

  // get three velocity components of scatterer
  double temp = 0.5*M_PROTON*pow(model.v_doppler(),2)/K_BOLTZ;
  double apar = 4.7e-3*sqrt(1e4/temp);
  u0 = voigt.Sample_U(p.xloc,apar);
  R10 =  gsl_rng_uniform(rangen);
  R11 =  gsl_rng_uniform(rangen);
  u1 = sqrt(-1.0*log(1.0-R11))*cos(2*PI*R10); 
  u2 = sqrt(-1.0*log(1.0-R11))*sin(2*PI*R10); 
  
  // parallel component
  uvec[0] =  u0*p.D[0];
  uvec[1] =  u0*p.D[1];
  uvec[2] =  u0*p.D[2];
  // perpindicular component 1
  uvec[0] +=    u1*p.D[1];
  uvec[1] += -1*u1*p.D[0];
  // perpindicular component 2
  uvec[0] += u2*p.D[0]*p.D[2];
  uvec[1] += u2*p.D[1]*p.D[2];
  uvec[2] += -1*u2*(p.D[0]*p.D[0] + p.D[1]*p.D[1]);
    
  //  magnitude of Doppler shift into scatterer frame
  double vd_inc =  (uvec[0]*p.D[0] + uvec[1]*p.D[1] + uvec[2]*p.D[2]);
  
   // choose new isotropic direction
  double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
  double phi = 2.0*PI*gsl_rng_uniform(rangen);
  double sin_theta = sqrt(1 - mu*mu);
  p.D[0] = sin_theta*cos(phi);
  p.D[1] = sin_theta*sin(phi);
  p.D[2] = mu;
  
  // magnitude of shift out of scatterer frame
  double vd_out =(uvec[0]*p.D[0] + uvec[1]*p.D[1] + uvec[2]*p.D[2]);
  
  // apply the doppler shift in and out
  p.xloc = p.xloc - vd_inc + vd_out; 
  
  // go back to wavelength
  p.lloc = lam*(1 + p.xloc*model.v_interact()/C_LIGHT);

  // now get change in observer frame wavelength
  rad = sqrt(p.r[0]*p.r[0] + p.r[1]*p.r[1] + p.r[2]*p.r[2]);
  double vdotD = Get_vdotD(rad, p.r, p.D);
  p.lam  = p.lloc*(1 - vdotD/C_LIGHT);
  
  // set proper magnitude of uvec
  uvec[0] *= model.v_doppler();
  uvec[1] *= model.v_doppler();
  uvec[2] *= model.v_doppler();
  
  // check for weirdness
  if (isnan(p.lam)) printf("Snan %e %e %e %e %e\n",p.lloc,p.xloc,u2,R11,rad);

}

///////////////////////////////////// /////////////////////////////////////
//  Escape probability
double P_esc(double lloc, double *r, int l_scat)
{
  double xloc, TOT_ESC, r_sq, step;
  int scatter,  flg, j, l;

  int n_lines = lines.n();

  TOT_ESC = 1.;

  // See if we escape in this line (Doppler allows for multiple scatterings)

  // See if there is a transition to the red
  //  This code was only proper when there was no Doppler velocity to consider
  /// flg = 0;
  /// if(l_scat >= 0) {
  ///   for (j=0; j<n_lines;j++) {
  ///     if (lines.lambda(j) > lines.lambda(l_scat) && j != l_scat) {flg=1;} // Yes there is
  ///   }
  ///   if (flg == 0) return TOT_ESC;  // No there isn't
  /// }
  // printf("# Got here\n");
  
  double xr[3];
  xr[0] = r[0];
  xr[1] = r[1];
  xr[2] = r[2];

  double rad = sqrt(xr[0]*xr[0] + xr[1]*xr[1] + xr[2]*xr[2]);
  double lam = lloc*(1 - Get_vdotD(rad, xr, D_obs)/C_LIGHT);

  int  flg_resonance[50];
  for (l=0;l<50;l++) flg_resonance[l] = 0;

  while (1)
    {
      
      // photon wavelength in local frame
      rad = sqrt(xr[0]*xr[0] + xr[1]*xr[1] + xr[2]*xr[2]);
      lloc = lam/(1 - Get_vdotD(rad, xr, D_obs)/C_LIGHT);
      
      // Variable step size (kpc)
      if (rad < 2.0) step = 1e-4;
      else {
	if (rad > 10.0) step = 0.1; else step = 0.01;
      }
      
      scatter = -1;
      l_scat = -1;
      
      // calculate random step size to each possible line scatter
      if (rad > model.r_inner) { 
	//l_step = VERY_LARGE_NUMBER;
      	//int l_scat = 0;
	for (j=0;j<lines.n();j++)
	  {
	    xloc = (lloc/lines.lambda(j) - 1)*C_LIGHT/model.v_interact();
	    // In resonance, 
	    if ((xloc*xloc) <  1 && flg_resonance[j] == 0) {
	      TOT_ESC *= (1-model.Covering(rad));
	      flg_resonance[j] = 1;
	      if(j == l_scat) printf("# Got here P\n");
	    }
	  }
      }
      
      // take the step
      xr[0] += D_obs[0]*step;
      xr[1] += D_obs[1]*step;
      xr[2] += D_obs[2]*step;

      // Did we escape?
      r_sq = xr[0]*xr[0] + xr[1]*xr[1] + xr[2]*xr[2];
      if (r_sq > model.r_outer*model.r_outer) return TOT_ESC;
    }
      
}


void Get_Rotated_Coords(double *r, double *f)
{
  double r0 = r[0] - model.r_outer;
  double r1 = r[1] - model.r_outer;
  double r2 = r[2] - model.r_outer;

  // apply rotation matrix 
  f[0] = cos_obs[0]*cos_obs[2]*r0 - cos_obs[3]*r1 + cos_obs[1]*cos_obs[2]*r2;
  f[1] = cos_obs[0]*cos_obs[3]*r0 + cos_obs[2]*r1 + cos_obs[1]*cos_obs[3]*r2;
  f[2] = -1*cos_obs[1]*r0 +  cos_obs[0]*r2;

  f[0] = f[0] + model.r_outer;
  f[1] = f[1] + model.r_outer;
  f[2] = f[2] + model.r_outer;
  
}


//------------------------------------------------------------
// General MPI helper functions
//------------------------------------------------------------
void MPI_Average_Array(double *array, int n_el)
{
  double *new_ptr;
  int i;          
  
  // allocate the memory for new pointer
  new_ptr = new double[n_el];
  if (new_ptr == NULL) printf("YO! MPI PROBLEM\n");
  // zero out array
  for (i=0;i<n_el;i++) new_ptr[i] = 0;
  MPI_Allreduce(array,new_ptr,n_el,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  // put back into place
  for (i=0;i<n_el;i++) array[i] = new_ptr[i];
  // free up the memory
  delete new_ptr;

  // get # of processors
  int size;
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  // divide by n_procs
  for (int i=0;i<n_el;i++) array[i] = array[i]/size;
 }


