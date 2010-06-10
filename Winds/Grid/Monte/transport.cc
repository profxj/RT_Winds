#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "locate_array.hh"
#include "voigt.hh"
#include "spectrum.hh"
#include "model.hh"
#include "lines.hh"

extern gsl_rng      *rangen;    // random number generator
extern int verbose;             // output parameter
extern MODEL model;
extern SPECTRUM spectrum;
extern VOIGT voigt;
extern double D_obs[3], cos_obs[4];
extern LINES lines;


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
  int    i, scatter;
  double tau_r, tau_x, vdotD;
  double l_step, step, max_step;
  double Pe, lobs, r_rot[3];
  PHOTON p;

  // functions to call
  void MPI_Average_Array(double *, int);
  double P_esc(double, double*);
  void Emit(PHOTON*);
  void Line_Scatter(PHOTON&, double);
  void Get_Rotated_Coords(double *, double *);

  // set the start timer 
  time_t start_tp,end_tp;
  time(&start_tp);

  // send the photons
  for (i=0;i<n_photons;i++)
  {
    Emit(&p);
    // Energy per photon
    p.E_p = model.L_tot/(1.0*n_photons); 

    // add in straight light
    Pe   = P_esc(p.lloc,p.r);
    lobs = p.lloc*(1 - model.vdotD(p.ic,D_obs)/C_LIGHT);
    Get_Rotated_Coords(p.r,r_rot);
    spectrum.Count(1,lobs,r_rot[0],r_rot[1],p.E_p*Pe);

    // propogate until escaped (unless nothing stands in the way)
    if (Pe < 0.99)
    while (1)
    {
      // locate zone and escape if so
      p.ic = model.Locate_Zone(p.r);
      if (p.ic < 0) {p.escaped = 1; break; }

      // photon wavelength in local frame
      p.lloc = p.lam/(1 - model.vdotD(p.ic, p.D)/C_LIGHT);

      // default step size
      max_step = model.dx*stepfrac;
      step = max_step;
      scatter = -1;

      // calculate random step size to each possible line scatter
      l_step = VERY_LARGE_NUMBER;
      int l_scat = 0;
      for (int j=0;j<lines.n();j++)
      {
	p.xloc = (p.lloc/lines.lambda(j) - 1)*C_LIGHT/model.v_doppler(p.ic);
	tau_r   =  -1.0*log(1.0 - gsl_rng_uniform(rangen));
	tau_x = model.line_opacity(p.ic)*lines.cs(j)*voigt.Profile(p.xloc)*KILOPARSEC;
	double this_step = tau_r/tau_x;
	if (tau_x == 0) this_step = VERY_LARGE_NUMBER;
	if (this_step < l_step) {l_step = this_step; l_scat = j; }
      }
      
      // see if this a line step
      if (l_step < max_step) {step = l_step; scatter = 1; }

      // take the step
      p.r[0] += p.D[0]*step;
      p.r[1] += p.D[1]*step;
      p.r[2] += p.D[2]*step;
      
      // Do line scatter if necessary
      if (scatter == 1) 
      {
	// scatter the photon
	Line_Scatter(p, lines.lambda(l_scat));

	// count scattered light
	lobs = p.lloc*(1 - model.vdotD(p.ic,D_obs)/C_LIGHT);
	Pe   = P_esc(p.lloc,p.r)*lines.P_scat(l_scat);

	Get_Rotated_Coords(p.r,r_rot);
	spectrum.Count(1,lobs,r_rot[0],r_rot[1],p.E_p*Pe); 

	// count branching lines
	for (int j=0;j<lines.n_branch(l_scat);j++) 
	{
	  lobs = lines.blam(l_scat,j)*(1 - model.vdotD(p.ic,D_obs)/C_LIGHT);
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




void Emit(PHOTON *p)
{
  // Get initial positions and direction
  double r_emit[3], lam_emit;

  int cell = model.Emit(r_emit, &lam_emit);
  p->ic   = cell;
  p->r[0] = r_emit[0];
  p->r[1] = r_emit[1];
  p->r[2] = r_emit[2];
  
  // local and observer frame wavelength
  p->lloc  = lam_emit;
  double vdotD = model.vdotD(p->ic, p->D);
  p->lam  = p->lloc*(1 - vdotD/C_LIGHT);

  // initial isotropic emission
  double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
  double phi = 2.0*PI*gsl_rng_uniform(rangen);
  double sin_theta = sqrt(1 - mu*mu);
  p->D[0] = sin_theta*cos(phi);
  p->D[1] = sin_theta*sin(phi);
  p->D[2] = mu;

  p->escaped = 0;

}


void Line_Scatter(PHOTON &p, double lam)
{
  double u0,u1,u2, R10, R11;
  double uvec[3];

  p.xloc = (p.lloc/lam - 1)*C_LIGHT/model.v_doppler(p.ic);

  // get three velocity components of scatterer
  u0 = voigt.Scatter_Velocity(p.xloc);
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
  p.lloc = lam*(1 + p.xloc*model.v_doppler(p.ic)/C_LIGHT);

  // now get change in observer frame wavelength
  double vdotD = model.vdotD(p.ic, p.D);
  p.lam  = p.lloc*(1 - vdotD/C_LIGHT);
  
  if (isnan(p.lam)) printf("Snan %e %e %e %e\n",p.lloc,p.xloc,u2,R11);

}

double P_esc(double lloc, double *r)
{
  double MAX_TAU = 20;
  double xloc;
  int ind;
  
  double dr = model.dx*stepfrac;
  
  double xr[3];
  xr[0] = r[0];
  xr[1] = r[1];
  xr[2] = r[2];

  int n_lines = lines.n();

  // get location
  ind = model.Locate_Zone(xr);
  if (ind < 0) return 1;
  double lam = lloc*(1 - model.vdotD(ind,D_obs)/C_LIGHT);


  // integrate along a ray to the observer
  double tau = 0;
  for (double z = 0;z < model.xmax; z+= dr)
  {
    lloc = lam/(1 - model.vdotD(ind,D_obs)/C_LIGHT);

    for (int j=0;j<n_lines;j++) 
    {
      xloc = (lloc/lines.lambda(j) - 1)*C_LIGHT/model.v_doppler(ind);
      tau += model.opac[ind]*lines.cs(j)*(dr*KILOPARSEC)*voigt.Profile(xloc);
    }

    if (isnan(tau))  
      printf("tnan %e %e %e %e %e\n",model.opac[ind],lloc,
	     lam, xloc,voigt.Profile(xloc));
	     
    if (tau > MAX_TAU) break;
    
    // take the step
    xr[0] += D_obs[0]*dr;
    xr[1] += D_obs[1]*dr;
    xr[2] += D_obs[2]*dr;

    // get location
    ind = model.Locate_Zone(xr);
    if (ind < 0) break;
  }
    
  double P;
  if (tau >= MAX_TAU) P = 0;
  else P = exp(-tau);
  return P;
}


void Get_Rotated_Coords(double *r, double *f)
{
  double r0 = r[0] - model.xmax/2.0;
  double r1 = r[1] - model.xmax/2.0;
  double r2 = r[2] - model.xmax/2.0;

  // apply rotation matrix 
  f[0] = cos_obs[0]*cos_obs[2]*r0 - cos_obs[3]*r1 + cos_obs[1]*cos_obs[2]*r2;
  f[1] = cos_obs[0]*cos_obs[3]*r0 + cos_obs[2]*r1 + cos_obs[1]*cos_obs[3]*r2;
  f[2] = -1*cos_obs[1]*r0 +  cos_obs[0]*r2;

  f[0] = f[0] + model.xmax/2.0;
  f[1] = f[1] + model.xmax/2.0;
  f[2] = f[2] + model.xmax/2.0;
  
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


