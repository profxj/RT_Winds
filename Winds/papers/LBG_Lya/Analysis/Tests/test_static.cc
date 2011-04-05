#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "locate_array.hh"
#include "voigt.hh"
#include "spectrum.hh"

// monte carlo parameters
double n_photons =  1e5;    // number of photons
double stepsize  = 0.01;   // maximum size of photon step 

// output spectrum parameters
double l_start   =  1205;     // beginning wavelength (Angstroms)
double l_stop    =  1225;     // ending wavelength (Angstroms)
double l_delta   =   0.1;     // wavelength resolution (Angstroms)
double F_cont    =    1;      // continuum flux level
int    n_mu      =    1;      // number of theta bins
int    n_phi     =    1;      // number of phi bins

// grid parameters
double r_inner   =  0.0;      // inner boundary radius, in kpc
double r_outer   = 10.0;      // outer boundary radius, in kpc
double r_emit    =  0.0;      // boundary to emit from
double n_0       =  1e0;     // number density at inner boundary
double n_law     =    2.;      // -1*power law exponent of density law
double v_max     =  4.0e7;    // velocity at outer boundary (cm/s)
double v_min     =  50.*1e5;    // velocity at inner boundary (cm/s)
double v_law     =    1.0;      // power law of velocity profile
double bipolar   =    0;      // degree of bipolarity

double dust_cs     = 3.33e-24;    // dust cross-section
double dust_dens   = 0.;          // density of NH/dust, 1 gives tau=1
double dust_albedo = 0.0;         // ratio of scattering to absorption


// line parameters
// --------------------------
 //number of lines to use
int    n_lines      = 1; 
// line center wavelengths
double lambda_0[]   = {1215.6701};
// line oscillator strengths
double f_lu[]       = {0.4164};
// abundances of element of line
double abun[]       = {1.};  // Solar metallicity and 1/10 down for dust
double metallicity       = 1.0;                 // 1 = Solar
 // lines doppler velocity in cm/s
double v_doppler    =   20*1e5;              

// parameters describing voigt profile
double Dnu   = v_doppler * 2.466e15 / 2.9979e10 ;   
double voigt_a   = 6.265e8 / (4 * 3.14159 * Dnu);   // gamma/(4 pi Dnu)
int    nvoigt    = 2000;
double voigt_x   =  200;  // Allows for over 2000 km/s
VOIGT voigt;


// globals
gsl_rng      *rangen;    // random number generator
int verbose;             // output parameter

//--------------------------------------------
// the main program
//--------------------------------------------
int main(int argc, char **argv)
{
  void Run_Monte_Carlo(char*);

  // initialize MPI for parallelism
  int my_rank, n_procs;
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs );
  if (my_rank == 0) verbose = 1; else verbose = 0;

  // setup and seed random number generator
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  gsl_rng_default_seed = (unsigned int)time(NULL) + my_rank;
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);
  
  // initialize the Voigt profile
  // printf("# Voigt a =  %.3e \n", voigt_a);
  voigt.New(nvoigt,voigt_x,voigt_a);
  // printf("# Finished Voigt!");


  // get photons per processor
  n_photons = n_photons/n_procs;
  if (verbose)
    printf("# Sending %.3e photons on %d procs (%.3e total photons)\n",
	   n_photons,n_procs,n_photons*n_procs);


  // Do the monte carlo calculation
  char outfile[1000];
  sprintf(outfile,"spec.dat");
  Run_Monte_Carlo(outfile);

  // finish up mpi
  MPI_Finalize();

}


//--------------------------------------------
// Function returns the velocity given a radius
//--------------------------------------------
double Get_Velocity(double *x, double r)
{
  //  return v_max*pow(r/r_outer,v_law);
  // return v_min + (r-r_inner)/(r_outer-r_inner) * (v_max-v_min);
  return 0.;  // Static medium
}


//--------------------------------------------
// Function returns the density given a radius
//--------------------------------------------
double Get_Density(double *x, double r)
{
  if (r == 0) return 0;
  if (r < r_inner) return 0;
  //  if (r > r_outer) return 0;
  // double mu = x[2]/r;
  return n_0;
}


//--------------------------------------------
// Calculate a spectrum using Monte Carlo
//--------------------------------------------
void Run_Monte_Carlo(char *outfile)
{
  // local variables
  int i, l, ind, scatter, dust_scatter, count_it, flg_scatter;
  double x, lam_loc, xloc,lam, lam_emit;
  double r[3], D[3];
  double mu,phi,sin_theta;
  double tau_r, tau_x, step, r_sq;
  double vd_inc, vd_out, l_step, d_step;
  double u0,u1,u2, R10, R11, rad, vel;
  double uvec[3];
  double nu_d, cross_sec; 
  double tau_const, atau, NHI, xcr;

  tau_const = 0.0062344885;   // Constants + f :: Needs N_HI and Dnu

  // functions to call
  void MPI_Average_Array(double *, int);
  void Emit(double *r, double *D, double r_inner);
  double Get_Velocity(double*, double);
  double Get_Density(double*, double);

  // set the start timer 
  time_t start_tp,end_tp;
  time(&start_tp);

  // define spectrum counter
  SPECTRUM spectrum;
  double lgrid[3];
  lgrid[0] = l_start;
  lgrid[1] = l_stop;
  lgrid[2] = l_delta;
  spectrum.Init(lgrid,n_mu,n_phi);
  spectrum.Set_Name(outfile);
 
  // Energy per photon
  int n_wave = (l_stop - l_start)/l_delta;
  double E_p = F_cont*n_wave*n_mu*n_phi/n_photons*l_delta;

  // send the photons
  for (i=0;i<n_photons;i++)
  {
    // Get initial positions and direction
    Emit(r,D,r_emit);

    // initial wavelength (Uniform continuum)
    lam = l_start + (l_stop-l_start)*gsl_rng_uniform(rangen);
    lam_emit = lam;
    flg_scatter = 0;

    // printf("#  Photon %d, lambda = %.4e \n", i, lam);

    // propogate until escaped
    while (1)
    {
      // photon wavelength in local frame
      rad = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
      vel = Get_Velocity(r,rad)/C_LIGHT;
      lam_loc = lam*(1 + vel*(r[0]*D[0] + r[1]*D[1] + r[2]*D[2])/rad);

      // printf("#  rad = %.4e, lambda = %.4e \n", rad, lam_loc);
      if (rad == 0) lam_loc = lam;

      // default step size
      step = stepsize;

      // calculate random step size to each possible line scatter
      scatter = -1;
      for (l=0;l<n_lines;l++)
      {
	// x parameter for this line
	xloc = (lam_loc/lambda_0[l] - 1)*C_LIGHT/v_doppler;
	// random optical depth to travel
	tau_r     =  -1.0*log(1 - gsl_rng_uniform(rangen));
	nu_d = (C_LIGHT/lambda_0[l]/ANGS_TO_CM)*(v_doppler/C_LIGHT);
	cross_sec = CLASSICAL_CS*f_lu[l]*voigt.Profile(xloc)/nu_d;
	tau_x     = KPARSEC*Get_Density(r,rad)*abun[l]*metallicity*cross_sec;
	l_step = tau_r/tau_x;
	if (tau_x == 0) l_step = VERY_LARGE_NUMBER;
	if (l_step < step) {step = l_step; scatter = l; }
      }

      // get distance to dust scatter/absorption
      tau_r = -1.0*log(1 - gsl_rng_uniform(rangen));
      tau_x = KPARSEC*dust_dens*dust_cs;
      d_step = tau_r/tau_x;
      if (tau_x == 0) d_step = VERY_LARGE_NUMBER;
      if (d_step < step) {step = d_step; scatter = -1; dust_scatter = 1; }
      else  dust_scatter = 0;
      
      // take the step
      r[0] += D[0]*step;
      r[1] += D[1]*step;
      r[2] += D[2]*step;
    
      // see if we've exited (faster without sqrt)
      r_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];

      if (r_sq > r_outer*r_outer) {count_it = 1; break;}
      // see if we've gone under inner boundary
      // if (r_sq < r_emit*r_emit) {count_it = 0; break; }

      // if we line scattered, do it
      if (scatter >= 0)
      {
	flg_scatter = 1;
	xloc = (lam_loc/lambda_0[scatter] - 1)*C_LIGHT/v_doppler;

	// Get three velocity components of scatterer
 	u0 = voigt.Scatter_Velocity(xloc);
 	R10 =  gsl_rng_uniform(rangen);
 	R11 =  gsl_rng_uniform(rangen);

	// Original code
 	// u1 = sqrt(-1.0*log(R11))*cos(2*PI*R10);  
 	// u2 = sqrt(-1.0*log(R11))*sin(2*PI*R10); 

	// Calculate a*tau [take N_HI = integrated value from the current position]
	NHI = (r_outer - sqrt(r_sq)) * KPARSEC * n_0;
	atau = voigt_a * tau_const * NHI / Dnu;
	if (atau > 60)     xcr = 0.02*exp(1.4*pow(log(atau),0.6));
	else if (atau >1) xcr = 0.02*exp(0.6*pow(log(atau),1.2));
	else xcr = 0;
	//printf("#  xcr %.4e, NHI = %.4e \n", xcr, log10(NHI));
	u1  = sqrt(xcr*xcr-1.0*log(1.0-R11))*cos(2*PI*R10);
	u2  = sqrt(xcr*xcr-1.0*log(1.0-R11))*sin(2*PI*R10);
	
 	// parallel component
 	uvec[0] =  u0*D[0];
 	uvec[1] =  u0*D[1];
 	uvec[2] =  u0*D[2];
 	// perpindicular component 1
 	uvec[0] +=    u1*D[1];
 	uvec[1] += -1*u1*D[0];
 	// perpindicular component 2
 	uvec[0] += u2*D[0]*D[2];
 	uvec[1] += u2*D[1]*D[2];
 	uvec[2] += -1*u2*(D[0]*D[0] + D[1]*D[1]);
	
 	//  magnitude of Doppler shift into scatterer frame
 	vd_inc =  (uvec[0]*D[0] + uvec[1]*D[1] + uvec[2]*D[2]);
	
 	// choose new isotropic direction
 	mu  = 1 - 2.0*gsl_rng_uniform(rangen);
 	phi = 2.0*PI*gsl_rng_uniform(rangen);
 	sin_theta = sqrt(1 - mu*mu);
 	D[0] = sin_theta*cos(phi);
 	D[1] = sin_theta*sin(phi);
 	D[2] = mu;
	
 	// magnitude of shift out of scatterer frame
 	vd_out =(uvec[0]*D[0] + uvec[1]*D[1] + uvec[2]*D[2]);
	
 	// apply the doppler shift in and out
 	xloc = xloc - vd_inc + vd_out; 

	// go back to wavelength
	lam_loc = lambda_0[scatter]*(1 + xloc*v_doppler/C_LIGHT);

	// now get change in observer frame wavelength
	rad = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
	vel = Get_Velocity(r,rad)/C_LIGHT;
	if (rad > 0)
	  lam = lam_loc/(1 + vel*(r[0]*D[0] + r[1]*D[1] + r[2]*D[2])/rad);
      }

      if (dust_scatter)
      {
	double z =  gsl_rng_uniform(rangen);
	if (z > dust_albedo) {count_it = 0; break; }
	// choose new isotropic direction
 	mu  = 1 - 2.0*gsl_rng_uniform(rangen);
 	phi = 2.0*PI*gsl_rng_uniform(rangen);
 	sin_theta = sqrt(1 - mu*mu);
 	D[0] = sin_theta*cos(phi);
 	D[1] = sin_theta*sin(phi);
 	D[2] = mu;
      }
   }	
 
    // Count spectrum if needed
    double l_obs,t_obs,E_obs,m_obs,p_obs;
    if (count_it) 
    {
      l_obs = lam;
      t_obs = 0;
      E_obs = E_p;
      m_obs = D[2];
      p_obs = atan2(D[0],D[1]);
      if (p_obs < 0) p_obs += 2*PI;
      spectrum.Count(t_obs,l_obs,m_obs,p_obs,E_obs);
    }

    // Count un-absorbed photons
    if (flg_scatter == 0) 
      {
      l_obs = lam_emit;
      E_obs = E_p;
      m_obs = D[2];
      p_obs = atan2(D[0],D[1]);
      if (p_obs < 0) p_obs += 2*PI;
      spectrum.Scatter(t_obs,l_obs,m_obs,p_obs,E_obs);
    }
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



//------------------------------------------------------------
// Emit a single photon from some source
//------------------------------------------------------------
void Emit(double *r, double *D, double rphot)
{

  // either emit isotropically from a point
  if (rphot == 0)
  {
    // initial location
    r[0] = 0;
    r[1] = 0;
    r[2] = 0;
    
    // initial isotropic emission
    double mu  = 1 - 2.0*gsl_rng_uniform(rangen);
    double phi = 2.0*PI*gsl_rng_uniform(rangen);
    double sin_theta = sqrt(1 - mu*mu);
    D[0] = sin_theta*cos(phi);
    D[1] = sin_theta*sin(phi);
    D[2] = mu;
  }

  // or emit from the surface of a photosphere
  else
  {
    // pick initial position on photosphere
    double phi_core = 2*PI*drand48(); 
    double cosp_core  = cos(phi_core);
    double sinp_core  = sin(phi_core);
    double cost_core  = 1 - 2.0*drand48();
    double sint_core  = sqrt(1-cost_core*cost_core);
    // pick photon propogation direction wtr to local normal
    double phi_loc = 2*PI*drand48(); 
    // choose sqrt(R) to get outward, cos(theta) emission
    double cost_loc  = sqrt(drand48());
    double sint_loc  = sqrt(1 - cost_loc*cost_loc);
    // local direction vector
    double D_xl = sint_loc*cos(phi_loc);
    double D_yl = sint_loc*sin(phi_loc);
    double D_zl = cost_loc;
    // real coordinates
    r[0] = rphot*sint_core*cosp_core;
    r[1] = rphot*sint_core*sinp_core;
    r[2] = rphot*cost_core;
    // apply rotation matrix to convert vector into overall frame
    D[0] = cost_core*cosp_core*D_xl - sinp_core*D_yl + sint_core*cosp_core*D_zl;
    D[1] = cost_core*sinp_core*D_xl + cosp_core*D_yl + sint_core*sinp_core*D_zl;
    D[2] = -sint_core*D_xl +  cost_core*D_zl;
  }
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


