#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "Lua.hh"
#include "locate_array.hh"
#include "spectrum.hh"
#include "model_lbg.hh"
#include "lines.hh"


// globals
gsl_rng  *rangen;    // random number generator
double   n_photons;
int      verbose;     
MODEL    model;
SPECTRUM spectrum;
double   D_obs[3], cos_obs[4];
LINES    lines;


//--------------------------------------------
// the main program
//--------------------------------------------
int main(int argc, char **argv)
{

  void Run_Monte_Carlo(double);
  void Run_UV_Transport(double, double, int);

  // initialize MPI for parallelism
  int my_rank, n_procs;
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs );
  if (my_rank == 0) verbose = 1; else verbose = 0;
  if (verbose) printf("# Using %d processors\n",n_procs);

  // setup and seed random number generator
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  gsl_rng_default_seed = (unsigned int)time(NULL) + my_rank;
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);
  
   // open up the lua parameter file
  Lua lua;
  std::string script_file = "param.lua";  // This default is *not* functional
  if( argc > 1 ) script_file = std::string( argv[ 1 ] ); else return 0;
  lua.init( script_file );

  // Photons and viewing angle
  n_photons =  lua.scalar<double>("n_photons");
  double O_theta = lua.scalar<double>("O_theta"); 
  double O_phi   = lua.scalar<double>("O_phi");
  D_obs[0]  = sin(O_theta/180.0*PI)*cos(O_phi/180.0*PI);
  D_obs[1]  = sin(O_theta/180.0*PI)*sin(O_phi/180.0*PI);
  D_obs[2]  = cos(O_theta/180.0*PI);
  cos_obs[0] = cos(O_theta/180.0*PI);
  cos_obs[1] = sin(O_theta/180.0*PI);
  cos_obs[2] = cos(O_phi/180.0*PI);
  cos_obs[3] = sin(O_phi/180.0*PI);

  // read lines
  std::string linefile = lua.scalar<std::string>("line_file");
  lines.Init(linefile.c_str());

  // initalize the model
  std::string infile = lua.scalar<std::string>("model_file");
  model.Init(infile.c_str());
  model.l_start = lua.scalar<double>("l_start");
  model.l_stop  = lua.scalar<double>("l_stop");
  //model.opac_fac  = lua.scalar<double>("opac_factor");
  if (verbose) printf("# Model Read\n");
  //model.Calculate_Opacities();
  if (verbose) printf("# Model set\n");
  if (verbose) printf("# L_tot = %e\n",model.L_tot);


  // define spectrum counter
  double lgrid[3];
  lgrid[0] = lua.scalar<double>("l_start");
  lgrid[1] = lua.scalar<double>("l_stop");
  lgrid[2] = lua.scalar<double>("l_delta");
  int l_nx = lua.scalar<int>("l_nx");
  double l_xmax = lua.scalar<double>("l_xmax");
  spectrum.Init(lgrid,l_nx,l_xmax);
  spectrum.Set_Name("spec.dat");

  // Run the photoionization calculation
  //double n_uv =  lua.scalar<double>("n_uv");
  //double L_uv =  lua.scalar<double>("L_uv");
  //int iter_uv =  lua.scalar<double>("iter_uv");
  //n_uv = n_uv/n_procs;
  //if (verbose) 
  //{
  //  if (n_uv > 0)
  //    printf("# Sending %e UV phots on %d procs (%e total photons)\n",
	//     n_uv,n_procs,n_uv*n_procs);
  //  else printf("# No UV photons sent\n");
  //}
  //if (n_uv > 0) Run_UV_Transport(n_uv, L_uv, iter_uv);


  // get photons per processor
  n_photons = n_photons/n_procs;
  if (verbose) 
  {
    if (n_photons > 0)
      printf("# Sending %e photons on %d procs (%e total photons)\n",
	     n_photons,n_procs,n_photons*n_procs);
    else printf("# No line photons sent\n");
  }
		
  // Do the monte carlo calculation
  if (n_photons > 0) Run_Monte_Carlo(n_photons);

  // finish up mpi
  MPI_Finalize();

}
