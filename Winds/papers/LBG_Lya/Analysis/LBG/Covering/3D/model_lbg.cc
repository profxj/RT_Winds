#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "model_lbg.hh"
#include "math.h"
#include "physical_constants.hh"
#include "Lua.hh"

void MODEL::Init(const char* fname)
{
  Lua lua;
  lua.init( fname );

  // Dimensions
  r_inner = lua.scalar<double>("r_inner");
  r_outer = lua.scalar<double>("r_outer");
  r_emit = lua.scalar<double>("r_emit");

  // Wind
  LBG_gamma = lua.scalar<double>("LBG_gamma");
  LBG_fc = lua.scalar<double>("LBG_fc");
  LBG_vmax = lua.scalar<double>("LBG_vmax");
  LBG_alpha = lua.scalar<double>("LBG_alpha");
  LBG_Reff = lua.scalar<double>("LBG_Reff");
  vdop = lua.scalar<double>("vdoppler");
  vinteract = lua.scalar<double>("vinteract");

  // Calcualte a few
  // LBG_A =  LBG_vmax*LBG_vmax * (1.-LBG_alpha);  // Wind parameter
  LBG_A =  LBG_vmax*LBG_vmax * (1.-LBG_alpha) / 
    ( pow(1.,1-LBG_alpha) - pow(LBG_Reff, 1-LBG_alpha));  // Wind parameter
  LBG_Aa =  sqrt(LBG_A/(1-LBG_alpha));

  // Normalize
  L_tot = 1;
}



double MODEL::Velocity(double r)
{
  if (r <= r_inner) return 0;
  return LBG_Aa * sqrt(1 - pow(r,1-LBG_alpha));  // Assumes r_inner=1kpc
}

double MODEL::DvDr(double r)  // cm/s per kpc
{
  if (r <= r_inner) return 0;
  double dvdr = LBG_Aa * 0.5 / sqrt(1-pow(r,1-LBG_alpha)) * (LBG_alpha-1) * pow(r,-1*LBG_alpha);
  return dvdr;
}

double MODEL::Covering(double r)
{
  if (r <= r_inner) return 0.;
  return LBG_fc * pow(r/r_inner, -1*LBG_gamma);
}



//------------------------------------------------------------
// General MPI helper functions
//------------------------------------------------------------
void MODEL::MPI_Sum_Float_Array(float *array, int n_el)
{
  float *new_ptr;
  int i;          
  
  // allocate the memory for new pointer
  new_ptr = new float[n_el];
  if (new_ptr == NULL) printf("YO! MPI Sum problem\n");
  // zero out array
  for (i=0;i<n_el;i++) new_ptr[i] = 0;
  MPI_Allreduce(array,new_ptr,n_el,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  // put back into place
  for (i=0;i<n_el;i++) array[i] = new_ptr[i];
  // free up the memory
  delete new_ptr;
}
