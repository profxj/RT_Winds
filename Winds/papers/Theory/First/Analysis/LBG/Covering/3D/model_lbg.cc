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
  model.r_inner = lua.scalar<double>("r_inner");
  model.r_outer = lua.scalar<double>("r_outer");
  model.r_emit = lua.scalar<double>("r_emit");

  // Wind
  model.LBG_gamma = lua.scalar<double>("LBG_gamma");
  model.LBG_fc = lua.scalar<double>("LBG_fc");
  model.LBG_vmax = lua.scalar<double>("LBG_vmax");
  model.LBG_alpha = lua.scalar<double>("LBG_alpha");
  model.vdop = lua.scalar<double>("vdoppler");

  // Calcualte a few
  model.LBG_A =  model.LBG_vmax*model.LBG_vmax * (1.-model.LBG_alpha);  // Wind parameter
  model.LBG_Aa =  sqrt(model.LBG_A/(1-model.LBG_alpha));
}



void MODEL::Calculate_Opacities(void)
{
  double nu_Lya   = C_LIGHT/(LAM_LYA*ANGS_TO_CM);
  double flu_Lya  = 0.42;
  double Lyman_CS = CLASSICAL_CS*flu_Lya/nu_Lya;
  double E_Lya    = 1.63422e-11;

  // For your reference
  // sigma(Lya) = CLASSICAL_CS*flu/nu_dop
  // or sigma(Lya) = 1e-13 cm^2 (T/10^4)^(-1/2)
  
  
  double N_e, N_H, Clya, this_e;
  int ind;
  for (int i=0;i<n_x;i++)
    for (int j=0;j<n_x;j++)
      for (int k=0;k<n_x;k++)
      {
	ind = i*n_x*n_x + j*n_x + k;
	indx[ind] = i;
	indy[ind] = j;
	indz[ind] = k;

	Clya = 2.41e-6*pow(temp[ind]/1e4,0.22)
	  *exp(-10.2/K_BOLTZ_EV/temp[ind])*pow(temp[ind],-0.5);
	
	N_H = dens[ind];
	opac[ind] = N_H*(C_LIGHT/vdop[ind]);
	opac[ind] *= opac_fac;

      }

}




int MODEL::Emit(double *r, double *l)
{
  int ind = emis.Sample();

  r[0] = dx*indx[ind] + dx*drand48();
  r[1] = dx*indy[ind] + dx*drand48();
  r[2] = dx*indz[ind] + dx*drand48();

  double off = vdop[ind]/C_LIGHT;
  double z = drand48();
  if (z > 0.5) off *= -1;

  *l = l_start + (l_stop - l_start)*drand48();

  return ind;
}


double MODEL::vdotD(int ind, double *D)
{
  double vd = D[0]*velx[ind] + D[1]*vely[ind] + D[2]*velz[ind];
  return vd;
}



void MODEL::Wipe_J_UV()
{
  for (int i=0;i<n_pts;i++) J_UV[i] = 0;
}


void MODEL::Compute_Photo_Opacity()
{
  double H_PI_CS = 6*1e-18;

  for (int i=0;i<n_pts;i++) 
    popac[i] = dens[i]*x_HI[i]*H_PI_CS;

}


double MODEL::photo_opacity(int ind, double l) 
{ 
  return popac[ind]*pow(l/912,-3.0); 
  
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
