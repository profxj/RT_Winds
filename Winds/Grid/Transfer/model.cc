#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "model.hh"
#include "math.h"
#include "physical_constants.hh"

void MODEL::Init(const char* fname)
{
  FILE *in = fopen(fname,"r");
  if (in == NULL) {
    printf("Can't open model file %s\n",fname);
    exit(0); }
  
  int n;
  double xx;
  
  fscanf(in,"%d %lf\n",&n,&xx);
  n_x = n;
  n_pts = n_x*n_x*n_x;
  dx   = xx;
  xmax = dx*n_x;
  
  dens = new float[n_pts];
  temp = new float[n_pts];
  x_HI = new float[n_pts];
  velx = new float[n_pts];
  vely = new float[n_pts];
  velz = new float[n_pts];
  opac = new float[n_pts];
  apar = new float[n_pts];
  vdop = new float[n_pts];
  J_UV = new float[n_pts];
  popac = new float[n_pts];
  emis.New(n_pts); 
  indx = new int[n_pts];
  indy = new int[n_pts];
  indz = new int[n_pts];


  L_tot = 0;
  double cell_vol = pow(dx*KILOPARSEC,3);

  double x1,x2,x3,x4,x5,x6;
  for (int i=0;i<n_pts;i++)
  {
    fscanf(in,"%lf %lf %lf %lf %lf %lf\n",&x1,&x2,&x3,&x4,&x5,&x6);
    dens[i] = x1;
    velx[i] = x2*1e5;
    vely[i] = x3*1e5;
    velz[i] = x4*1e5;
    emis.Put(i,x5*cell_vol);
    vdop[i] = x6*1e5;
    if (x6 == 0) vdop[i] = 1e5;
    temp[i] = 0.5*M_PROTON*vdop[i]*vdop[i]/K_BOLTZ;
    apar[i] = 4.7e-3*sqrt(1e4/temp[i]);

    L_tot += emis.Value(i);
  }
  emis.Normalize();
  L_tot=1;
}


int MODEL::Locate_Zone(double *xr)
{
  int x = (int)((xr[0])/dx);
  int y = (int)((xr[1])/dx);
  int z = (int)((xr[2])/dx);
  if ((x < 0)||(x >= n_x)) return -1;
  if ((y < 0)||(y >= n_x)) return -1;
  if ((z < 0)||(z >= n_x)) return -1;
  int c = x*n_x*n_x + y*n_x + z;
  
  return c;

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

void MODEL::Compute_Ionization_State()
{
  for (int i=0;i<n_pts;i++)
  {
  }
}


void MODEL::Print_J_UV()
{
  int i = n_x/2;
  int j = n_x/2; //for (int j=0;j<n_x;j++)
  for (int k=0;k<n_x;k++)
  {
    int ind = i*n_x*n_x + j*n_x + k;
    printf("%e %e\n",(k - n_x/2)*dx,J_UV[ind]);
  }
}

void MODEL::Reduce_J_UV()
{
  // Sum up over all procs
  MPI_Sum_Float_Array(J_UV,n_pts);

  double cvol = pow(dx*KILOPARSEC,3);
  cvol = cvol/KILOPARSEC; // length scale
  
  // get # of processors
  int size;
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  // divide by n_procs
  for (int i=0;i<n_pts;i++) J_UV[i] = J_UV[i]/size/cvol;
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
