#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "physical_constants.hh"
#include "spectrum.hh"

//***************************************************************
// Constructors
//***************************************************************

SPECTRUM::SPECTRUM()
{
  strcpy(name,DEFAULT_NAME);
  n_elements = 0;
  n_escaped = 0;
}


void SPECTRUM::Init(double *t, double *l, int nmu, int nphi)
{
  
  // initalize grids
  lambda_grid.New(l[0],l[1],l[2]);
  time_grid.New(t[0],t[1],t[2]);
  mu_grid.Uniform_New(-1,2.0/nmu,nmu);
  phi_grid.Uniform_New(0,2*PI/nphi,nphi);
  
  n_times  = time_grid.Size();
  n_lambda = lambda_grid.Size();
  n_mu     = mu_grid.Size();
  n_phi    = phi_grid.Size();
  
  a3 = phi_grid.Size();
  a2 = mu_grid.Size()*a3;            ;
  a1 = lambda_grid.Size()*a2;
  n_elements = time_grid.Size()*a1;
  
  // allocate
  click  = new double[n_elements];
  count  = new double[n_elements];

  // clear 
  Wipe();
}

void SPECTRUM::Init(double *l, int nmu, int nphi)
{
  double t[] = {1,1,1};
  Init(t,l,nmu,nphi);
}


void SPECTRUM::Set_Name(char *n)
{
  strcpy(name,n);
}

//***************************************************************
// Functional procedure: Wipe
//***************************************************************
void SPECTRUM::Wipe()
{
  n_escaped = 0;
  for (int i=0;i<n_elements;i++) {
    click[i]  = 0;
    count[i]  = 0; }
}


//***************************************************************
// handles the indexing: should be called in this order
//    time, wavelength, mu, phi
//***************************************************************
int SPECTRUM::Index(int t, int l, int m, int p)
{
  return t*a1 + l*a2 + m*a3 + p;
}

//***************************************************************
//--------------------------------------------------------------
//********************** COUNT FUNCTIONS ************************
//--------------------------------------------------------------
//***************************************************************

void SPECTRUM::Count(double t, double l, double mu, double phi, double E)
{
  n_escaped += E;;

  // locate bin number in all dimensions
  int t_bin = time_grid.Locate(t);
  int l_bin = lambda_grid.Locate(l);
  int m_bin = mu_grid.Locate(mu);
  int p_bin = phi_grid.Locate(phi);

  // if off the grids, just return without counting
  if ((t_bin < 0)||(l_bin < 0)||(m_bin < 0)||(p_bin < 0)) return;
  
  // add to counters
  int index      = Index(t_bin,l_bin,m_bin,p_bin);
  count[index]  += E;
  click[index]  += 1;
}


//***************************************************************
//--------------------------------------------------------------
//********************** PRINT FUNCTIONS *************************
//--------------------------------------------------------------
//***************************************************************

void SPECTRUM::Print_Spectrum()
{
  char basename[10000],filename[10000];
  sprintf(filename,"%s",name);
  FILE *out;
  out = fopen(filename,"w");

  int id;

  for (int i=0;i<n_times;i++) 
    for (int m=0;m<n_mu;m++) 
      for (int p=0;p<n_phi;p++)
	for (int j=0;j<n_lambda;j++) 
	{
	  id = Index(i,j,m,p);
	  if (n_times > 1) fprintf(out,"%15.5e ",time_grid.Center(i));
	  if (n_mu > 1) fprintf(out,"%15.5e ",mu_grid.Center(m));
	  if (n_phi > 1) fprintf(out,"%15.5e ",phi_grid.Center(p));
	  fprintf(out,"%15.5e %15.5e %15.5e\n",
	    lambda_grid.Center(j),count[id],click[id]);
	}
  
  fclose(out);
}



void SPECTRUM::Print_Lightcurve()
{
  int i,j;
  double sum;
  
  char filename[10000];
  sprintf(filename,"%s.lc",name);
  FILE *out = fopen(filename,"w");
  for (i=0;i<n_times;i++)
  {  
    sum = 0;
    for (j=0;j<n_lambda;j++) sum += count[Index(i,j,0,0)];
    fprintf(out,"%12.5e %12.5e\n",time_grid.Center(i),sum);
  }
  fclose(out);

}

//***************************************************************
//--------------------------------------------------------------
//******************* MAGNITUDE FUNCTIONS ***********************
//--------------------------------------------------------------
//***************************************************************



//***************************************************************
//--------------------------------------------------------------
//********************** MPI FUNCTIONS **************************
//--------------------------------------------------------------
//***************************************************************

void SPECTRUM::Normalize()
{
  double norm_factor =  1.0/(1.0*mu_grid.Size()*phi_grid.Size());
  
  // renormalize flux
   for (int i=0;i<time_grid.Size();i++) 
    for (int j=0;j<lambda_grid.Size();j++) 
      for (int m=0;m<mu_grid.Size();m++) 
	for (int p=0;p<phi_grid.Size();p++)
	  count[Index(i,j,m,p)] /= 
	    time_grid.Delta(i)*lambda_grid.Delta(j)*norm_factor;
}


void SPECTRUM::Rescale(double r)
{
  for (int i=0;i<n_elements;i++) count[i] *= r;
}


void SPECTRUM::MPI_Sum_All()
{
  MPI_Allreduce_Array(count);
  MPI_Allreduce_Array(click);
}

void SPECTRUM::MPI_Average_All()
{
  MPI_Sum_All();

  int mpi_procs;
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );
  for (int i=0;i<n_elements;i++) {
    count[i] /= mpi_procs; }
}


void SPECTRUM::MPI_Allreduce_Array(double *arr)
{
  double *new_ptr,*this_ptr;
  int j;          

  // do it in a1 chunks
  int chunk = a1;
  // allocate the memory for new pointer
  new_ptr = new double[chunk];
         
  for (int i=0;i<time_grid.Size();i++)
  {
    // zero out array
    for (j=0;j<chunk;j++) new_ptr[j] = 0;
    // reduce the stuff
    this_ptr = &(arr[i*chunk]);
    MPI_Allreduce(this_ptr,new_ptr,chunk,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    // put back into place
    for (j=0;j<chunk;j++) arr[i*chunk + j] = new_ptr[j];
  }
  // free up the memory
  delete new_ptr;
}
