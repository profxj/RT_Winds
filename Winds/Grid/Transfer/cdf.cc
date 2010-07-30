#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "cdf.hh"

CDF::~CDF()
{
 delete y;
}

//--------------------------------------------------------
// Initialize by allocating memory, making y uniform
//--------------------------------------------------------
void CDF::New(int n) 
{
  n_elements = n; 
  y = new double[n];
  for (int i=0;i<n;i++) y[i] = 1; 
}


//------------------------------------------------------
// return the actual y value, not the integrated
//------------------------------------------------------
double CDF::Value(int i)   {
  if (i==0) return y[0];
  else return (y[i] - y[i-1]);  }

//------------------------------------------------------
// Set Y from an entire given array
//------------------------------------------------------
void CDF::Set(double *f)
{
  double sum = 0;
  for (int i=0;i<n_elements;i++) {
    sum  += f[i];
    y[i]  = sum; }
  Normalize();
}


//------------------------------------------------------
// Set Y from an entire given array and weights
//------------------------------------------------------
void CDF::Set(double *f, double *w)
{
  double sum = 0;
  for (int i=0;i<n_elements;i++) {
    sum  += f[i]*w[i];
    y[i]  = sum; }
  Normalize();
}
  
//------------------------------------------------------
// Set a uniform distribution up to cell c
//------------------------------------------------------
void CDF::Uniform_UpTo(int c)
{
  int i;
  for (i=0;i<c;i++) y[i] = i+1;
  for (i=c;i<n_elements;i++) y[i] = c;
  Normalize();
}

//------------------------------------------------------
// Normalize such that the last entry is 1.0
//------------------------------------------------------
void CDF::Normalize() 
{
  int i;
  if (y[n_elements-1] == 0)
    for (i=1;i<n_elements;i++) y[i] = y[i-1] + 1.0;  
  for (i=0;i<n_elements;i++) y[i] /= y[n_elements-1]; 
}


//---------------------------------------------------------
// Sample the probability distribution using binary search
// A form of locate from numerical recipes.  Return
// the bin number
//---------------------------------------------------------
int CDF::Sample()
{
  // random choose a value
  double z = drand48();
  //printf("%f ",z);

  // mid, lower, and upper points
  int bm;                         // mid point
  int bl = 0;                     // lower bound
  int bu = n_elements - 1;        // upper bound
  // see if we are off the top or bottom
  if (z > y[bu]) return bu;
  if (z < y[bl]) return bl;
  // search
  while (bu-bl > 1)
  {
    bm = (bu+bl)/2;
    if (y[bm] <= z) bl = bm;
    else bu = bm;
  }
  return bu;
}


//------------------------------------------------------
// Simple printout
//------------------------------------------------------
void CDF::Print() {
  for (int i=0;i<n_elements;i++) 
    printf("%5d %10.4e %10.4e\n",i,Value(i),y[i]);
}
  
//------------------------------------------------------
// Clear the arrays
//------------------------------------------------------
void CDF::Wipe()
{
  int i;
  for (i=0;i<n_elements;i++)  y[i] = 0;
}
  
//------------------------------------------------------
// MPI Reduce this class
//------------------------------------------------------
void CDF::MPI_Combine() {
  MPI_Sum_Array(y);
}

//--------------------------------------------------------
// MPI Reduce helper
//--------------------------------------------------------
void CDF::MPI_Sum_Array(double *array)
{
    double *new_ptr;
    int i;          
    
    // allocate the memory for new pointer
    new_ptr = new double[n_elements];
    // zero out array
    for (i=0;i<n_elements;i++) new_ptr[i] = 0;
    MPI_Allreduce(array,new_ptr,n_elements,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    // put back into place
    for (i=0;i<n_elements;i++) array[i] = new_ptr[i];
    // free up the memory
    delete new_ptr;
}   

//****************************************************************
// CDF_FLOAT functions (I know this is bad programming to repeat)
//****************************************************************
CDF_FLOAT::~CDF_FLOAT()
{
 delete y;
}

//--------------------------------------------------------
// Initialize by allocating memory, making y uniform
//--------------------------------------------------------
void CDF_FLOAT::New(int n) 
{
  n_elements = n; 
  y = new float[n];
  for (int i=0;i<n;i++) y[i] = 1; 
}


//------------------------------------------------------
// return the actual y value, not the integrated
//------------------------------------------------------
float CDF_FLOAT::Value(int i)   {
  if (i==0) return y[0];
  else return (y[i] - y[i-1]);  }

//------------------------------------------------------
// Set Y from an entire given array
//------------------------------------------------------
void CDF_FLOAT::Set(double *f)
{
  float sum = 0;
  for (int i=0;i<n_elements;i++) {
    sum  += f[i];
    y[i]  = sum; }
  Normalize();
}


//------------------------------------------------------
// Set Y from an entire given array and weights
//------------------------------------------------------
void CDF_FLOAT::Set(double *f, double *w)
{
  float sum = 0;
  for (int i=0;i<n_elements;i++) {
    sum  += f[i]*w[i];
    y[i]  = sum; }
  Normalize();
}
  
//------------------------------------------------------
// Set a uniform distribution up to cell c
//------------------------------------------------------
void CDF_FLOAT::Uniform_UpTo(int c)
{
  int i;
  for (i=0;i<c;i++) y[i] = i+1;
  for (i=c;i<n_elements;i++) y[i] = c;
  Normalize();
}

//------------------------------------------------------
// Normalize such that the last entry is 1.0
//------------------------------------------------------
void CDF_FLOAT::Normalize() 
{
  int i;
  if (y[n_elements-1] == 0)
    for (i=1;i<n_elements;i++) y[i] = y[i-1] + 1.0;  
  for (i=0;i<n_elements;i++) y[i] /= y[n_elements-1]; 
}


//---------------------------------------------------------
// Sample the probability distribution using binary search
// A form of locate from numerical recipes.  Return
// the bin number
//---------------------------------------------------------
int CDF_FLOAT::Sample()
{
  // random choose a value
  double z = drand48();
  //printf("%f ",z);

  // mid, lower, and upper points
  int bm;                         // mid point
  int bl = 0;                     // lower bound
  int bu = n_elements - 1;        // upper bound
  // see if we are off the top or bottom
  if (z > y[bu]) return bu;
  if (z < y[bl]) return bl;
  // search
  while (bu-bl > 1)
  {
    bm = (bu+bl)/2;
    if (y[bm] <= z) bl = bm;
    else bu = bm;
  }
  return bu;
}


//------------------------------------------------------
// Simple printout
//------------------------------------------------------
void CDF_FLOAT::Print() {
  for (int i=0;i<n_elements;i++) 
    printf("%5d %10.4e %10.4e\n",i,Value(i),y[i]);
}
  
//------------------------------------------------------
// Clear the arrays
//------------------------------------------------------
void CDF_FLOAT::Wipe()
{
  int i;
  for (i=0;i<n_elements;i++)  y[i] = 0;
}
  
//------------------------------------------------------
// MPI Reduce this class
//------------------------------------------------------
void CDF_FLOAT::MPI_Combine() {
  MPI_Sum_Float_Array(y);
}

//--------------------------------------------------------
// MPI Reduce helper
//--------------------------------------------------------
void CDF_FLOAT::MPI_Sum_Float_Array(float *array)
{
  float *new_ptr;
  int i;          
    
  // allocate the memory for new pointer
  new_ptr = new float[n_elements];
  // zero out array
  for (i=0;i<n_elements;i++) new_ptr[i] = 0;
  MPI_Allreduce(array,new_ptr,n_elements,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  // put back into place
  for (i=0;i<n_elements;i++) array[i] = new_ptr[i];
  // free up the memory
  delete new_ptr;
}   


