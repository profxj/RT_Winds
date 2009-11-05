#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "locate_array.hh"

using namespace std;

//---------------------------------------------------------
// Just allocation the memory for this
//---------------------------------------------------------
void LOCATE_ARRAY::New(int n) 
{
  if (n<=0) n_elements = 1;
  else n_elements = n;
  low    = new double[n_elements];
  high   = new double[n_elements];
}

//---------------------------------------------------------
// Initialize read from file
//---------------------------------------------------------
void LOCATE_ARRAY::New(const char* filen)
{
  // Read in data from file
  double tmp;
  vector<double> bounds;
  FILE *in = fopen(filen,"r");
  if (in == NULL) {
    printf("# ERROR: Can't open array file %s\n",filen);
    return; }

  while(!feof(in)) {
    fscanf(in,"%lf\n",&tmp);
    bounds.push_back(tmp); }
  fclose(in);
  
  New(bounds.size()-1);
  
  for (int i=0;i<n_elements;i++)
  {
    low[i]    = bounds[i];
    high[i]   = bounds[i+1];
    if (high[i]<=low[i]) printf("# Warning, Locate Array is non-monotonic\n");
  }
}

//---------------------------------------------------------
// Initialize with start, stop and delta
//---------------------------------------------------------
void LOCATE_ARRAY::New(double start, double stop, double delta)
{
  int n = (int)((stop-start)/delta);
  if (n <= 1) {n = 1; delta = stop-start;}
  if (delta < 0) delta = 0;
  New(n);
  for (int i=0;i<n_elements;i++) {
    low[i]    = start + i*delta;
    high[i]   = start + (i+1)*delta;}
}

//---------------------------------------------------------
// Initialize with uniform spacing
//---------------------------------------------------------
void LOCATE_ARRAY::Uniform_New(double start, double delta, int n)
{
  New(n);
  if (n == 1) delta = 0;
  for (int i=0;i<n_elements;i++) {
    low[i]    = start + i*delta;
    high[i]   = start + (i+1)*delta;}
}

//---------------------------------------------------------
// Initialize with logarithmic spacing
//---------------------------------------------------------
void LOCATE_ARRAY::Log_New(double start, double delta, int n)
{
  New(n);
  if (n == 1) delta = 0;
  for (int i=0;i<n_elements;i++) {
    low[i]  = start*pow(10.0,i*delta);
    high[i] = start*pow(10.0,(i+1)*delta);    }
}


void LOCATE_ARRAY::New(int n, double *bounds)
{
  New(n);
  
  if (n == 0) {
    low[0]  = bounds[0];
    high[0] = bounds[0]; }
  
  else for (int i=0;i<n_elements;i++)
  {
    low[i]    = bounds[i];
    high[i]   = bounds[i+1];
    if (high[i] <= low[i]) printf("Warning, Locate Array is non-monotonic\n");
  }
}


// Copy
void LOCATE_ARRAY::Copy(LOCATE_ARRAY la)
{
  New(la.n_elements);
  for (int i=0;i<n_elements;i++) {
    low[i]  = la.low[i];
    high[i] = la.high[i];}
}

void LOCATE_ARRAY::Renormalize(double d)
{
  for (int i=0;i<n_elements;i++) {
    low[i]  = low[i]*d;
    high[i] = high[i]*d;}
}


//---------------------------------------------------------
// LOCATE 
// If off the boundaries of the array, return negative
// numbers as flag
//---------------------------------------------------------
int LOCATE_ARRAY::Locate(double x)
{
  // the degenerate case always returns 0
  if (n_elements == 1) return 0;
  
  // a form of locate from numerical recipes
  int bm;                             // mid point
  int bl = 0;                         // lower bound
  int bu = n_elements-1;              // upper bound
  
  // check if we are off the ends of the array
  if (x >= high[bu]) return -2;
  if (x <  low[bl])  return -1;
  // check for first bin
  if (x < high[bl])  return 0;
  
  // search the array for this index
  while (bu-bl > 1)
  {
    bm = (bu+bl)/2;
    if (high[bm] <= x) bl = bm;
    else bu = bm;
  }
  return bu;  
} 

//---------------------------------------------------------
// LOCATE_Bounded
// If off the boundaries of the array, return the
// boundary value
//---------------------------------------------------------
int LOCATE_ARRAY::Locate_Bounded(double x)
{
  // the degenerate case always returns 0
  if (n_elements == 1) return 0;
  
  // a form of locate from numerical recipes
  int bm;                             // mid point
  int bl = 0;                         // lower bound
  int bu = n_elements-1;              // upper bound

  // check if we are off the ends of the array
  if (x >= high[bu]) return n_elements-1;
  if (x < high[bl])  return 0;
    
  // search the array for this index
  while (bu-bl > 1)
  {
    bm = (bu+bl)/2;
    if (high[bm] <= x) bl = bm;
    else bu = bm;
  }
  return bu;  
} 

//---------------------------------------------------------
// Linear Interpolation of a passed array
//---------------------------------------------------------
double LOCATE_ARRAY::Interpolate(double *y, double x)
{
  int i = Locate_Bounded(x);
  double slope;
 
  if (i == n_elements-1)  slope = (y[i]-y[i-1])/(Center(i) - Center(i-1));
  else if (i == 0)        slope = (y[i+1]-y[i])/(Center(i+1) - Center(i));
  else if (x < Center(i)) slope = (y[i]-y[i-1])/(Center(i) - Center(i-1));
  else                    slope = (y[i+1]-y[i])/(Center(i+1) - Center(i));

  double v = y[i] + slope*(x - Center(i));
  return v;
}



double LOCATE_ARRAY::Sample(int i)
{
  return Low(i) + drand48()*Delta(i);
}


void LOCATE_ARRAY::Print()
{
  printf("# Print Locate Array; n_elements = %d\n",n_elements);
  for (int i=0;i<n_elements;i++) Print(i);
}
  
void LOCATE_ARRAY::Print(int i)
{
  if ((i<0)||(i>=n_elements)) printf("Print: Index out of Array bounds\n");
  else printf("%4d %12.4f %12.4f %12.4f %12.4f\n",
	      i,Center(i),Low(i),High(i),Delta(i));
}
  
