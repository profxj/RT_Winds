#include <time.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf.h>
#include "emit_line.hh"


void EMIT_LINE::New(int n, double EW, double sigma, double wave_start, double wave_end)
{
  n_x = n;
  EW_x = EW;
  sigma_x = sigma*1215.6701/2.99792e10;  // Ang

  // Create wavelength array
  wave = new double[n_x];
  double dwv = (wave_end-wave_start)/float(n_x-1);
  for (int i=0;i<n_x;i++) {
    wave[i]    = wave_start + i*dwv;
  }

  // Fill in
  Compute_Cumul();
}

void EMIT_LINE::Compute_Cumul()
{
  int i, j;
  double uval;

  double u2 = (wave[n_x-1] - 1215.6701) / sqrt(2) / sigma_x;
  double u1 = (wave[0] - 1215.6701) / sqrt(2) / sigma_x;
  double norm = EW_x * ( gsl_sf_erf(u2) - gsl_sf_erf(u1)) + (wave[n_x-1]-wave[0]);

  double intg[10000];

  // Integrate
  for (i=0; i<n_x; i++) {
    uval = (wave[i] - 1215.6701) / sqrt(2) / sigma_x;
    intg[i] = EW_x * ( gsl_sf_erf(uval) - gsl_sf_erf(u1)) + (wave[i]-wave[0]);
    intg[i] = intg[i]/norm;
    // printf("# wave %.7e,  Cumul %.3e \n", wave[i], intg[i]);
  }


  // Initialize locate array
  cumul.New(n_x, intg); 


}


double EMIT_LINE::Get_Wave(double y)
{
  int i = cumul.Locate(y);
  if (i < 0) return 0;
  if (i >= cumul.Size()-1) return wave[i];

  // linear interpolation
  double dy = (y - cumul.Center(i))/(cumul.Center(i+1) - cumul.Center(i));
  double val = wave[i] + dy*(wave[i+1] - wave[i]);
  printf("# y %.7e,  i %d, wave %.6e \n", y, i, val);
  return val;
}

