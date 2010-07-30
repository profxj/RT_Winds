#ifndef _MODEL_H
#define _MODEL 1
#include "cdf.hh"
#define LAM_LYA  1215.0

class MODEL
{

private:
  
  void Read_File(float *, const char*, char*);

public:

  int n_x;
  int n_pts;

  double opac_fac;
  double xmax, dx;
  float *dens, *temp, *x_HI;
  float *velx, *vely, *velz;
  float *opac, *vdop, *apar;
  float *J_UV, *popac;
  int   *indx, *indy, *indz;
  CDF emis;

  double l_start, l_stop;
  double L_tot;

  void   Init(const char *);
  void   Calculate_Opacities();
  int    Locate_Zone(double*);
  int    Emit(double *, double*);
  double P_esc(double, int*);

  double vdotD(int, double*);
  double v_doppler(int ind)    { return vdop[ind]; }
  double line_opacity(int ind) { return opac[ind]; }

  void Wipe_J_UV();
  void Reduce_J_UV();
  void Print_J_UV();
  void Compute_Photo_Opacity();
  void Compute_Ionization_State();
  double photo_opacity(int ind, double l);

  void MPI_Sum_Float_Array(float *, int);


};


#endif
