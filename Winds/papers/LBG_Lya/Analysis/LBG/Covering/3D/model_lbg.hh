#ifndef _MODEL_H
#define _MODEL 1
#include "cdf.hh"
#define LAM_LYA  1215.0

class MODEL
{

private:
  
  void Read_File(float *, const char*, char*);

public:

  float vdop, vinteract;
  CDF emis;

  double l_start, l_stop;
  double L_tot;

  // LBG wind parameters
  double r_inner;  // inner boundary radius, in kpc [Should always be 1 !!]
  double r_outer;  // outer boundary radius, in kpc
  double r_emit;   // boundary to emit from
  
  double LBG_gamma;          // Covering fraction parameter
  double LBG_fc;             // Maximum covering fraction (for MgII 2796)
  double LBG_vmax;           // cm/s
  double LBG_alpha;          // Velocity field parameter
  double LBG_A;              // Wind parameter
  double LBG_Aa;


  void   Init(const char *);
  void   Calculate_Opacities();
  int    Locate_Zone(double*);
  int    Emit(double *, double*);
  double P_esc(double, int*);

  double Velocity(double);
  double DvDr(double);
  double Covering(double);
  double v_doppler()    { return vdop; }
  double v_interact()    { return vinteract; }

  void Wipe_J_UV();
  void Reduce_J_UV();
  void Print_J_UV();
  void Compute_Photo_Opacity();
  void Compute_Ionization_State();
  double photo_opacity(int ind, double l);

  void MPI_Sum_Float_Array(float *, int);


};


#endif
