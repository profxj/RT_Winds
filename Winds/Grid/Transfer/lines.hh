#ifndef _LINES_H
#define _LINES_H
#include "physical_constants.hh"

struct SINGLE_LINE
{
  double lam, nu, fX, f_lu, gl,gu, X, A_ul;
  int atom, ion;
  double P_scat;
  int    n_branch;
  double *blam;
  double *bprob;

};

class LINES
{
  int n_lines;
  SINGLE_LINE *l;

public:

  int n() { return n_lines; }
  double lambda(int i)   {return  l[i].lam; }
  int    n_branch(int i) {return  l[i].n_branch; }
  double P_scat(int i)   {return  l[i].P_scat;}
  double cs(int i)       { return l[i].fX; }
  double blam(int i, int j)  { return l[i].blam[j]; }
  double bprob(int i, int j) { return l[i].bprob[j]; }

  void Init(const char *filen)
  {

    double solar_abun[]  = {0.0,
        1.000000000,    0.0977237,  1.44544e-11,  1.41254e-11, 3.98107e-10,
        0.000354814,  9.33254e-05,  0.000741310,  3.63078e-08, 0.000123027,
        2.13796e-06,  3.80189e-05,  2.95121e-06,  3.54814e-05, 2.81838e-07,
        1.62181e-05,  3.16228e-07,  3.63078e-06,  1.31826e-07, 2.29087e-06,
        1.47911e-09,  1.04713e-07,  1.00000e-08,  4.67735e-07, 2.45471e-07,
        3.16228e-05,  8.31764e-08,  3.31131e-06,  1.62181e-08, 3.98107e-08,
        7.58578e-10,  1.69824e-09,  3.98107e-10,  7.94328e-10, 1.73780e-10,
        3.98107e-10,  2.63027e-11,  1.31826e-11,  1.34896e-10, 1.65959e-11};


    FILE *in = fopen(filen,"r");
    if (in == NULL) { 
      printf("Can't open line file %s\n",filen);
      return; }

    // count the number of lines
    n_lines = 0;
    int    d1,d2,d3;
    double x1,x2,x3,x4;
    while (!feof(in)) 
    {
      fscanf(in,"%lf %d %d %lf %lf %lf %d\n",&x1,&d1,&d2,&x2,&x3,&x4,&d3);
      for (int j=0;j<d3;j++) fscanf(in,"%lf %lf\n",&x1,&x2);
      n_lines++; 
    }
    fclose(in);

    // read in the lines
    l = new SINGLE_LINE[n_lines];
    in = fopen(filen,"r");
    for (int i=0;i<n_lines;i++) 
    {
      fscanf(in,"%lf %d %d %lf %lf %lf %d\n",&x1,&d1,&d2,&x2,&x3,&x4,&d3);
      l[i].lam  = x1;
      l[i].A_ul = x2;
      l[i].gl   = x3;
      l[i].gu   = x4;
      l[i].atom = d1;
      l[i].ion  = d2;
      l[i].n_branch = d3;

      // oscillator strength (Rutten page 24)
      l[i].f_lu = l[i].A_ul*l[i].gu/l[i].gl/6.67e13*pow(l[i].lam*ANGS_TO_NM,2.0);
      // frequency of transition
      l[i].nu   = C_LIGHT/(l[i].lam*ANGS_TO_CM);
      // abundance
      l[i].X    = 1.0; 
      l[i].fX   = l[i].X*CLASSICAL_CS*l[i].f_lu/l[i].nu;
      
      // read in branching lines
      l[i].blam  = new double[l[i].n_branch];
      l[i].bprob = new double[l[i].n_branch];

      double norm = l[i].A_ul;
      for (int j=0;j<l[i].n_branch;j++) {
	fscanf(in,"%lf %lf\n",&x1,&x2);
	l[i].blam[j]  = x1;
	l[i].bprob[j] = x2;
	norm += x2;  }
      
      // renormalize
      for (int j=0;j<l[i].n_branch;j++) l[i].bprob[j] /= norm;
      l[i].P_scat = l[i].A_ul/norm;

      // printout to check
      //printf("%e %e %e %e\n",l[i].lam,l[i].A_ul,l[i].P_scat,l[i].f_lu);
      //for (int j=0;j<l[i].n_branch;j++) printf("%e %e\n",l[i].blam[j],l[i].bprob[j]);
   }
  }



};

#endif
