#ifndef _SPECTRUM_H
#define _SPECTRUM_H 1

#include <string>
#include "locate_array.hh"
using std::string;

// default values
#define DEFAULT_NAME "optical"

class SPECTRUM {
 
private:

  // spectrum name
  char name[1000];
  string datadir;

  // number of elements and index helpers
  int n_elements;
  int n_times,n_lambda,n_x,n_y;
  int a1,a2,a3;

  // bin arrays
  LOCATE_ARRAY time_grid;
  LOCATE_ARRAY lambda_grid;
  LOCATE_ARRAY x_grid;
  LOCATE_ARRAY y_grid;
  
  // counting arrays
  double *click,*count; 
  double n_escaped;




  // Return total index given index of time,lambda,mu,phi
  int Index(int,int,int,int);

  // Used for MPI
  void MPI_Allreduce_Array(double *arr);
    
public:

  // constructors
  SPECTRUM();
  
  // Initialize
  void Init(double*, double*, int, double);
  void Init(double*, int, double);
  void Set_Name(char *n);
  void Set_Data_Directory(string f) {  datadir = f; }

  void Count(double, double, double, double, double);
    
  double lambda_min() { return lambda_grid.Min(); }
  double lambda_max() { return lambda_grid.Max(); }
  double lambda_size() { return lambda_grid.Size(); }

  // MPI functions
  void MPI_Sum_All();
  void MPI_Average_All();
  
  // Count and Normalize counted packets
  //void Count(PACKET p);
  void Normalize();
  void Rescale(double);
  void Get_Magnitude(int t, int m, double *result);
  void Bolometric_Luminosity(int t, int m, double *result);

  // Print out
  void Print();
  void Print_Spectrum();
  void Print_Lightcurve();

  double& operator()(int m, int p, int l, int s);
  double& operator()(int m, int p, int l);
  
  void Wipe();
  double N_Escaped() {return n_escaped;}
};

#endif
