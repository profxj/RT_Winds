#ifndef _CDF_H
#define _CDF_H 1

//**********************************************************
// CDF == Comulative Distribution Function
//
// This simple class just holds a vector which should be
// monitonically increasing and sums up to one
// We can sample from it using a binary search.
// The class has just a few functions;
//   void New(int n) 
//   void Normalize() 
//   int  Sample();
// And these are defined right below.
//**********************************************************

class CDF 
{

private:
  
  // data stored, just a vector
  int n_elements;
  double *y;
  // MPI Reduce helper
  void MPI_Sum_Array(double *array);
  
public:

  // constructor
  CDF() {n_elements = 0;}
  CDF(int n) {New(n);}
  ~CDF();

  void New(int n);

  double Value(int i);                   // Get the local value
  double Get(int i)  {return y[i];}      // Get local CDF value

  void Set(int c, double f)  {y[c] = f;} // Set cell CDF value 
  void Set(double *f);                   // Set cell with array
  void Set(double *f,double *w);         // Set cell with array&weight
  void Uniform_UpTo(int c);              // Set uniform

  void Put(int c, double f)
  { 
    if (c == 0)   y[0] = f;
    else y[c] = y[c-1] + f; 
  }

  void Normalize(); 
  int Sample();
  void Print(); 
  void Wipe();
  void MPI_Combine(); 

};

class CDF_FLOAT 
{

private:
  
  // data stored, just a vector
  int n_elements;
  float *y;
  // MPI Reduce helper
  void MPI_Sum_Float_Array(float *array);
  
public:

  // constructor
  CDF_FLOAT() {n_elements = 0;}
  CDF_FLOAT(int n) {New(n);}
  ~CDF_FLOAT();

  void New(int n);

  float Value(int i);                   // Get the local value
  float Get(int i)  {return y[i];}      // Get local CDF value

  void Set(int c, double f)  {y[c] = f;} // Set cell CDF value 
  
  void Put(int c, double f)
  { 
    if (c == 0)   y[0] = f;
    else y[c] = y[c-1] + f; 
  }

  void Set(double *f);                   // Set cell with array
  void Set(double *f,double *w);         // Set cell with array&weight
  void Uniform_UpTo(int c);              // Set uniform

  void Normalize(); 
  int Sample();
  void Print(); 
  void Wipe();
  void MPI_Combine(); 

};
  
#endif
