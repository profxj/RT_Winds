#ifndef _LOCATE_ARRAY_H
#define _LOCATE_ARRAY_H 1

class LOCATE_ARRAY {
 
  int n_elements;
  double *low;
  double *high;

public:

  // constructors
  LOCATE_ARRAY()  {n_elements = 0;}
  LOCATE_ARRAY(int n) {New(n);}
  LOCATE_ARRAY(int n,double *bounds) {New(n,bounds);}

  // Return
  int    Size()        {return n_elements;}
  double Center(int i) {return 0.5*(low[i] + high[i]);}
  double Low(int i)    {return low[i];}
  double High(int i)   {return high[i];}
  double Max()         {return high[n_elements-1];}
  double Min()         {return low[0];}
  double Span()        {return high[n_elements-1] - low[0];} 
  double Delta(int i)  {if (n_elements == 1) return 1.0;
                        else return high[i]-low[i]; }

  void New(int);
  void New(const char* filen);
  void New(double,double,double);
  void Uniform_New(double,double,int);
  void Log_New(double,double,int);
  void New(int,double*);
  void Copy(LOCATE_ARRAY);
  void Renormalize(double);
  int  Locate(double);
  int  Locate_Bounded(double);
  double Interpolate(double*,double);
  double Sample(int);
  void Print();
  void Print(int);
};

#endif
