/************************************************************************
This function computes the optical depth 
in each cell for a line of a given transition
at a given frequency "frequency"

Each gas cell has a line centre frequency which is corrected for 
proper motion, assuming the reference frame is zero at the center
of the galaxy.

*************************************************************************/


//Parameter structure used by the integration function
struct voig_params{ 
  double a; 
  double u;
};


//Integral part of the Voig function in the approximation by zaghloul et al. 2007
double voigt (double x, void * params) {
  struct voig_params * p = (struct voig_params *)params;
  double f =exp(-(p->u*p->u-x*x))*sin(2*p->a*(p->u-x));
  return f;
}

/************************
******MAIN FUNCTION******
************************/

void tau_line(sim_spectra linetau, vector<double> * tau, double frequency, line_cons line, long ncell){
  
  //cout<<"TAU_LINE: Computing optical depth in each cell..."<<endl;
  
  //set constants (use cgs everywhere)
  double s=2.654E-2*line.oscill;
  double clight=2.9979E10;
  double rest_nu=clight/line.rest_lambda*1E8;
  double kbol=1.3807E-16;

  //some variable
  double lineceter_nu,delta_nud,phi;
  double Hau,error,ffactor;
  
  //initilize integral
  gsl_integration_workspace * w=gsl_integration_workspace_alloc (1000);
       
  //intialize function
  gsl_function F;
  F.function = &voigt;
  struct voig_params p;

  //loop over the cells
  for(long i=0; i<ncell; i++){
    
    //this is the frequency of the center line corrected for proper motion of the cell
    lineceter_nu=linetau.losvel[i]/(clight*1E-5)*rest_nu+rest_nu;
    //define constants for this cell. Note that Delta(nu_D) is at rest frequency, NOT linecenter
    delta_nud=(rest_nu/clight)*sqrt(2*kbol*linetau.temp[i]/line.mass+line.turbul*line.turbul);
    p.a=(double) line.gamma/(12.566371*delta_nud);
    p.u=(double) (frequency-lineceter_nu)/delta_nud;
        
    //Use the approximation by zaghloul 2007
    //I checked with MW3 at z=2.3 against Voigt in IDL and I have a 
    //relative error less than 4E-5.

    if(p.a < 200){
      //compute the inegral term for H(a,u)
      F.params=&p;
      //gsl_integration_qags(&F,0.,p.u,0,1e-8,1000,w,&Hau,&error); 
      gsl_integration_qag(&F,0.,p.u,0,1e-8,1000,6,w,&Hau,&error); 
      //compute the factor in the front (2 cases)
      if(p.a < 26.6){
	//use full function
	ffactor=exp(-p.u*p.u)*cos(2*p.a*p.u)*exp(p.a*p.a)*erfc(p.a);
	Hau=ffactor+1.12837917*Hau;
      } else {
	//use series expansion
	ffactor=exp(-p.u*p.u)*cos(2*p.a*p.u)/(1.7724538*p.a)*
	  (1.-(0.5/pow(p.a,2.))+(0.75/pow(p.a,4.))-(1.875/pow(p.a,6.))+
	   (6.5625/pow(p.a,8.))-(29.53125/pow(p.a,10.))+(162.421875/pow(p.a,12.)));
	Hau=ffactor+1.12837917*Hau;
	  }
    } else {
      //for large a, use approximate equation
      Hau=p.a/(p.a*p.a+p.u*p.u)/1.7724538;
    }

    //get optical depth
    phi=Hau/(delta_nud*1.77245);
    //apply the ionization correction "boost"
    tau->push_back(phi*s*linetau.dens[i]*linetau.x_RT[i]*
		   linetau.cell[i]*linetau.boost[i]*3.0857E18);
  }
  
  //free integral
  gsl_integration_workspace_free (w);
  
  //cout<<"TAU_LINE: Done with opacity"<<endl;
  
  return;
}


