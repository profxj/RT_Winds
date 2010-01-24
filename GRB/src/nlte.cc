// some useful transforamtions when reading in the line list
void Read_Line()
{

  // get energy difference from ionization potentials
  double dEv = level[lu].E_excute - level[ll].E_excite;
  double nu  = dEv*EV_TO_ERGS/H_PLANCK;
  double lambda   = C_LIGHT/nu*CM_TO_ANGS; // in angstroms
  
  // ratio of statistical weights
  gu_over_gl = level[lu].Get_Weight()/level[ll].Get_Weight();

  // Einstein Coeficients, given A_ul
  B_ul = A_ul*C_LIGHT*C_LIGHT/2.0/H_PLANCK/nu/nu/nu; 
  B_lu = B_ul*gu_over_gl;

  // Oscillator strength, see Rutten page 24 
  f_lu = A_ul*gu_over_gl*(lambda/10)*(lambda/10)/6.67e13;
}



int NLTE_ATOM::Solve_NLTE()
{
  int n_levels;
  int n_lines;

  // rate matrix
  gsl_matrix *M_nlte;
  //   vector of level populations - to be solved for
  gsl_vector *x_nlte;
  // vector of right hand side of equation = 0 except for last row = 1
  gsl_vector *b_nlte;
  // matrix used by the solver
  gsl_permutation *p_nlte;

  // allocate matrix and vectors
  M_nlte = gsl_matrix_calloc(n_levels,n_levels);
  b_nlte = gsl_vector_calloc(n_levels);
  x_nlte = gsl_vector_calloc(n_levels);
  p_nlte = gsl_permutation_alloc(n_levels);

  // zero out matrix and vectors
  gsl_matrix_set_zero(M_nlte);
  gsl_vector_set_zero(b_nlte);
  gsl_vector_set_zero(x_nlte);
  gsl_permutation_init(p_nlte);

  // setup the rate matrix
  for (i=0;i<n_lines;i++) 
  {
    // mean intensity of the radiation field integrated
    // over line profile
    double J; 

    // we calculate 2 rates per line transition -- up and down

    // rate down
    int l1 = line[i].l_up;
    int l2 = line[i].l_dn;
    dNdT  = line[l].A_ul + line[l].B_ul*J;
    R12  = gsl_matrix_get(M_nlte,l1,l2) + dNdt;
    R22  = gsl_matrix_get(M_nlte,l2,l2) - dNdt;
    gsl_matrix_set(M_nlte,l1,l2,R12);
    gsl_matrix_set(M_nlte,l2,l2,R22);

    // rate up
    int l1 = line[i].l_dn;
    int l2 = line[i].l_up;
    dNdT  = line[l].B_lu*J;
    R12  = gsl_matrix_get(M_nlte,l1,l2) + dNdt;
    R22  = gsl_matrix_get(M_nlte,l2,l2) - dNdt;
    gsl_matrix_set(M_nlte,l1,l2,R12);
    gsl_matrix_set(M_nlte,l2,l2,R22);
  }
    
  // last row expresses number conservation
  // otherwise matrix is singular
  for (i=0;i<n_levels;i++) gsl_matrix_set(M_nlte,n_levels-1,i,1.0);
  gsl_vector_set(b_nlte,n_levels-1,1.0);
  
  // solve matrix
  int status;
  gsl_linalg_LU_decomp(M_nlte, p_nlte, &status);
  gsl_linalg_LU_solve(M_nlte, p_nlte, b_nlte, x_nlte);

  // print out solved for levels
  double pop;
  for (i=0;i<n_levels;i++)
  {
    pop = gsl_vector_get(x_nlte,i);
    printf("%d %e\n",i,pop);
  }
  
  
}
