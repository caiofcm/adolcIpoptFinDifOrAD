
#include "cfcmINLP.h"

// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: cfcmINLP.cpp 2005 2011-06-06 12:55:16Z stefan $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include <cassert>
#include <iostream>

using namespace Ipopt;

// OptInterface::OptInterface(uint nn, uint mm = 0){
//   n = nn;
//   m = mm;
//   for (uint i=0; i<n; i++)
//     x_l[i] = popt->xLB[i];
//   xLB =
// }
// OptInterface::OptInterface()
baseIOPt::baseIOPt(uint n, uint m){
  xLB.assign(n, -1e19);
  xUB.assign(n, 1e19);
  gLB.assign(n, -1e19);
  gUB.assign(n, 1e19);
  this->n = n;
  this->m = m;
}

OptInterface::OptInterface(uint n, uint m) : baseIOPt(n, m){
}

ad_OptInterface::ad_OptInterface(uint n, uint m ) : baseIOPt(n, m){
}

int createAppAndRun(SmartPtr<TNLP> &myadolc_nlp){
  SmartPtr<IpoptApplication> app = new IpoptApplication();
  app->Options()->SetStringValue("hessian_approximation", "limited-memory"); 
  createAppAndRun(myadolc_nlp, app);
  return 0;
}
int createAppAndRun(SmartPtr<TNLP> &myadolc_nlp, SmartPtr<IpoptApplication> &app){
  // Create an instance of your nlp...
  // SmartPtr<TNLP> myadolc_nlp = new HS071_NLP();

  // Create an instance of the IpoptApplication
  // if (app == NULL)
  // SmartPtr<IpoptApplication> app = new IpoptApplication();

  // Options for IPOPT
  // app->Options()->SetStringValue("hessian_approximation", "limited-memory");

  // Initialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\n\n*** Error during initialization!\n");
    return (int) status;
  }

  status = app->OptimizeTNLP(myadolc_nlp);

  if (status == Solve_Succeeded) {
    // Retrieve some statistics about the solve
    Index iter_count = app->Statistics()->IterationCount();
    printf("\n\n*** The problem solved in %d iterations!\n", iter_count);

    Number final_obj = app->Statistics()->FinalObjective();
    printf("\n\n*** The final value of the objective function is %e.\n", final_obj);
  }

  return (int) status;
}


// ######################################
// ######################################
// ######################################
// ######################################

// constructor
cfcmINLP::cfcmINLP(OptInterface *popt_)
{
  this->popt = popt_;
  potm = this->popt;
}
cfcmINLP::cfcmINLP(ad_OptInterface *adopt_)
{
  this->adopt = adopt_;
  potm = this->adopt;
}

//destructor
cfcmINLP::~cfcmINLP(){
  // if (xcpy != NULL)
  //   delete[] xcpy;
  // xcpy = NULL;
}

// returns the size of the problem
bool cfcmINLP::get_nlp_info(Index &n, Index &m, Index &nnz_jac_g,
                            Index &nnz_h_lag, IndexStyleEnum &index_style)
{
  // The problem has kg and kb
  n = potm->n;

  // number of constraints
  m = potm->m;

  // (considering as dense)
  nnz_jac_g = n * m;

  // the hessian is also dense and has 16 total nonzeros, but we
  // only need the lower left corner (since it is symmetric)
  nnz_h_lag = (n) * (1 + (n)) / 2;

  // use the C style indexing (0-based)
  index_style = TNLP::C_STYLE;

  // Initialize auxiliary variables:
  if (this->popt != NULL){
    xcpy = new double[n];
    gauxM = new double[m];
    gauxm = new double[m];
    htol = 1e-7; 
  }
  Jac = new double*[m];
  for(Index i=0;i<m;i++)
    Jac[i] = new double[n];   

  if (this->adopt != NULL)
    generate_tapes(n, m);

  return true;
}

// returns the variable bounds
bool cfcmINLP::get_bounds_info(Index n, Number *x_l, Number *x_u,
                               Index m, Number *g_l, Number *g_u){
  // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
  // If desired, we could assert to make sure they are what we think they are.
  assert(n == potm->n);
  assert(m == potm->m);

  // Setting Bounds
  for (Index i = 0; i < n; i++)
    x_l[i] = potm->xLB[i];
  for (Index i = 0; i < n; i++)
    x_u[i] = potm->xUB[i];
  for (Index i = 0; i < m; i++)
    g_l[i] = potm->gLB[i];
  for (Index i = 0; i < m; i++)
    g_u[i] = potm->gUB[i];
  return true;
}

// returns the initial point for the problem
bool cfcmINLP::get_starting_point(Index n, bool init_x, Number *x,
                                  bool init_z, Number *z_L, Number *z_U,
                                  Index m, bool init_lambda,
                                  Number *lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // initialize to the given starting point
  for (Index i = 0; i < n; i++)
    x[i] = potm->x0[i];
  return true;
}

// returns the value of the objective function
bool cfcmINLP::eval_f(Index n, const Number *x, bool new_x, Number &obj_value)
{
  assert(n == potm->n);
  // CrystRawllignsModel* pDeriv = dynamic_cast<CrystRawllignsModel*>(mdl);
  // pDeriv->kb = x[0];
  // pDeriv->kg = x[1];
  // for (uint j=0;j<mdl->ny;j++)
  //   mdl->pyiter[j] = mdl->y0[j];
  // mdl->sim(mdl->pk1, mdl->pyiter, mdl->prhs, mdl->ny, mdl->tvector, mdl->pmvspan, mdl->pdsim);

  // // Least Squares
  // uint nt = mdl->tvector.size();
  // double Jaux = 0.0, Ccalc, Cmes;
  // double sigma = 1.0; // check where to put it.
  // for (uint j=1;j<nt;j++){
  //   Ccalc = mdl->pdsim[basefunc::ele(j, 0, mdl->ny)];
  //   Cmes = prest->dataEXP[j];
  //   Jaux += 1/sigma * (Ccalc - Cmes) * (Ccalc - Cmes);
  // }
  // obj_value = Jaux;
  if (this->popt != NULL) //stupid error -> Requires the brackets!
    popt->eval_f(n, const_cast<double*>(x), obj_value);
  else if (this->adopt != NULL)
    function(tag_f,1,n,const_cast<double*>(x), &obj_value);

  return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool cfcmINLP::eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f)
{
  assert(n == potm->n);

  // Finite Difference method:
  if (this->popt != NULL){
    double J_M = 0.0, J_m = 0.0;
    std::copy(x, x + popt->n, xcpy);
    for (uint j = 0; j < popt->n; j++){
      xcpy[j] += htol;
      eval_f(n, xcpy, new_x, J_M);
      xcpy[j] = x[j];
      xcpy[j] -= htol;
      eval_f(n, xcpy, new_x, J_m);
      grad_f[j] = (J_M - J_m) / (2.0 * htol);
      xcpy[j] = x[j];
    }
  }
  else if (this->adopt != NULL){
    gradient(tag_f,n,x,grad_f);
  }

  return true;
}

// return the value of the constraints: g(x)
bool cfcmINLP::eval_g(Index n, const Number *x, bool new_x, Index m, Number *g){
  assert(n == potm->n);
  // assert(m == 2);

  if (this->popt != NULL)
    popt->eval_g(n, const_cast<double*>(x), m, const_cast<double*>(g));
  else if (this->adopt != NULL)
    function(tag_g,m,n,const_cast<double*>(x), g);

  return true;
}

// return the structure or values of the jacobian
bool cfcmINLP::eval_jac_g(Index n, const Number *x, bool new_x,
                          Index m, Index nele_jac, Index *iRow, Index *jCol,
                          Number *values)
{
  if (values == NULL)
  {
    // return the structure of the jacobian,
    // assuming that the Jacobian is dense
    Index idx = 0;
    for (Index i = 0; i < m; i++)
      for (Index j = 0; j < n; j++)
      {
        iRow[idx] = i;
        jCol[idx++] = j;
      }
  }
  else
  {
    // Finite Difference method:
    if (this->popt != NULL)
    {
      // double Jacij = 0.0;
      std::copy(x, x + potm->n, xcpy);
      for (uint j = 0; j < potm->n; j++)
      {
          xcpy[j] += htol;
          eval_g(n, xcpy, new_x, m, gauxM);
          xcpy[j] = x[j];
          xcpy[j] -= htol;
          eval_g(n, xcpy, new_x, m, gauxm);
          xcpy[j] = x[j];
          for (uint i = 0; i < potm->m; i++){
            values[i*n + j] = (gauxM[i] - gauxm[i]) / (2.0 * htol);
          }
      }
    }
    else if (this->adopt != NULL){
      jacobian(tag_g,m,n,x,Jac);
      Index idx = 0;
      for(Index i=0; i<m; i++)
        for(Index j=0; j<n; j++)
      values[idx++] = Jac[i][j];
    }
  }

  return true;
}

//return the structure or values of the hessian
bool cfcmINLP::eval_h(Index n, const Number *x, bool new_x,
                      Number obj_factor, Index m, const Number *lambda,
                      bool new_lambda, Index nele_hess, Index *iRow,
                      Index *jCol, Number *values){
  if (values == NULL)
  {
    // return the structure. This is a symmetric matrix, fill the lower left
    // triangle only.

    // the hessian for this problem is actually dense
    Index idx = 0;
    for (Index row = 0; row < n; row++)
    {
      for (Index col = 0; col <= row; col++)
      {
        iRow[idx] = row;
        jCol[idx] = col;
        idx++;
      }
    }

    assert(idx == nele_hess);
  }
  else
  {
    if (this->adopt != NULL)
    {
      obj_lam[0] = obj_factor;
      for (Index i = 0; i < m; i++)
        obj_lam[1 + i] = lambda[i];

      set_param_vec(tag_L, m + 1, obj_lam);
      hessian(tag_L, n, const_cast<double *>(x), Hess);
      Index idx = 0;
      for (Index i = 0; i < n; i++)
        for (Index j = 0; j <= i; j++)
          values[idx++] = Hess[i][j];
    }
  }

  return true;
}

void cfcmINLP::finalize_solution(SolverReturn status,
                                 Index n, const Number *x, const Number *z_L, const Number *z_U,
                                 Index m, const Number *g, const Number *lambda,
                                 Number obj_value,
                                 const IpoptData *ip_data,
                                 IpoptCalculatedQuantities *ip_cq){
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  // For this example, we write the solution to the console
  std::cout << std::endl
            << std::endl
            << "Solution of the primal variables, x" << std::endl;
  for (Index i = 0; i < n; i++)
  {
    std::cout << "x[" << i << "] = " << x[i] << std::endl;
  }

  // std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
  // for (Index i=0; i<n; i++) {
  //   std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
  // }
  // for (Index i=0; i<n; i++) {
  //   std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  // }

  std::cout << std::endl
            << std::endl
            << "Objective value" << std::endl;
  std::cout << "f(x*) = " << obj_value << std::endl;

  if (m > 0){
    std::cout << std::endl
              << "Final value of the constraints:" << std::endl;
  }
  for (Index i = 0; i < m; i++){
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  }

  // Saving otm result:
  potm->xOtm.resize(n);
  potm->gOtm.resize(m);
  for (uint i = 0; i < potm->n; i++)
    potm->xOtm[i] = x[i];
  for (uint i = 0; i < potm->m; i++)
    potm->gOtm[i] = g[i];
  potm->FOtm = obj_value;

  // Deallocations:
  if (this->popt != NULL){
    delete[] xcpy;  
    delete[] gauxM;  
    delete[] gauxm;
  }
  for(Index i=0;i<m;i++)
    delete[] Jac[i];
  delete[] Jac; 
  if (this->adopt != NULL){
    for(Index i=0;i<n;i++)
      delete[] Hess[i];
    delete[] Hess;   
  }
}

//***************    ADOL-C part ***********************************

void cfcmINLP::generate_tapes(Index n, Index m)
{
  Number *xp    = new double[n];
  Number *lamp  = new double[m];
  Number *zl    = new double[m];
  Number *zu    = new double[m];

  adouble *xa   = new adouble[n];
  adouble *g    = new adouble[m];
  double *lam   = new double[m];
  double sig;
  adouble obj_value;
  
  double dummy;

  Jac = new double*[m];
  for(Index i=0;i<m;i++)
    Jac[i] = new double[n];

  obj_lam   = new double[m+1];

  Hess = new double*[n];
  for(Index i=0;i<n;i++)
    Hess[i] = new double[i+1];

  get_starting_point(n, 1, xp, 0, zl, zu, m, 0, lamp);

  trace_on(tag_f);
    
    for(Index i=0;i<n;i++)
      xa[i] <<= xp[i];

    adopt->eval_f(n,xa,obj_value);

    obj_value >>= dummy;

  trace_off();
  
  trace_on(tag_g);
    
    for(Index i=0;i<n;i++)
      xa[i] <<= xp[i];

    adopt->eval_g(n,xa,m,g);


    for(Index i=0;i<m;i++)
      g[i] >>= dummy;

  trace_off();

   trace_on(tag_L);
    
    for(Index i=0;i<n;i++)
      xa[i] <<= xp[i];
    for(Index i=0;i<m;i++)
      lam[i] = 1.0;
    sig = 1.0;

    adopt->eval_f(n,xa,obj_value);

    obj_value *= mkparam(sig);
    adopt->eval_g(n,xa,m,g);
 
    for(Index i=0;i<m;i++)
      obj_value += g[i]*mkparam(lam[i]);

    obj_value >>= dummy;

  trace_off();

  delete[] xa;
  delete[] xp;
  delete[] g;
  delete[] lam;
  delete[] lamp;
  delete[] zu;
  delete[] zl;

}
