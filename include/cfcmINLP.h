// Copyright (C) 2005, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: hs071_nlp.hpp 1861 2010-12-21 21:34:47Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-09
#pragma once
#include "IpTNLP.hpp"
#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
// #include "crystmodel.h"

#include <coin/IpTNLP.hpp>
#include <adolc/adolc.h>

#define tag_f 1
#define tag_g 2
#define tag_L 3

using namespace Ipopt;

typedef unsigned int uint;
typedef	std::vector<double> vec;

struct baseIOPt{
	baseIOPt(uint n, uint m);
	~baseIOPt() {}  
	vec x0;
	vec xOtm, gOtm;
	double FOtm;
	vec dataEXP;
  vec xLB, xUB, gLB, gUB;
  uint n = 0, m = 0;
};

struct OptInterface : public baseIOPt{
	OptInterface(uint n, uint m);
	~OptInterface() {}

	// vec x0;
	// vec xOtm, gOtm;
	// double FOtm;
	// vec dataEXP;
  // vec xLB, xUB, gLB, gUB;
  // uint n = 0, m = 0;

  bool useFiniteDifference = true;

  /** Method to return the objective value */
  virtual void eval_f(uint n, double *x, double& obj_value) {}

  /** Method to return the constraint residuals */
  virtual void eval_g(uint n, double* x, uint m, double* g) {}                 
};

struct ad_OptInterface : public baseIOPt{
  ad_OptInterface(uint n, uint m);
  ~ad_OptInterface() {}
  int b = 2;

 /** Method to return the objective value */
 virtual void eval_f(uint n, adouble *x, adouble& obj_value) {}  

  /** Method to return the constraint residuals */
  virtual void eval_g(uint n, adouble* x, uint m, adouble* g) {}   
};

int createAppAndRun(SmartPtr<TNLP> &myadolc_nlp);
int createAppAndRun(SmartPtr<TNLP> &myadolc_nlp, SmartPtr<IpoptApplication> &app_);


/** C++ Example NLP for interfacing a problem with IPOPT.
 *  HS071_NLP implements a C++ example of problem 71 of the
 *  Hock-Schittkowski test suite. This example is designed to go
 *  along with the tutorial document and show how to interface
 *  with IPOPT through the TNLP interface. 
 *
 * Problem hs071 looks like this
 *
 *     min   x1*x4*(x1 + x2 + x3)  +  x3
 *     s.t.  x1*x2*x3*x4                   >=  25
 *           x1**2 + x2**2 + x3**2 + x4**2  =  40
 *           1 <=  x1,x2,x3,x4  <= 5
 *
 *     Starting point:
 *        x = (1, 5, 5, 1)
 *
 *     Optimal solution:
 *        x = (1.00000000, 4.74299963, 3.82114998, 1.37940829)
 *
 *
 */
class cfcmINLP : public TNLP
{
public:

  // My Dynamic system:
  OptInterface *popt = NULL;
  ad_OptInterface *adopt = NULL;
  baseIOPt *potm = NULL;
  double *xcpy, *gauxM, *gauxm;
  double htol;

  // ParamEst *prest;

  /** Constructor */
  cfcmINLP(OptInterface *popt_);
  cfcmINLP(ad_OptInterface *adopt_);

  /** default destructor */
  virtual ~cfcmINLP();

  /**@name Overloaded from TNLP */
  //@{
  /** Method to return some info about the nlp */
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style);

  /** Method to return the bounds for my problem */
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u);

  /** Method to return the starting point for the algorithm */
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda);

  /** Method to return the objective value */
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value);

  /** Method to return the gradient of the objective */
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f);

  /** Method to return the constraint residuals */
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g);

  /** Method to return:
   *   1) The structure of the jacobian (if "values" is NULL)
   *   2) The values of the jacobian (if "values" is not NULL)
   */
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values);

  /** Method to return:
   *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
   *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
   */
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values);

  //@}

  /** @name Solution Methods */
  //@{
  /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
				 const IpoptData* ip_data,
				 IpoptCalculatedQuantities* ip_cq);
  //@}

//***************    start ADOL-C part ***********************************

 /** Method to generate the required tapes */
  virtual void generate_tapes(Index n, Index m);

//***************    end   ADOL-C part ***********************************

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *  
   */
  //@{
  //  cfcmINLP();
  cfcmINLP(const cfcmINLP&);
  cfcmINLP& operator=(const cfcmINLP&);
  //@}

  double **Jac; 
  double *obj_lam;
  double **Hess;  
};
