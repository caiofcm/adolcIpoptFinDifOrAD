/*----------------------------------------------------------------------------
 ADOLC-C with Ipopt Wrapper: Finite Differences or Automatic Differentiation
 File: mainNLPTest.cpp
 Contents: example of using the ADOC/IPOPT wrapper for ease optimation with
 finite diffences or using automatic differentiation. It was used the command
 line tool for convert a regular double based file (obj function and constraint)
 to the ADOLC type adouble. See `convertAdolc.py`
 Author: Caio Marcellos

 Based on example file from ADOLC/IPOPT:
 ADOL-C -- Automatic Differentiation by Overloading in C++
 File:     cpp_example.cpp
 Revision: $$
 Contents: example for class myADOLC_NPL for interfacing with Ipopt
 Copyright (c) Andrea Walther
 This code is based on the file corresponding file cpp_example.cpp contained 
 in the Ipopt package with the authors:  Carl Laird, Andreas Waechter   
----------------------------------------------------------------------------*/

#include "optexample1.h"
#include "ad_optexample1.h"
#include "cfcmINLP.h"

using namespace Ipopt;
using namespace std;

int main(int argv, char* argc[]){

  // ####################################
  // Ipopt Finite Differences version:
  // ####################################
  myExample1 aow;
  OptInterface *optI = &aow;
  SmartPtr<TNLP> nlp = new cfcmINLP(optI);

  // Define ipopt options:
  SmartPtr<IpoptApplication> app = new IpoptApplication();
  app->Options()->SetStringValue("hessian_approximation", "limited-memory"); //obeyed for finite diff.

  createAppAndRun(nlp, app);

  // ####################################
  // ADOLC Ipopt version:
  // ####################################
  ad_myExample1 ad_aow;
  ad_OptInterface *ad_optI = &ad_aow;
  SmartPtr<TNLP> ad_nlp = new cfcmINLP(ad_optI);

  // Define ipopt options:
  SmartPtr<IpoptApplication> ad_app = new IpoptApplication();
  // app->Options()->SetStringValue("hessian_approximation", "limited-memory");

  createAppAndRun(ad_nlp, ad_app);
}