#include "cfcmINLP.h"
#include <math.h>

/*
ADCONVERT: @AD@ ad_myExample1 @AD@  ad_OptInterface
*/

class ad_myExample1 : public ad_OptInterface{

public:
	ad_myExample1();
	~ad_myExample1() {}

	int a;
	double ajsdijds; int aks;

	void eval_f(uint n, adouble *x, adouble& obj_value);

	void eval_g(uint n, adouble* x, uint m, adouble* g);   

};
