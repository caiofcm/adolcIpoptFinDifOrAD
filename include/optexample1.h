#include "cfcmINLP.h"
#include <math.h>

/*
ADCONVERT: @AD@ myExample1 @AD@  OptInterface
*/

class myExample1 : public OptInterface{

public:
	myExample1();
	~myExample1() {}

	int a;
	double ajsdijds; int aks;

	void eval_f(uint n, double/*a*/ *x, double/*a*/& obj_value);

	void eval_g(uint n, double/*a*/* x, uint m, double/*a*/* g);   

};