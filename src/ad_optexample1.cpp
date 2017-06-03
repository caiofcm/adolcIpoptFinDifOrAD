
#include "ad_optexample1.h"

/*
ADCONVERT: @AD@ ad_myExample1 @AD@  ad_OptInterface
*/

ad_myExample1::ad_myExample1() : ad_OptInterface(2, 2){
	x0 = {1.234, 5.678};
	xLB[1] = 0.0;
	gLB[0] = -2e19;
	gLB[1] = -2e19;
	gUB[0] = 0.0;
	gUB[1] = 0.0;
}

void ad_myExample1::eval_f(uint n, adouble *x, adouble& obj_value){
	obj_value = sqrt(x[1]);
	return;
}
void ad_myExample1::eval_g(uint n, adouble* x, uint m, adouble* g){
	double a1 = 2.0, b1 = 0.0, a2 = -1.0, b2 = 1.0;
	g[0] = pow((a1*x[0] + b1), 3) - x[1];
	g[1] = pow((a2*x[0] + b2), 3) - x[1];
	return;
}
