#include "helmholtz.h"

double force0 ( double x, double y ) {
	return 2*((1-6*x*x)*y*y*(1-y*y)+(1-6*y*y)*x*x*(1-x*x));
}
double exact0 ( double x, double y ) {
	return (x*x-x*x*x*x)*(y*y*y*y-y*y);
}

double force1 ( double x, double y ) {
	double expy = exp(x+y);
	double eymx = exp(y-x);
	return -M_PI*expy*(cos(M_PI*expy/2)-M_PI*expy*sin(M_PI*expy/2)/2) + M_PI*eymx*(sin(M_PI*eymx/2)+M_PI*eymx*cos(M_PI*eymx/2)/2);
}
double exact1 ( double x, double y ) {
	return cos(M_PI*exp(y-x)/2) + sin(M_PI*exp(x+y)/2);
}

double force2 ( double x, double y ) {
	return (x==0) ? 1.0 : 0.0;
}
double exact2 ( double x, double y ) {
	return 0.0;
}