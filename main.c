#define HELMHOLTZ_2D

#include "helmholtz.h"

int main () {

	printf ( "2d helmholtz multigrid test:\n" );

	printf ( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	printf ( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	printf ( "-u''(x,y) - k^2 u(x,y) = 2*((1-6*x*x)*y*y*(1-y*y)+(1-6*y*y)*x*x*(1-x*x)), for -1<x<1, -1<y<1\n" );
	printf ( "Solution is u(x,y) = (x*x-x*x*x*x)*(y*y*y*y-y*y)\n" );
	problem (
		force0,
		0,//3125,
		exact0
	);

	// printf ( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	// printf ( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	// printf ( "-u''(x,y) - k^2 u(x,y) = -pi*exp(x+y)*(cos(pi*exp(x+y)/2)-pi*exp(x+y)*sin(pi*exp(x+y)/2)/2) + pi*exp(y-x)*(sin(pi*exp(y-x)/2)+pi*exp(y-x)*cos(pi*exp(y-x)/2)/2)\n" );
	// printf ( "Solution is u(x,y) = cos(pi*exp(y-x)/2) + sin(pi*exp(x+y)/2)\n" );
	// problem (
	// 	force1,
	// 	0.0,//3125,
	// 	exact1
	// );

	return 0;

}