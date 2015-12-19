#include "2d.h"

int main ( );

void problem (
	double force (
		double x,
		double y
	),
	double k,
	double exact (
		double x,
		double y
	)
);

bool trial (
	int intervals,
	char coarse_gamma,
	int parallel_size,
	int ddm_steps,
	int relax_steps,
	double force (
		double x,
		double y
	),
	double k,
	int relax (
		int n,
		char resolution,
		int max_step,
		double step_size,
		double ** A,
		double * b,
		double * x_next,
		double * delta_l2,
		double * time
	),
	char resolution,
	double step_size,
	double exact (
		double x,
		double y
	)
);

double force0 ( double x, double y );
double exact0 ( double x, double y );

double force1 ( double x, double y );
double exact1 ( double x, double y );

double force2 ( double x, double y );
double exact2 ( double x, double y );