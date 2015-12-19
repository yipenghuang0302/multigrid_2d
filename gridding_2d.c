#include "2d.h"

void gridding_2d (
	double span,
	int intervals,
	double * x,
	double * y,
	double force ( double x, double y ),
	double exact ( double x, double y ),
	double ** r,
	double * setup_time
) {

	double h = span / ( double ) intervals;

	struct timespec tsi={0,0}, tsf={0,0};
	clock_gettime (CLOCK_MONOTONIC, &tsi);

	// boundary conditions
	for ( int row = 0; row < intervals+1; row++ ) {
		r[row][0] = h * h * exact ( x[row], y[0] );
		r[row][intervals] = h * h * exact ( x[row], y[intervals] );
	}
	for ( int col = 0; col < intervals+1; col++ ) {
		r[0][col] = h * h * exact ( x[0], y[col] );
		r[intervals][col] = h * h * exact ( x[intervals], y[col] );
	}

	for ( int row = 1; row < intervals; row++ ) {
		for ( int col = 1; col < intervals; col++ ) {
			r[row][col] = h * h * force ( x[row], y[col] );
		}
	}

	clock_gettime(CLOCK_MONOTONIC, &tsf);
	*setup_time += ((double)tsf.tv_sec + 1.0e-9*tsf.tv_nsec) - ((double)tsi.tv_sec + 1.0e-9*tsi.tv_nsec);

	return;
}