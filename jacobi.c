#include "smoothers.h"

int jacobi_float (
	int n,
	char resolution,
	int max_step,
	double step_size,
	double ** A,
	double * b,
	double * x_next,
	double * delta_l2,
	double * time
) {

	double *x_curr = (double*) malloc (n*sizeof(double));
	double *residual = (double*) malloc (n*sizeof(double));

	struct timespec tsi={0,0}, tsf={0,0};
	clock_gettime (CLOCK_MONOTONIC, &tsi);

	int step = 0;
	do {
		for (int row=0; row<n; row++) x_curr[row] = x_next[row];
		for (int row=0; row<n; row++) {
			double sigma = 0.0;
			residual[row] = b[row];
			for (int col=0; col<n; col++) {
				if (row!=col) sigma += A[row][col] * x_curr[col];
				residual[row] -= A[row][col]*x_curr[col];
			}
			x_next[row] = (b[row] - sigma) / A[row][row];
		}
		*delta_l2 = vec_vec_mult (n,residual,residual);
		step++;
	} while (
		*delta_l2>DBL_EPSILON
		&& step<max_step
	);

	clock_gettime(CLOCK_MONOTONIC, &tsf);
	*time = ((double)tsf.tv_sec + 1.0e-9*tsf.tv_nsec) - ((double)tsi.tv_sec + 1.0e-9*tsi.tv_nsec);

	free (residual);
	free (x_curr);
	return step;
}

int jacobi_fixed (
	int n,
	char resolution,
	int max_step,
	double step_size,
	double ** A,
	double * b,
	double * x_next,
	double * delta_l2,
	double * time
) {

	// create fixed version of numbers
	long ** A_fix = (long**) malloc (n*sizeof(long*));
	long *b_fix = (long*) malloc (n*sizeof(long));
	long *x_fix_curr = (long*) malloc (n*sizeof(long));
	long *x_fix_next = (long*) malloc (n*sizeof(long));
	for (int row=0; row<n; row++) {
		b_fix[row] = to_fixed(resolution,b[row]);
		x_fix_next[row] = to_fixed(resolution,x_next[row]);
		A_fix[row] = (long*) malloc (n*sizeof(long));
		for (int col=0; col<n; col++) {
			A_fix[row][col] = to_fixed(resolution,A[row][col]);
		}
	}
	long *residual = (long*) malloc (n*sizeof(long));

	struct timespec tsi={0,0}, tsf={0,0};
	clock_gettime (CLOCK_MONOTONIC, &tsi);

	int step = 0;
	do {
		for (int row=0; row<n; row++) x_fix_curr[row] = x_fix_next[row];
		for (int row=0; row<n; row++) {
			long sigma = 0;
			residual[row] = b_fix[row];
			for (int col=0; col<n; col++) {
				if (row!=col) sigma += A_fix[row][col] * x_fix_curr[col];
				residual[row] -= (A_fix[row][col]*x_fix_curr[col]) >> resolution;
			}
			x_fix_next[row] = ((b_fix[row]<<resolution) - sigma) / A_fix[row][row];
		}
		step++;
	} while (
		vec_vec_mult_fixed (n,residual,residual) > 0
		&& step<max_step
	);

	clock_gettime(CLOCK_MONOTONIC, &tsf);
	*time = ((double)tsf.tv_sec + 1.0e-9*tsf.tv_nsec) - ((double)tsi.tv_sec + 1.0e-9*tsi.tv_nsec);

	*delta_l2 = (double) vec_vec_mult_fixed (n,residual,residual) / (double)(1<<2*resolution);
	for (int row=0; row<n; row++) x_next[row] = to_float(resolution,x_fix_curr[row]);

	free (residual);
	for (int row=0; row<n; row++) free(A_fix[row]);
	free (x_fix_next);
	free (x_fix_curr);
	free (b_fix);
	free (A_fix);
	// printf("*delta_l2=%f\n",*delta_l2);

	return step;

}