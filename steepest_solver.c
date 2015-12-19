#include "solvers.h"

int steepest_solver_float (
	int n,
	char resolution,
	int max_step,
	double step_size,
	double ** A,
	double * b,
	double * x,
	double * delta_l2,
	double * time
) {

	double *residual = (double*) malloc (n*sizeof(double));
	*delta_l2 = DBL_MAX;

	struct timespec tsi={0,0}, tsf={0,0};
	clock_gettime (CLOCK_MONOTONIC, &tsi);

	int step = 0;
	while (
		*delta_l2>DBL_EPSILON
		&& step<max_step
	) {
		for (int row=0; row<n; row++) {
			residual[row] = b[row];
			for (int col=0; col<n; col++)
				residual[row] -= A[row][col]*x[col];
		}
		*delta_l2 = vec_vec_mult (n,residual,residual);
		double alpha = *delta_l2 / (vec_mat_vec_mult (n, residual, A)+DBL_EPSILON) ;
		for (int row=0; row<n; row++)
			x[row] += alpha * residual[row];
		step++;
	}

	clock_gettime(CLOCK_MONOTONIC, &tsf);
	*time = ((double)tsf.tv_sec + 1.0e-9*tsf.tv_nsec) - ((double)tsi.tv_sec + 1.0e-9*tsi.tv_nsec);

	free (residual);

	return step;

}

int steepest_solver_fixed (
	int n,
	char resolution,
	int max_step,
	double step_size,
	double ** A,
	double * b,
	double * x,
	double * delta_l2,
	double * time
) {

	// create fixed version of numbers
	long ** A_fix = (long**) malloc (n*sizeof(long*));
	long *b_fix = (long*) malloc (n*sizeof(long));
	long *x_fix = (long*) malloc (n*sizeof(long));
	for (int row=0; row<n; row++) {
		b_fix[row] = to_fixed(resolution,b[row]);
		x_fix[row] = to_fixed(resolution,x[row]);
		A_fix[row] = (long*) malloc (n*sizeof(long));
		for (int col=0; col<n; col++) {
			A_fix[row][col] = to_fixed(resolution,A[row][col]);
		}
	}
	long *residual = (long*) malloc (n*sizeof(long));

	// delta = rTr;
	long delta = 1;

	struct timespec tsi={0,0}, tsf={0,0};
	clock_gettime (CLOCK_MONOTONIC, &tsi);

	int step = 0;
	while (
		delta > 0
		&& step<max_step
	) {
		for (int row=0; row<n; row++) {
			residual[row] = b_fix[row];
			// printf("residual[%d]=%ld\n",row,residual[row]);
			for (int col=0; col<n; col++) {
				residual[row] -= (A_fix[row][col]*x_fix[col]) >> resolution;
				// printf("residual[%d]=%ld\n",row,residual[row]);
			}
		}
		delta = vec_vec_mult_fixed (n,residual,residual);
		// printf("delta=%ld\n",delta);
		for (int row=0; row<n; row++) {
			x_fix[row] += ((delta*residual[row])<<resolution) / (vec_mat_vec_mult_fixed (n, residual, A_fix)+1);
			// printf("x_fix[%d]=%ld\n",row,x_fix[row]);
		}
		step++;
	}

	clock_gettime(CLOCK_MONOTONIC, &tsf);
	*time = ((double)tsf.tv_sec + 1.0e-9*tsf.tv_nsec) - ((double)tsi.tv_sec + 1.0e-9*tsi.tv_nsec);

	*delta_l2 = ((double)delta) / (double)(1<<2*resolution);
	for (int row=0; row<n; row++) x[row] = to_float(resolution,x_fix[row]);

	free (residual);
	for (int row=0; row<n; row++) free(A_fix[row]);
	free (x_fix);
	free (b_fix);
	free (A_fix);
	// printf("*delta_l2=%f\n",*delta_l2);

	return step;

}