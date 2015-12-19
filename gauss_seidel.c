#include "smoothers.h"

int gauss_seidel (
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
	double *residual = (double*) malloc (n*sizeof(double));

	struct timespec tsi={0,0}, tsf={0,0};
	clock_gettime (CLOCK_MONOTONIC, &tsi);

	int step = 0;
	do {
		for (int row=0; row<n; row++) {
			double sigma = 0.0;
			residual[row] = b[row];
			for (int col=0; col<n; col++) {
				if (row!=col) sigma += A[row][col] * x_next[col];
				residual[row] -= A[row][col]*x_next[col];
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
	return step;
}