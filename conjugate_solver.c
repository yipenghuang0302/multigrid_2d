#include "solvers.h"

int conjugate_solver (
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

	// printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nINSIDE CG SOLVER\nA=\n");
	// print_sq_mat ( n, A );
	// printf("\nb=\n");
	// print_vec ( n, b );
	// printf("\nx_next=\n");
	// print_vec ( n, x_next );
	// printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

	int row, col;
	double *residual, *dee;

	residual = (double*) malloc (n*sizeof(double));
	dee = (double*) malloc (n*sizeof(double));

	// x is the vector that minimizes f(x) = 0.5(x^T)Ax - ((b)^T)x
	// initialize

	struct timespec tsi={0,0}, tsf={0,0};
	clock_gettime (CLOCK_MONOTONIC, &tsi);

	int step = 0;
	for (row=0; row<n; row++) {
		residual[row] = b[row];
		// residual[row] -= (2.0-kh*kh)*x_next[row];
		// if (0<row) residual[row] -= -1.0*x_next[row-1];
		// if (row<n-1) residual[row] -= -1.0*x_next[row+1];
		for (col=0; col<n; col++)
			residual[row] -= A[row][col]*x_next[col]; // in HCDC, x elements are non constant all but row number of terms are constant
		dee[row] = residual[row];
	}

	// delta = rTr;
	*delta_l2 = vec_vec_mult (n,residual,residual);

	do {
		double alpha = *delta_l2 / (vec_mat_vec_mult (n, dee, A) + DBL_EPSILON);
		for (row=0; row<n; row++)
			x_next[row] += alpha * dee[row];
		for (row=0; row<n; row++) {
			residual[row] = b[row];
			// residual[row] -= (2.0-kh*kh)*x_next[row];
			// if (0<row) residual[row] -= -1.0*x_next[row-1];
			// if (row<n-1) residual[row] -= -1.0*x_next[row+1];
			for (col=0; col<n; col++)
				residual[row] -= A[row][col]*x_next[col]; // in HCDC, x elements are non constant all but row number of terms are constant
		}
		double delta_old = *delta_l2;
		*delta_l2 = vec_vec_mult (n,residual,residual);
		double beta = *delta_l2 / (delta_old);// + DBL_EPSILON);
		for (row=0; row<n; row++) {
			dee[row] = residual[row] + beta*dee[row];
		}
		step++;
		// printf("step = %d\n",step++);
	} while (
		*delta_l2>DBL_EPSILON
		&& step<max_step
	);

	clock_gettime(CLOCK_MONOTONIC, &tsf);
	*time = ((double)tsf.tv_sec + 1.0e-9*tsf.tv_nsec) - ((double)tsi.tv_sec + 1.0e-9*tsi.tv_nsec);

	free (residual);
	free (dee);

	return step;

}