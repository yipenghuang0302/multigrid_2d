#include "2d.h"

void ftoc_2d (
	int node_fine,
	double kh,
	double ** b_fine,
	double ** u_fine,
	int node_coarse,
	double ** u_coarse,
	double ** residual_coarse,
	double * setup_time
) {

	// printf("ftoc_2d\n");fflush(stdout);

	// Restrict the residual r_h = b_h − A_h u_h to the coarse grid by r_2h = R^2h_h r_h .
	// the residual r_h = b_h − A_h u_h
	double ** residual_fine = (double **) malloc (node_fine * sizeof(double*));
	for (int row=0; row<node_fine; row++) {
		residual_fine[row] = (double *) malloc (node_fine * sizeof(double));
	}

	struct timespec tsi={0,0}, tsf={0,0};
	clock_gettime (CLOCK_MONOTONIC, &tsi);

	for (int row=0; row<node_fine; row++) {
		for (int col=0; col<node_fine; col++) {
			residual_fine[row][col] = b_fine[row][col];
			residual_fine[row][col] -= (4.0-kh*kh) * u_fine[row][col];
			if (0<row) residual_fine[row][col] += u_fine[row-1][col];
			if (row<node_fine-1) residual_fine[row][col] += u_fine[row+1][col];
			if (0<col) residual_fine[row][col] += u_fine[row][col-1];
			if (col<node_fine-1) residual_fine[row][col] += u_fine[row][col+1];
		}
	}
	// printf("\nresidual_fine\n");
	// print_sq_mat(node_fine, residual_fine);
	// printf("\n````````\n");

	// take care of u_coarse
	for (int row_coarse=0; row_coarse<node_coarse; row_coarse++)
		for (int col_coarse=0; col_coarse<node_coarse; col_coarse++)
			u_coarse[row_coarse][col_coarse] = 0.0;

	// restrict to residual_coarse
	// to the coarse grid by r_2h = R^2h_h r_h
	for (int col_coarse=0; col_coarse<node_coarse; col_coarse++) {
		residual_coarse[0][col_coarse] = 0.0;
		residual_coarse[node_coarse-1][col_coarse] = 0.0;
	}

	for (int row_coarse=1; row_coarse<node_coarse-1; row_coarse++) {
		residual_coarse[row_coarse][0] = 0.0;//residual_fine[0]*2.0/3.0 + residual_fine[1]/3.0;
		residual_coarse[row_coarse][node_coarse-1] = 0.0;//residual_fine[node_fine-1]*2.0/3.0 + residual_fine[node_fine-2]/3.0;

		int row_fine = 2*row_coarse;
		for (int col_coarse=1; col_coarse<node_coarse-1; col_coarse++) {
			int col_fine = 2*col_coarse;
			residual_coarse[row_coarse][col_coarse] =
				4.0 * residual_fine [row_fine-1] [col_fine-1] / 16.0 +
				4.0 * residual_fine [row_fine-1] [col_fine  ] / 08.0 +
				4.0 * residual_fine [row_fine-1] [col_fine+1] / 16.0 +
				4.0 * residual_fine [row_fine  ] [col_fine-1] / 08.0 +
				4.0 * residual_fine [row_fine  ] [col_fine  ] / 04.0 +
				4.0 * residual_fine [row_fine  ] [col_fine+1] / 08.0 +
				4.0 * residual_fine [row_fine+1] [col_fine-1] / 16.0 +
				4.0 * residual_fine [row_fine+1] [col_fine  ] / 08.0 +
				4.0 * residual_fine [row_fine+1] [col_fine+1] / 16.0 ;
		}
	}

	// printf("\nresidual_coarse\n");
	// print_sq_mat(node_coarse, residual_coarse);
	// printf("\n````````\n");

	clock_gettime(CLOCK_MONOTONIC, &tsf);
	*setup_time += ((double)tsf.tv_sec + 1.0e-9*tsf.tv_nsec) - ((double)tsi.tv_sec + 1.0e-9*tsi.tv_nsec);

	for (int row=0; row<node_fine; row++) free ( residual_fine[row] );
	free (residual_fine);

	return;
}