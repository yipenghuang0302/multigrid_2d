#include "2d.h"

// - U''(X) = F(X) for A < X < B,
// with boundary conditions U(A) = UA, U(B) = UB.
bool multigrid_2d (
	int interval_fine,
	double k,
	char coarse_gamma,
	int parallel_size,
	int ddm_steps,
	int relax_steps,
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
	double ** b_fine,
	double ** u_fine,
	int * multigrid_iterations,
	double * delta_l2,
	double * relax_time,
	double * setup_time
) {

	int node_fine = interval_fine + 1;

	// printf("\nmultigrid_2d b_fine\n");
	// print_vec(node_fine, b_fine);
	// printf("\n````````\n");
	// printf("\nmultigrid_2d u_fine\n");
	// print_vec(node_fine, u_fine);
	// printf("\n````````\n");

	* delta_l2 = DBL_MAX;
	while (*delta_l2 > DBL_EPSILON) {// && (*relax_time+*setup_time)<(1.0/(double)(1<<1))) {

		ddm (
			node_fine,
			parallel_size,
			ddm_steps,
			relax_steps,
			relax,
			resolution,
			step_size,
			k/(double)interval_fine,
			b_fine,
			u_fine,
			delta_l2,
			relax_time,
			setup_time
		);

		// recurse if coarser available
		if (interval_fine>2) for (char coarse=0; coarse<coarse_gamma; coarse++) {

			// printf("interval_fine = %d\n", interval_fine);
			// printf("coarse = %d\n", coarse);
			// printf("coarse_gamma = %d\n", coarse_gamma);

			int node_coarse = interval_fine/2 + 1;
			double ** u_coarse = (double **) malloc ((node_coarse) * sizeof(double *));
			double ** residual_coarse = (double **) malloc ((node_coarse) * sizeof(double *));
			for (int row=0; row<node_coarse; row++) {
				u_coarse[row] = (double *) malloc ((node_coarse) * sizeof(double));
				residual_coarse[row] = (double *) malloc ((node_coarse) * sizeof(double));
			}

			// printf("done malloc\n");fflush(stdout);

			ftoc_2d (
				node_fine,
				k/(double)interval_fine,
				b_fine,
				u_fine,
				node_coarse,
				u_coarse,
				residual_coarse,
				setup_time
			);

			// Solve A_2h E_2h = r_2h
			// multigrid recursions
			multigrid_2d (
				interval_fine/2,
				k,
				coarse_gamma,
				parallel_size,
				ddm_steps,
				relax_steps,
				relax,
				resolution,
				step_size,
				residual_coarse,
				u_coarse,
				multigrid_iterations,
				delta_l2,
				relax_time,
				setup_time
			);

			// Interpolate E_2h back to E_h = I^h_2h E_2h.
			// Add E_h to u_h.
			struct timespec tsi={0,0}, tsf={0,0};
			clock_gettime (CLOCK_MONOTONIC, &tsi);

			ctof_2d (
				node_coarse,
				u_coarse,
				u_fine
			);

			clock_gettime(CLOCK_MONOTONIC, &tsf);
			*setup_time += ((double)tsf.tv_sec + 1.0e-9*tsf.tv_nsec) - ((double)tsi.tv_sec + 1.0e-9*tsi.tv_nsec);

			for (int row=0; row<node_coarse; row++) {
				free (residual_coarse[row]);
				free (u_coarse[row]);
			}
			free (residual_coarse);
			free (u_coarse);

		}

		// Iterate on A_h u = b_h starting from the improved u_h + E_h.
		ddm (
			node_fine,
			parallel_size,
			ddm_steps,
			relax_steps,
			relax,
			resolution,
			step_size,
			k/(double)interval_fine,
			b_fine,
			u_fine,
			delta_l2,
			relax_time,
			setup_time
		);

		// printf("\nRESULT: u_fine\n");
		// print_vec(node_fine, u_fine);
		// printf("\n````````\n");
		(*multigrid_iterations)++;
	}

	// for (int row=0; row<node_fine; row++) free (A_fine[row]);
	// free (A_fine);

	return (*delta_l2 < DBL_EPSILON);

}