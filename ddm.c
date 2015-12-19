 #include "ddm.h"

int	ddm (
	int node_fine,
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
	double kh,
	double ** b_fine,
	double ** u_next,
	double * delta_l2,
	double * relax_time,
	double * setup_time
) {

	// double block_time = 0.0;
	// relax (
	// 	node_fine,
	// 	resolution,
	// 	relax_steps,
	// 	A_fine,
	// 	b_fine,
	// 	u_next,
	// 	delta_l2,
	// 	&block_time
	// );
	// *relax_time += block_time;


	int blocks = (node_fine-2) / parallel_size;
	// double * u_curr = (double *) malloc (node_fine*sizeof(double));
	double dummy_delta;

	// printf("blocks = %d\n",blocks);
 
	// printf("\nddm: A_fine\n");
	// print_sq_mat(node_fine, A_fine);
	// printf("\n````````\n");

	// printf("\nddm: b_fine\n");
	// print_sq_mat(node_fine, b_fine);
	// printf("\n````````\n");

	// printf("\nddm: u_next\n");
	// print_sq_mat(node_fine, u_next);
	// printf("\n````````\n");

	// iterate in gauss-seidel fashion
	for (int ddm_step=0 ; ddm_step<ddm_steps ; ddm_step++) {

		// for (int row=0; row<node_fine; row++) {
		// 	u_curr[row] = u_next[row];
		// }

		int leftover = node_fine-2 - blocks*parallel_size;
		// whole rows
		for (int block_row=0 ; block_row<blocks ; block_row++) {
			// whole cols
			for (int block_col=0 ; block_col<blocks ; block_col++) {
				ddm_block (
					block_row*parallel_size+1,
					block_col*parallel_size+1,
					parallel_size,
					parallel_size,
					relax_steps,
					relax,
					resolution,
					step_size,
					kh,
					b_fine,
					// u_curr,
					u_next,
					&dummy_delta,
					relax_time,
					setup_time
				);
			} // done with whole cols
			// incomplete cols
			ddm_block (
				block_row*parallel_size+1,
				blocks*parallel_size+1,
				parallel_size,
				leftover,
				relax_steps,
				relax,
				resolution,
				step_size,
				kh,
				b_fine,
				// u_curr,
				u_next,
				&dummy_delta,
				relax_time,
				setup_time
			);
		} // done with whole rows
		// incomplete rows
		// whole cols
		for (int block_col=0 ; block_col<blocks ; block_col++) {
			ddm_block (
				blocks*parallel_size+1,
				block_col*parallel_size+1,
				leftover,
				parallel_size,
				relax_steps,
				relax,
				resolution,
				step_size,
				kh,
				b_fine,
				// u_curr,
				u_next,
				&dummy_delta,
				relax_time,
				setup_time
			);
		} // done with whole cols
		// incomplete cols
		ddm_block (
			blocks*parallel_size+1,
			blocks*parallel_size+1,
			leftover,
			leftover,
			relax_steps,
			relax,
			resolution,
			step_size,
			kh,
			b_fine,
			// u_curr,
			u_next,
			&dummy_delta,
			relax_time,
			setup_time
		);

	}

	// printf("\nresidual calc: A_fine\n");
	// print_sq_mat(node_fine, A_fine);
	// printf("\n````````\n");

	// printf("\nresidual calc: b_fine\n");
	// print_sq_mat(node_fine, b_fine);
	// printf("\n````````\n");

	// printf("\nresidual calc: u_next\n");
	// print_sq_mat(node_fine, u_next);
	// printf("\n````````\n");

	double *residual = (double*) malloc (node_fine*node_fine*sizeof(double));
	for (int row=0; row<node_fine; row++) {
		for (int col=0; col<node_fine; col++) {
			residual[row*node_fine+col] = 0.0;
		}
	}
	for (int row=1; row<node_fine-1; row++) {
		for (int col=1; col<node_fine-1; col++) {
			// printf("\nresidual calc: residual, row=%d, col=%d, row*node_fine+col=%d, b_fine[%d][%d]=%f\n", row, col, row*node_fine+col, row, col, b_fine[row][col]);
			// print_vec(node_fine*node_fine, residual);
			// printf("\n````````\n");
			residual[row*node_fine+col] = b_fine[row][col];
			// for (int col=0; col<node_fine; col++) {
			residual[row*node_fine+col] -= (4.0-kh*kh)*u_next[row][col];
			residual[row*node_fine+col] -= -1.0*u_next[row+1][col];
			residual[row*node_fine+col] -= -1.0*u_next[row-1][col];
			residual[row*node_fine+col] -= -1.0*u_next[row][col+1];
			residual[row*node_fine+col] -= -1.0*u_next[row][col-1];
			// }
		}
	}
	// printf("\nresidual calc: residual\n");
	// print_vec(node_fine*node_fine, residual);
	// printf("\n````````\n");
	*delta_l2 = vec_vec_mult (node_fine*node_fine,residual,residual);
	// printf("\n*delta_l2 = %e\n",*delta_l2);
	// printf("\n````````\n");
	free (residual);
	// free (u_next);

	return 0;
}

void ddm_block (
	int row_start_index,
	int col_start_index,
	int row_parallel_size,
	int col_parallel_size,
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
	// double ** A_fine,
	double kh,
	double ** b_fine,
	// double * u_curr,
	double ** u_next,
	double * delta_l2,
	double * relax_time,
	double * setup_time
) {

	// generate stencil stupid way
	double ** A_block = (double **) malloc (row_parallel_size*col_parallel_size*sizeof(double*));
	for (int row=0; row<row_parallel_size*col_parallel_size; row++) {
		A_block[row] = (double*) malloc (row_parallel_size*col_parallel_size*sizeof(double));
	}
	double * b_block = (double *) malloc (row_parallel_size*col_parallel_size*sizeof(double));
	double * x_block = (double *) malloc (row_parallel_size*col_parallel_size*sizeof(double));

	struct timespec tsi={0,0}, tsf={0,0};
	clock_gettime (CLOCK_MONOTONIC, &tsi);

	// create A matrix
	for (int kron_row=0; kron_row<row_parallel_size; kron_row++) {
		for (int kron_col=0; kron_col<row_parallel_size; kron_col++) {
			for (int sub_row=0; sub_row<col_parallel_size; sub_row++) {
				for (int sub_col=0; sub_col<col_parallel_size; sub_col++) {
					// printf("row_parallel_size=%d,col_parallel_size=%d,kron_row=%d,kron_col=%d,sub_row=%d,sub_col=%d\n",row_parallel_size,col_parallel_size,kron_row,kron_col,sub_row,sub_col);fflush(stdout);
					int row_index = kron_row*col_parallel_size + sub_row;
					int col_index = kron_col*col_parallel_size + sub_col;
					if (kron_row==kron_col) {
						// 4s down the diagonal
						if (sub_row==sub_col) A_block[row_index][col_index] = (4.0-kh*kh);
						// 1s to the immediate left
						else if (sub_row==sub_col-1) A_block[row_index][col_index] = -1.0;
						// 1s to the immediate right
						else if (sub_row==sub_col+1) A_block[row_index][col_index] = -1.0;
						else A_block[row_index][col_index] = 0.0;
					// 1s at col_parallel_size away
					} else if (kron_row==kron_col-1 || kron_row==kron_col+1) {
						if (sub_row==sub_col) A_block[row_index][col_index] = -1.0;
						else A_block[row_index][col_index] = 0.0;
					} else {
						A_block[row_index][col_index] = 0.0;
					}
				}
			}
		}
	}
	// printf("\nA_block\n");
	// print_sq_mat(row_parallel_size*col_parallel_size, A_block);
	// printf("\n````````\n");

	// create b vector
	for (int b_row=0; b_row<row_parallel_size; b_row++) {
		for (int b_col=0; b_col<col_parallel_size; b_col++) {
			int vector_index = col_parallel_size * b_row + b_col;
			b_block[vector_index] = b_fine[row_start_index + b_row][col_start_index + b_col];
			// top
			if (b_row==0) b_block[vector_index] -=
				-1.0 * u_next[row_start_index-1][col_start_index+b_col];
			// bottom
			if (b_row==row_parallel_size-1) b_block[vector_index] -=
				-1.0 * u_next[row_start_index+row_parallel_size][col_start_index+b_col];
			// left
			if (b_col==0) b_block[vector_index] -=
				-1.0 * u_next[row_start_index+b_row][col_start_index-1];
			// right
			if (b_col==col_parallel_size-1) b_block[vector_index] -=
				-1.0 * u_next[row_start_index+b_row][col_start_index+col_parallel_size];
		}
	}
	// printf("\nb_block\n");
	// print_vec(row_parallel_size*col_parallel_size, b_block);
	// printf("\n````````\n");

	// create x_block
	for (int x_row=0; x_row<row_parallel_size; x_row++) {
		for (int x_col=0; x_col<col_parallel_size; x_col++) {
			int vector_index = col_parallel_size * x_row + x_col;
			x_block[vector_index] = u_next[row_start_index+x_row][col_start_index+x_col];
		}
	}
	// printf("\nx_block\n");
	// print_vec(row_parallel_size*col_parallel_size, x_block);
	// printf("\n````````\n");

	clock_gettime(CLOCK_MONOTONIC, &tsf);
	*setup_time += ((double)tsf.tv_sec + 1.0e-9*tsf.tv_nsec) - ((double)tsi.tv_sec + 1.0e-9*tsi.tv_nsec);

	double block_time = 0.0;
	relax (
		row_parallel_size*col_parallel_size,
		resolution,
		relax_steps,
		step_size,
		A_block,
		b_block,
		x_block,
		delta_l2,
		&block_time
	);
	*relax_time += block_time;

	// printf("\nx_block\n");
	// print_vec(row_parallel_size*col_parallel_size, x_block);
	// printf("\n````````\n");

	// create x_block
	for (int x_row=0; x_row<row_parallel_size; x_row++) {
		for (int x_col=0; x_col<col_parallel_size; x_col++) {
			int vector_index = col_parallel_size * x_row + x_col;
			u_next[row_start_index+x_row][col_start_index+x_col] = x_block[vector_index];
		}
	}

	free (x_block);
	free (b_block);
	for (int row=0; row<row_parallel_size*col_parallel_size; row++) free (A_block[row]);
	free (A_block);
}