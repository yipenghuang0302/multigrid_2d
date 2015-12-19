#include <stdlib.h>
#include <time.h>
#include "blas.h"

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
);

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
	double kh,
	double ** b_fine,
	// double * u_curr,
	double ** u_next,
	double * delta_l2,
	double * relax_time,
	double * setup_time
);