#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "ddm.h"
#include "../relax/solvers.h"
#include "../relax/smoothers.h"

double * spacing (
  int n,
  double a,
  double b
);

void gridding_2d (
	double span,
	int intervals,
	double * x,
	double * y,
	double force ( double x, double y ),
	double exact ( double x, double y ),
	double ** r,
	double * setup_time
);

void ftoc_2d (
	int node_fine,
	double kh,
	double ** b_fine,
	double ** u_fine,
	int node_coarse,
	double ** u_coarse,
	double ** residual_coarse,
	double * setup_time
);

void ctof_2d (
  int nodes_coarse,
  double ** u_coarse,
  double ** u_fine
);

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
		double * relax_time
	),
	char resolution,
	double step_size,
	double ** b_fine,
	double ** u_fine,
	int * multigrid_iterations,
	double * delta_l2,
	double * relax_time,
	double * setup_time
);