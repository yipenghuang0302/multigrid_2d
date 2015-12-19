#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <time.h>

#include "blas.h"

int conjugate_solver_3d (
	int n,
	char resolution,
	int max_step,
	double step_size,
	double * b,
	double * x_next,
	double * delta_l2,
	double * time
);

int conjugate_solver_2d (
	int n,
	char resolution,
	int max_step,
	double step_size,
	double * b,
	double * x_next,
	double * delta_l2,
	double * time
);

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
);

int steepest_solver_float (
	int n,
	char resolution,
	int max_step,
	double step_size,
	double ** A,
	double * b,
	double * x_next,
	double * delta_l2,
	double * time
);

int steepest_solver_fixed (
	int n,
	char resolution,
	int max_step,
	double step_size,
	double ** A,
	double * b,
	double * x_next,
	double * delta_l2,
	double * time
);