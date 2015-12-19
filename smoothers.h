#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <time.h>

#include "blas.h"

int sor (
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

// sor_fixed infeasible because sor not efficient in dda

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
);

// gauss_seidel infeasible because gauss_seidel not efficient in dda

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
);

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
);

// stochastic unreasonable because matrix is already sparse

int smoother_float (
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

int smoother_fixed (
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