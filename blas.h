#include <stdio.h>
#include <math.h>

void transpose (
	int n,
	double **A,
	double **AT
);

double vec_mat_vec_mult (
	int n,
	double *x,
	double **A
);

double vec_mat_vec_mult_1d (
	int n,
	double *x
);

double vec_mat_vec_mult_2d (
	int n,
	int lin_size,
	double *x
);

double vec_mat_vec_mult_3d (
	int n,
	int lin_size,
	double *x
);

long vec_mat_vec_mult_fixed (
	int n,
	long *x,
	long **A
);

void mat_vec_mult (
	int n,
	double **A,
	double *x,
	double *result
);

double vec_vec_mult (
	int n,
	double *a,
	double *b
);

long vec_vec_mult_fixed (
	int n,
	long *a,
	long *b
);

void mat_mat_mult (
	int n,
	double **A,
	double **B,
	double **result
);

void inverse (
	int n,
	double ** inp,
	double ** out
);

void inverse_2x2 (
	double ** inp,
	double ** out
);

void inverse_3x3 (
	double ** inp,
	double ** out
);

double two_norm (
	int n,
	double *a,
	double *b
);

double energy_norm (
	int n,
	double ** A,
	double *x,
	double *y
);

void print_vec (
	int n,
	double *b
);

void print_mat (
	int n,
	int m,
	double **A
);

void print_sq_mat (
	int n,
	double **A
);

double to_float (
	char resolution,
	long fixed_num
);

long to_fixed (
	char resolution,
	double float_num
);