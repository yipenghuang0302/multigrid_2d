#include "norm.h"

#include "blas.h"

int norm (
	int n,
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
	double ** A,
	double * b,
	double * x_next
) {

	int row;
	double delta_l2, **AT, *ATb, **ATA;

	AT = (double**) malloc (n*sizeof(double*));
	ATb = (double*) malloc (n*sizeof(double));
	ATA = (double**) malloc (n*sizeof(double*));
	
	// allocate
	for (row=0; row<n; row++) {
		AT[row] = (double*) malloc(n*sizeof(double));
		ATA[row] = (double*) malloc(n*sizeof(double));
	}

	// aux variables
	transpose (n, A, AT);
	mat_mat_mult (n, AT, A, ATA);
	mat_vec_mult (n, AT, b, ATb);

	double time;
	int step = relax ( n, 16, 1<<16, 0.125, ATA, ATb, x_next, &delta_l2, &time );
	printf("norm took %.5f seconds\n", time);

	for (row=0; row<n; row++) {
		free(AT[row]);
		free(ATA[row]);
	}
	free (AT);
	free (ATA);
	free (ATb);

	return step;

}