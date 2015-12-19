#include "helmholtz.h"

bool trial (
	int intervals,
	char coarse_gamma,
	int parallel_size,
	int ddm_steps,
	int relax_steps,
	double force (
		double x,
		double y
	),
	double k,
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
	double exact (
		double x,
		double y
	)
) {

	// count # multigrid iterations needed
	int multigrid_iter = 0;
	double relax_time = 0.0;
	double setup_time = 0.0;

	double * x = spacing ( intervals + 1, -1.0, 1.0 );
	double ** u = ( double ** ) malloc ( ( intervals + 1 ) * sizeof ( double * ) );
	double ** r = ( double ** ) malloc ( ( intervals + 1 ) * sizeof ( double * ) );
	for ( int row=0 ; row<intervals+1 ; row++ ) {
		u[row] = ( double * ) malloc ( ( intervals + 1 ) * sizeof ( double ) );
		r[row] = ( double * ) malloc ( ( intervals + 1 ) * sizeof ( double ) );
		for ( int col=0 ; col<intervals+1 ; col++ ) u[row][col]=0.0;
	}

	gridding_2d (
		1.0-(-1.0),
		intervals,
		x,
		x,
		force,
		exact,
		r,
		&setup_time
	);

	double delta_l2;
	bool success = multigrid_2d (
		intervals,
		k,
		coarse_gamma,
		parallel_size,
		ddm_steps,
		relax_steps,
		relax,
		resolution,
		step_size,
		r,
		u,
		&multigrid_iter,
		&delta_l2,
		&relax_time,
		&setup_time
	);

	if (success&&false) {
		printf ( "\n" );
		printf ( "row\tcol\tx(row)        \ty(col)        \tu(x,y)        \texact(x,y)\n" );
		printf ( "\n" );
		for ( int row = 0; row < intervals + 1; row++ ) {
			for ( int col = 0; col < intervals + 1; col++ ) {
				printf ( "%d\t%d\t%e\t%e\t%e\t%e\n", row, col, x[row], x[col], u[row][col], exact ( x[row], x[col] ) );
			}
		}
	}

	if (success) printf ("%e\t%e\t%e\t%d\t%e\t%d\t%d\t%d\t%d\t%e\t%d\t%d\t%e\t%e",
		u[0][0],
		u[intervals/2][intervals/2],
		u[intervals][intervals],
		intervals,
		k,
		coarse_gamma,
		parallel_size,
		ddm_steps,
		relax_steps,
		step_size,
		resolution,
		multigrid_iter,
		// delta_l2,
		relax_time,
		setup_time
	);

	// printf ( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );
	// printf ( "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" );

	for (int row=0; row<intervals+1; row++) {
		free (r[row]);
		free (u[row]);
	}
	free (r);
	free (u);
	free (x);

	return (success);
}