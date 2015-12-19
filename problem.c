#include "helmholtz.h"

void problem (
	double force (
		double x,
		double y
	),
	double k,
	double exact (
		double x,
		double y
	)
) {


	printf("u(0)        \tu(n/2)        \tu(n)        \tintrvls\tk        \tgamma\tprlll\tddm\trelax\tstep_size\trsltn\tmg itrs\trelax seconds\tsetup seconds\tsolver\n");
	// intervals 1 - 65536
	for (int interval_pow=0 ; interval_pow<9 ; interval_pow++) {
		int intervals = 1<<interval_pow;

		// choice of how many W cycles to undergo
		for (int coarse_gamma=1 ; coarse_gamma<2 ; coarse_gamma++ ){

			// ddm jacobi / gauss-seidel / sor size 1 - 65536
			for (int parallel_pow=0 ; parallel_pow<interval_pow+1; parallel_pow++ ) {
				int parallel_size = 1<<parallel_pow;

				// ddm jacobi / gauss-seidel / sor steps 1 - 65536
				for (int ddm_pow=0 ; ddm_pow<1; ddm_pow++ ) {
					int ddm_steps = 1<<ddm_pow;

	// printf("u(0)        \tu(n/2)        \tu(n)        \tintrvls\tk        \tgamma\tprlll\tddm\trelax\tstep_size\trsltn\tmg itrs\trelax seconds\tsetup seconds\tsolver\n");

					// linear solver allowed to run 1 - 65536 steps
					for (int relax_pow=8 ; relax_pow<9; relax_pow++ ) {
						int relax_steps = 1<<(2*relax_pow);

							if (trial ( intervals, coarse_gamma, parallel_size, ddm_steps, relax_steps, force, k, conjugate_solver, 0, 1.0, exact ) ) {
								printf ("\tconjugate\n");
							}
					}
				}
			}
		}
	}
}