#include "blas.h"

void transpose ( int n, double **A, double **AT ) {
	for (int row=0; row<n; row++)
		for (int col=0; col<n; col++)
			AT[row][col] = A[col][row];
}

double vec_mat_vec_mult ( int n, double *x, double **A ) {
	double sum = 0.0;
	for (int row=0; row<n; row++)
		for (int col=0; col<n; col++)
			sum += x[row] * A[row][col] * x[col];
	return sum;
}

double vec_mat_vec_mult_1d ( int n, double *x ) {
	double sum = 0.0;
	for (int row=0; row<n; row++) {
		sum += x[row] * (2.0) * x[row];
		if (0<row) sum += x[row] * (-1.0) * x[row-1];
		if (row<n-1) sum += x[row] * (-1.0) * x[row+1];
		// for (int col=0; col<n; col++) {
		// 	sum += x[row] * A[row][col] * x[col];
		// }
	}
	return sum;
}

double vec_mat_vec_mult_2d ( int unused, int lin_size, double *x ) {
	double sum = 0.0;
	for (int grid_row=0; grid_row<lin_size; grid_row++) {
		for (int grid_col=0; grid_col<lin_size; grid_col++) {
			int index = grid_row*lin_size + grid_col;
			sum += x[index] * (4.0) * x[index];
			if (0<grid_row) sum += x[index] * (-1.0) * x[index-lin_size];
			if (grid_row<lin_size-1) sum += x[index] * (-1.0) * x[index+lin_size];
			if (0<grid_col) sum += x[index] * (-1.0) * x[index-1];
			if (grid_col<lin_size-1) sum += x[index] * (-1.0) * x[index+1];
		}
	}
	return sum;
}

double vec_mat_vec_mult_3d ( int unused, int lin_size, double *x ) {
	double sum = 0.0;
	for (int grid_row=0; grid_row<lin_size; grid_row++) {
		for (int grid_col=0; grid_col<lin_size; grid_col++) {
			for (int grid_sli=0; grid_sli<lin_size; grid_sli++) {
				int index = grid_row*lin_size*lin_size + grid_col*lin_size + grid_sli;
				sum += x[index] * (8.0) * x[index];
				if (0<grid_row) sum += x[index] * (-1.0) * x[index-lin_size*lin_size];
				if (grid_row<lin_size-1) sum += x[index] * (-1.0) * x[index+lin_size*lin_size];
				if (0<grid_col) sum += x[index] * (-1.0) * x[index-lin_size];
				if (grid_col<lin_size-1) sum += x[index] * (-1.0) * x[index+lin_size];
				if (0<grid_sli) sum += x[index] * (-1.0) * x[index-1];
				if (grid_sli<lin_size-1) sum += x[index] * (-1.0) * x[index+1];
			}
		}
	}
	return sum;
}

long vec_mat_vec_mult_fixed ( int n, long *x, long **A ) {
	long sum = 0;
	for (int row=0; row<n; row++) {
		// printf("vec_mat_vec_mult_fixed x[%d]=%ld\n",row,x[row]);
		for (int col=0; col<n; col++) {
			// printf("vec_mat_vec_mult_fixed A[%d][%d]=%ld\n",row,col,A[row][col]);
			sum += x[row] * A[row][col] * x[col];
		}
	}
	// printf("vec_mat_vec_mult_fixed sum=%ld\n",sum);
	return sum;
}

void mat_vec_mult ( int n, double **A, double *x, double *result ) {
	for (int row=0; row<n; row++) {
		result[row] = 0.0;
		for (int col=0; col<n; col++)
			result[row] += A[row][col] * x[col];
	}
}

double vec_vec_mult ( int n, double *a, double *b ) {
	double sum = 0.0;
	for (int row=0; row<n; row++)
		sum += a[row] * b[row];
	return sum;
}

long vec_vec_mult_fixed ( int n, long *a, long *b ) {
	long sum = 0;
	for (int row=0; row<n; row++)
		sum += a[row] * b[row];
	return sum;
}

void mat_mat_mult ( int n, double **A, double **B, double **result ) {
	for (int row=0; row<n; row++)
		for (int col=0; col<n; col++) {
			result[row][col] = 0.0;
			for (int iter=0; iter<n; iter++)
				result[row][col] += A[row][iter] * B[iter][col];
		}
}

void inverse ( int n, double ** inp, double ** out ) {
	switch (n) {
		case 2: inverse_2x2 (inp, out); break;
		case 3: inverse_3x3 (inp, out); break;
		default: break;
	}
}

void inverse_2x2 ( double ** inp, double ** out ) {
	double ad = (inp[0][0])*(inp[1][1]);
	double bc = (inp[0][1])*(inp[1][0]);
	double det = ad-bc;
	double invdet = 1 / ( det );
	out[0][0] = (inp[1][1])*(invdet);
	out[0][1] = -(inp[0][1])*(invdet);
	out[1][0] = -(inp[1][0])*(invdet);
	out[1][1] = (inp[0][0])*(invdet);
}

void inverse_3x3 (double ** a, double ** out) {
	int i, j;
	double det = 0.0;
	for (i=0;i<3;i++) {
		det = det + (
			a[0][i] * (
				a[1][(i+1)%3] * a[2][(i+2)%3] -
				a[1][(i+2)%3] * a[2][(i+1)%3]
			)
		);
	}
	for (i=0;i<3;i++) {
		for (j=0;j<3;j++)
			out [j][i] = (
				a[(i+1)%3][(j+1)%3] * a[(i+2)%3][(j+2)%3] - 
				a[(i+1)%3][(j+2)%3] * a[(i+2)%3][(j+1)%3]
			) / det;
	}
}

double two_norm ( int n, double *a, double *b ) {
	double sum = 0.0;
	for ( int row=0; row<n; row++ )
		sum += (a[row]-b[row]) * (a[row]-b[row]);
	return sqrt(sum);
}

double energy_norm ( int n, double ** A, double *x, double *y ) {
	double sum = 0.0;
	for ( int row=0; row<n; row++ )
		for ( int col=0; col<n; col++ )
			sum += (x[row]-y[row]) * A[row][col] * (x[col]-y[col]);
	return sqrt(fabs(sum));
}

void print_vec ( int n, double *b ) {
	if (n>0) printf("[%f", b[0]);
	for (int row=1; row<n; row++) printf(", %f", b[row]);
	printf("]");
}

void print_mat ( int n, int m, double **A ) {
	printf("[");
	if (n>0) print_vec(m, A[0]);
	for (int row=1; row<n; row++) {
		printf(",\n");
		print_vec(m, A[row]);
	}
	printf("]");
}

void print_sq_mat ( int n, double **A ) {
	printf("[");
	if (n>0) print_vec(n, A[0]);
	for (int row=1; row<n; row++) {
		printf(",\n");
		print_vec(n, A[row]);
	}
	printf("]");
}

double to_float (char resolution, long fixed_num) {
	return ((double)fixed_num) / (double)(1<<resolution);
}

long to_fixed (char resolution, double float_num) {
	// printf("float_num = %f\n", float_num);
	// printf("to_fixed 1<<resolution = %d\n", 1<<resolution);
	// printf("float_num * (double)(1<<resolution) = %f\n", float_num * (double)(1<<resolution));
	// printf("round(float_num * (double)(1<<resolution)) = %f\n", round(float_num * (double)(1<<resolution)));
	// printf("(long) round(float_num * (double)(1<<resolution)) = %ld\n", (long) round(float_num * (double)(1<<resolution)));
	return round(float_num * (double)(1<<resolution));
}