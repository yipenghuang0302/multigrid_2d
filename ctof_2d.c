void ctof_2d (
	int nodes_coarse,
	double ** u_coarse,
	double ** u_fine
) {

	int row_coarse, row_fine;
	int col_coarse, col_fine;

	// copy row
	for (row_coarse = 0; row_coarse < nodes_coarse; row_coarse++) {
		row_fine = 2 * row_coarse;
		// copy col
		for (col_coarse = 0; col_coarse < nodes_coarse; col_coarse++) {
			col_fine = 2 * col_coarse;
			u_fine[row_fine][col_fine] += u_coarse[row_coarse][col_coarse];
		}

		// interpolate col
		for (col_coarse = 0; col_coarse < nodes_coarse - 1; col_coarse++) {
			col_fine = 2 * col_coarse + 1;
			u_fine[row_fine][col_fine] += 0.5 * ( u_coarse[row_coarse][col_coarse] + u_coarse[row_coarse][col_coarse+1] );
		}
	}

	// interpolate row
	for (row_coarse = 0; row_coarse < nodes_coarse - 1; row_coarse++) {
		row_fine = 2 * row_coarse + 1;
		// copy col
		for (col_coarse = 0; col_coarse < nodes_coarse; col_coarse++) {
			col_fine = 2 * col_coarse;
			u_fine[row_fine][col_fine] += 0.5 * ( u_coarse[row_coarse][col_coarse] + u_coarse[row_coarse+1][col_coarse] );
		}

		// interpolate col
		for (col_coarse = 0; col_coarse < nodes_coarse - 1; col_coarse++) {
			col_fine = 2 * col_coarse + 1;
			u_fine[row_fine][col_fine] += 0.25 * ( u_coarse[row_coarse][col_coarse] + u_coarse[row_coarse][col_coarse+1] + u_coarse[row_coarse+1][col_coarse] + u_coarse[row_coarse+1][col_coarse+1] );
		}
	}

	return;
}