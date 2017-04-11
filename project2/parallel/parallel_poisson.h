#ifndef PARALLEL_POISSON_H
#define PARALLEL_POISSON_H

typedef double real;
typedef int bool;

// Function prototypes
double parallel_poisson(int n_threads_omp, int iterations, int size, int rank);
void commit_vector_types(int local_n, int remainder);
real *mk_1D_array(size_t n, bool zero);
real **mk_2D_array(size_t n1, size_t n2, bool zero);
void transpose(real **bt, real **b, int local_n);
real rhs(real x, real y);

// Functions implemented in FORTRAN in fst.f and called from C.
// The trailing underscore comes from a convention for symbol names, called name
// mangling: if can differ with compilers.
void fst_(real *v, int *n, real *w, int *nn);
void fstinv_(real *v, int *n, real *w, int *nn);

// MPI datatypes used for transpose operation and ??
MPI_Datatype
	row, col, 
	row_remainder, col_remainder,
	block_send, block_receive, 
	block_send_remainder_right, block_receive_remainder_right,
	block_send_remainder_bottom, block_receive_remainder_bottom,
	corner_send, corner_receieve;

// global variables
int global_size, global_rank, global_matrix_size, global_threads;

/*
	The main function
*/
double parallel_poisson(int n_threads_omp, int n, int size, int rank) {
	global_size = size;
	global_rank = rank;
	global_threads = n_threads_omp;

	int m = n - 1;
    real h = 1.0 / n;

    global_matrix_size = m;

    int local_n = n / size;
    int remainder = n - 1 - local_n * (size-1);
    int offset = rank * local_n;

	commit_vector_types(local_n, remainder);

	int local_size = local_n;
	if (rank == size-1) {
    	local_size = remainder;
    }

    // TODO: Check if all of the stuff here is still used?

	/*
     * Grid points are generated with constant mesh size on both x- and y-axis.
     */
    real *grid = mk_1D_array(n+1, false);
    for (size_t i = 0; i < n+1; i++) {
        grid[i] = i * h;
    }

    /*
     * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
     * defined Chapter 9. page 93 of the Lecture Notes.
     * Note that the indexing starts from zero here, thus i+1.
     */
    real *diag = mk_1D_array(m, false);
    for (size_t i = 0; i < m; i++) {
        diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
    }

    /*
     * Allocate the matrices b and bt which will be used for storing value of
     * G, \tilde G^T, \tilde U^T, U as described in Chapter 9. page 101.
     */
    real **b = mk_2D_array(local_size, m, false);
    real **bt = mk_2D_array(local_size, m, false);

    // TESTING:::::::::::::::::::::::::::::

    int count = global_rank*local_n*global_matrix_size;

    for (size_t i = 0; i < local_size; i++) {
        for (size_t j = 0; j < global_matrix_size; j++) {
            b[i][j] = count;
            bt[i][j] = count;
            count++;
        }
    }

    printf("\nPrinting bt before transpose (rank=%d) : \n", global_rank);

    for (size_t i = 0; i < local_size; i++) {
        for (size_t j = 0; j < global_matrix_size; j++) {
        	printf(" %2.0f ", bt[i][j]);
        }
        printf(" (rank=%d)\n", global_rank);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    transpose(bt, b, local_n);

    printf("\n\nPrinting bt after transpose (rank=%d) : \n", global_rank);

    for (size_t i = 0; i < local_size; i++) {
        for (size_t j = 0; j < global_matrix_size; j++) {
        	printf(" %2.0f ", bt[i][j]);
        }
        printf(" (rank=%d)\n", global_rank);
    }

    // DONE TESTING ::::::::::::::::::::::::

	double u_max = 2;
	return u_max;
}

/*
	Used to create the vector types used when sending data between processors
*/
void commit_vector_types(int local_n, int remainder) {

	// printf("Global size: %d, Global rank: %d\n", global_size, global_rank);

	int column_stride = global_matrix_size;

	// Datatypes for row and column for whole blocks

	MPI_Datatype temp_row, temp_col;

	MPI_Type_vector(local_n, 1, 1, MPI_DOUBLE, &temp_row);
	MPI_Type_create_resized(temp_row, 0, sizeof(double), &row);
	MPI_Type_commit(&row);

    MPI_Type_vector(local_n, 1, column_stride, MPI_DOUBLE, &temp_col);
    MPI_Type_create_resized(temp_col, 0, sizeof(double), &col);
    MPI_Type_commit(&col);

    // Datatypes for whole blocks

	MPI_Datatype 
		temp_block_send, temp_block_receive, 
		temp_block_send_remainder_right, temp_block_receive_remainder_right,
		temp_block_send_remainder_bottom, temp_block_receive_remainder_bottom;

	MPI_Type_vector(local_n, 1, column_stride, row, &temp_block_send);
	MPI_Type_create_resized(temp_block_send, 0, sizeof(double), &block_send);
	MPI_Type_commit(&block_send);

	MPI_Type_vector(local_n, 1, 1, col, &temp_block_receive);
	MPI_Type_create_resized(temp_block_receive, 0, sizeof(double), &block_receive);
	MPI_Type_commit(&block_receive);

	// Datatypes for row and column for remainder

	MPI_Datatype temp_row_remainder, temp_col_remainder;

	MPI_Type_vector(remainder, 1, 1, MPI_DOUBLE, &temp_row_remainder);
	MPI_Type_create_resized(temp_row_remainder, 0, sizeof(double), &row_remainder);
	MPI_Type_commit(&row_remainder);

    MPI_Type_vector(remainder, 1, column_stride, MPI_DOUBLE, &temp_col_remainder);
    MPI_Type_create_resized(temp_col_remainder, 0, sizeof(double), &col_remainder);
    MPI_Type_commit(&col_remainder);

    // Datatypes for whole remainder blocks

	MPI_Type_vector(local_n, 1, column_stride, row_remainder, &temp_block_send_remainder_right);
	MPI_Type_create_resized(temp_block_send_remainder_right, 0, sizeof(double), &block_send_remainder_right);
	MPI_Type_commit(&block_send_remainder_right);

	MPI_Type_vector(local_n, 1, 1, col_remainder, &temp_block_receive_remainder_right);
	MPI_Type_create_resized(temp_block_receive_remainder_right, 0, sizeof(double), &block_receive_remainder_right);
	MPI_Type_commit(&block_receive_remainder_right);

	MPI_Type_vector(remainder, 1, column_stride, row, &temp_block_send_remainder_bottom);
	MPI_Type_create_resized(temp_block_send_remainder_bottom, 0, sizeof(double), &block_send_remainder_bottom);
	MPI_Type_commit(&block_send_remainder_bottom);

	MPI_Type_vector(remainder, 1, 1, col, &temp_block_receive_remainder_bottom);
	MPI_Type_create_resized(temp_block_receive_remainder_bottom, 0, sizeof(double), &block_receive_remainder_bottom);
	MPI_Type_commit(&block_receive_remainder_bottom);

	// Datatypes for last corner where there is remainder in both directions

	MPI_Datatype temp_corner_send, temp_corner_receive;

	MPI_Type_vector(remainder, 1, column_stride, row_remainder, &temp_corner_send);
	MPI_Type_create_resized(temp_corner_send, 0, sizeof(double), &corner_send);
	MPI_Type_commit(&corner_send);

	MPI_Type_vector(remainder, 1, 1, col_remainder, &temp_corner_receive);
	MPI_Type_create_resized(temp_corner_receive, 0, sizeof(double), &corner_receieve);
	MPI_Type_commit(&corner_receieve);

}

void transpose(real **bt, real **b, int local_n) {
	
	int *send_counts = malloc(global_size * sizeof(int));
	int *send_displacements = malloc(global_size * sizeof(int));
	MPI_Datatype *send_datatypes = malloc(global_size * sizeof(MPI_Datatype));

    int *receive_counts = malloc(global_size * sizeof(int));
    int *receive_displacements= malloc(global_size * sizeof(int));
    MPI_Datatype *receive_datatypes = malloc(global_size * sizeof(MPI_Datatype));

    for (size_t i = 0; i < global_size; i++) {
    	send_counts[i] = 1;
    	receive_counts[i] = 1;
    	send_displacements[i] = (int) (i * local_n * sizeof(double));
    	receive_displacements[i] = (int) (i * local_n * sizeof(double));

    	if (global_rank == global_size - 1 && i == global_size - 1) {
    		send_datatypes[i] = corner_send;
    		receive_datatypes[i] = corner_receieve;
    		// maybe just swap them instead of sending ? for all diagonals i==rank
    		// send_counts[i] = 0;
    		// receive_counts[i] = 0;
    	}
    	else if (global_rank == global_size - 1) {
    		send_datatypes[i] = block_receive_remainder_right;
    		receive_datatypes[i] = block_send_remainder_bottom;
    		// send_counts[i] = 0;
    		// receive_counts[i] = 0;
    	}
    	else if (i == global_size - 1) {
    		send_datatypes[i] = block_receive_remainder_bottom;
    		receive_datatypes[i] = block_send_remainder_right;
    		// send_counts[i] = 0;
    		// receive_counts[i] = 0;
    	} else {
    		// should otherwise send and receieve whole blocks
    		send_datatypes[i] = block_send;
    		receive_datatypes[i] = block_receive;
    	}

    }

    MPI_Alltoallw(
    	b[0], send_counts, send_displacements, send_datatypes,
        bt[0], receive_counts, receive_displacements, receive_datatypes, 
        MPI_COMM_WORLD);

    free(send_counts);
    free(send_displacements);
    free(send_datatypes);

    free(receive_counts);
    free(receive_displacements);
    free(receive_datatypes);
}

/*
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */

real rhs(real x, real y) {
    return 2 * (y - y*y + x - x*x);
}

 /*
 * The allocation of a vectore of size n is done with just allocating an array.
 * The only thing to notice here is the use of calloc to zero the array.
 */

real *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (real *)calloc(n, sizeof(real));
    }
    return (real *)malloc(n * sizeof(real));
}

/*
 * The allocation of the two-dimensional array used for storing matrices is done
 * in the following way for a matrix in R^(n1*n2):
 * 1. an array of pointers is allocated, one pointer for each row,
 * 2. a 'flat' array of size n1*n2 is allocated to ensure that the memory space
 *   is contigusous,
 * 3. pointers are set for each row to the address of first element.
 */

real **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    // 1
    real **ret = (real **)malloc(n1 * sizeof(real *));

    // 2
    if (zero) {
        ret[0] = (real *)calloc(n1 * n2, sizeof(real));
    }
    else {
        ret[0] = (real *)malloc(n1 * n2 * sizeof(real));
    }
    
    // 3
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}

#endif