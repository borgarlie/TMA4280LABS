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
real solution(real x, real y);
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

    // OPEN MP TEST
    // #pragma omp parallel for num_threads(global_threads)
    // for (int i=0; i<10; i++) {
    // 	int tid = omp_get_thread_num();
    // 	int nthreads = omp_get_num_threads();
    // 	printf("Hello World from thread = %d\n", tid);
    // 	printf("Number of threads = %d\n", nthreads);
    // }

	/*
     * Grid points are generated with constant mesh size on both x- and y-axis.
     */
    real *grid = mk_1D_array(n+1, false);
    #pragma omp parallel for num_threads(global_threads)
    for (size_t i = 0; i < n+1; i++) {
        grid[i] = i * h;
    }

    /*
     * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
     * defined Chapter 9. page 93 of the Lecture Notes.
     * Note that the indexing starts from zero here, thus i+1.
     */
    real *diag = mk_1D_array(m, false);
    #pragma omp parallel for num_threads(global_threads)
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

    // int count = global_rank*local_n*global_matrix_size;

    // for (size_t i = 0; i < local_size; i++) {
    //     for (size_t j = 0; j < global_matrix_size; j++) {
    //         b[i][j] = count;
    //         bt[i][j] = count;
    //         count++;
    //     }
    // }

    // printf("\nPrinting bt before transpose (rank=%d) : \n", global_rank);

    // for (size_t i = 0; i < local_size; i++) {
    //     for (size_t j = 0; j < global_matrix_size; j++) {
    //     	printf(" %2.0f ", bt[i][j]);
    //     }
    //     printf(" (rank=%d)\n", global_rank);
    // }

    // MPI_Barrier(MPI_COMM_WORLD);

    // transpose(bt, b, local_n);

    // printf("\n\nPrinting bt after transpose (rank=%d) : \n", global_rank);

    // for (size_t i = 0; i < local_size; i++) {
    //     for (size_t j = 0; j < global_matrix_size; j++) {
    //     	printf(" %2.0f ", bt[i][j]);
    //     }
    //     printf(" (rank=%d)\n", global_rank);
    // }

    // DONE TESTING ::::::::::::::::::::::::

    /*
     * This vector will holds coefficients of the Discrete Sine Transform (DST)
     * but also of the Fast Fourier Transform used in the FORTRAN code.
     * The storage size is set to nn = 4 * n, look at Chapter 9. pages 98-100:
     * - Fourier coefficients are complex so storage is used for the real part
     *   and the imaginary part.
     * - Fourier coefficients are defined for j = [[ - (n-1), + (n-1) ]] while 
     *   DST coefficients are defined for j [[ 0, n-1 ]].
     * As explained in the Lecture notes coefficients for positive j are stored
     * first.
     * The array is allocated once and passed as arguments to avoid doings 
     * reallocations at each function call.
     */
    int nn = 4 * n;
    // Need to change z to a 2d array such that the different loops running
    // in parallel dont get race conditions - each thread has its own array
    // index with: omp_get_thread_num()
    real **z = mk_2D_array(global_threads, nn, false);

    /*
     * Initialize the right hand side data for a given rhs function.
     * Note that the right hand-side is set at nodes corresponding to degrees
     * of freedom, so it excludes the boundary (bug fixed by petterjf 2017).
     * 
     */
    #pragma omp parallel for collapse(2) num_threads(global_threads)
    for (size_t i = 0; i < local_size; i++) {
        for (size_t j = 0; j < global_matrix_size; j++) {
        	// must use offset here because of local_size
            b[i][j] = h * h * rhs(grid[i+1+offset], grid[j+1]);
        }
    }

    /*
     * Compute \tilde G^T = S^-1 * (S * G)^T (Chapter 9. page 101 step 1)
     * Instead of using two matrix-matrix products the Discrete Sine Transform
     * (DST) is used.
     * The DST code is implemented in FORTRAN in fsf.f and can be called from C.
     * The array zz is used as storage for DST coefficients and internally for 
     * FFT coefficients in fst_ and fstinv_.
     * In functions fst_ and fst_inv_ coefficients are written back to the input 
     * array (first argument) so that the initial values are overwritten.
     */
    #pragma omp parallel for num_threads(global_threads)
    for (size_t i = 0; i < local_size; i++) {
        fst_(b[i], &n, z[omp_get_thread_num()], &nn);
    }
    transpose(bt, b, local_n);
    #pragma omp parallel for num_threads(global_threads)
    for (size_t i = 0; i < local_size; i++) {
        fstinv_(bt[i], &n, z[omp_get_thread_num()], &nn);
    }

    /*
     * Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
     */
    #pragma omp parallel for collapse(2) num_threads(global_threads)
    for (size_t i = 0; i < local_size; i++) {
        for (size_t j = 0; j < global_matrix_size; j++) {
        	// Need to append offset to diag[i] since we use local size
            bt[i][j] = bt[i][j] / (diag[i + offset] + diag[j]);
        }
    }

    /*
     * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
     */
    #pragma omp parallel for num_threads(global_threads)
    for (size_t i = 0; i < local_size; i++) {
        fst_(bt[i], &n, z[omp_get_thread_num()], &nn);
    }
    transpose(b, bt, local_n);
    #pragma omp parallel for num_threads(global_threads)
    for (size_t i = 0; i < local_size; i++) {
        fstinv_(b[i], &n, z[omp_get_thread_num()], &nn);
    }

    /*
     * Compute maximal value of solution for convergence analysis in L_\infty
     * norm.
     */
    double u_max = 0.0;
    #pragma omp parallel for reduction(max: u_max) collapse(2) num_threads(global_threads)
    for (size_t i = 0; i < local_size; i++) {
        for (size_t j = 0; j < global_matrix_size; j++) {
            u_max = u_max > b[i][j] ? u_max : b[i][j];
        }
    }

    // Find the maximum error accross all ranks
    double u_max_reduced = 0.0;
    MPI_Reduce(&u_max, &u_max_reduced, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // printf("U_max : %e (Rank=%d) \n", u_max, global_rank);

    /*
    *	Compute the error as described in appendix B
    */
    double local_error = 0.0;
    for (size_t i = 0; i < local_size; i++) {
        for (size_t j = 0; j < global_matrix_size; j++) {
            real correct = solution(grid[(i + offset) + 1], grid[j+1]);
            real current_error = fabs(correct - b[i][j]);
            local_error = local_error > current_error ? local_error : current_error;
        }
    }

    double global_error = 0.0;
    MPI_Reduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // printf("Local error : %e (Rank=%d) \n", local_error, global_rank);
    if (rank == 0) {
    	printf("Global error: %e \n", global_error);
    }

	return u_max_reduced;
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

    #pragma omp parallel for num_threads(global_threads)
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
 * Solution used for verification as described in Appendix B
 */
real solution(real x, real y) {
    return sin(PI * x) * sin(2 * PI * y);
}

/*
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */
real rhs(real x, real y) {
    // return 2 * (y - y*y + x - x*x);
    return 5 * PI * PI * sin(PI * x) * sin(2 * PI * y);
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