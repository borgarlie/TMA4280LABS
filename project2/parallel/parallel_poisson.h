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
	block_send, block_receive, 
	block_send_remainder_right, block_receive_remainder_right,
	block_send_remainder_bottom, block_receive_remainder_bottom;

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

    // if (rank == size-1) {
    // 	local_n = remainder;
    // }

	commit_vector_types(local_n, remainder);

	double u_max = 2;
	return u_max;
}

/*
	Assumes that there is a squere matrix that always is dividable by number of processors
*/
void commit_vector_types(int local_n, int remainder) {

	// printf("Global size: %d, Global rank: %d\n", global_size, global_rank);

	MPI_Datatype row, col;

	MPI_Datatype 
		temp_block_send, temp_block_receive, 
		temp_block_send_remainder_right, temp_block_receive_remainder_right,
		temp_block_send_remainder_bottom, temp_block_receive_remainder_bottom;

	int column_stride = global_matrix_size;

	MPI_Type_vector(local_n, local_n, column_stride, MPI_FLOAT, &temp_block_send);
	MPI_Type_create_resized(temp_block_send, 0, sizeof(float), &block_send);
	MPI_Type_commit(&block_send);

	MPI_Type_vector(local_n, local_n, column_stride, MPI_FLOAT, &temp_block_receive);
	MPI_Type_create_resized(temp_block_receive, 0, sizeof(float), &block_receive);
	MPI_Type_commit(&block_receive);

}


void transpose(real **bt, real **b, int local_n) {

}

/*
 * Write the transpose of b a matrix of R^(m*m) in bt.
 * In parallel the function MPI_Alltoallv is used to map directly the entries
 * stored in the array to the block structure, using displacement arrays.
 */

// void transpose(real **bt, real **b, size_t m)
// {
//     for (size_t i = 0; i < m; i++) {
//         for (size_t j = 0; j < m; j++) {
//             bt[i][j] = b[j][i];
//         }
//     }
// }

#endif