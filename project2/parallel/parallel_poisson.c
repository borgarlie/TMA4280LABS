#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../common/common.h"
#include <mpi.h>
#include "parallel_poisson.h"
#include <omp.h>

int main(int argc, char **argv){
    
    if(argc != 3){
        printf("Usage: %s <n_threads_omp> <n_iterations>\n", argv[0]);
        exit(-1);
    }

    // initialise MPI
    int size, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int nthreads = atoi(argv[1]);
    int n = atoi(argv[2]);

    // check that n is a power of two
    // assert(isPowTwo(n) == 1);

    double start_time;
    if (rank == 0) {
        start_time = MPI_Wtime();
    }

    // this is now u_max_reduced
    double result = parallel_poisson(nthreads, n, size, rank);

    if (rank == 0) {
        double end_time = MPI_Wtime();
        double total_time = end_time - start_time;
        printf("Result: %e \n", result);
        printf("Total time: %f sec \n", total_time);
    }

    // Finalize the MPI environment.
    MPI_Finalize();

    exit(EXIT_SUCCESS);
}
