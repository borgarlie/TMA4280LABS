#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../common/common.h"
#include <mpi.h>
#include "zeta1hybrid.h"
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

    // assert that the number of processes are correct
    // COMMENT: I think this is an ugly way to handle this.
    assert(isPowTwo(size) == 1);

    double result = zeta1_hybrid(nthreads, n, size, rank);

    if (rank == 0) {
        printf("Result: %.50f \n", result);
    }

    // Finalize the MPI environment.
    MPI_Finalize();

    exit(EXIT_SUCCESS);
}
