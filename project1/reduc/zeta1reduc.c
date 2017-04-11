#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../common/common.h"
#include <mpi.h>
#include "zeta1reduc.h"

int main(int argc, char **argv){
    
    if(argc != 2){
        printf("Usage: %s <n_iterations>\n", argv[0]);
        exit(-1);
    }

    // initialise MPI
    int size, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n = atoi(argv[1]);

    assert(isPowTwo(size) == 1);

    double result = zeta1reduc(n, size, rank);
    printf("Result at rank %d: %.50f \n", rank, result);

    // Finalize the MPI environment.
    MPI_Finalize();

    exit(EXIT_SUCCESS);
}
