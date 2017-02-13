#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../common/common.h"
#include <mpi.h>

// function prototypes
void compute();
void gather();
void print_values();

// variables
int size; // number of processes
int rank; // rank of process
double *values; // global array containing the final values
double *local_values; // local array containing the final values
int n; // number of total iterations
int local_count; // numer of local iterations

int main(int argc, char **argv){

    // get number of iterations
	if(argc != 2){
        printf("Usage: %s <n_iterations>\n", argv[0]);
        exit(-1);
    }
    n = atoi(argv[1]);

    // initialise MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // find number of local iterations
    // TODO: should any n be handled and evenly distributed?
    local_count = n / size;
    printf("Local count: %d \n", local_count);

    // assert that the number of processes are correct
    // COMMENT: I think this is an ugly way to handle this.
    assert(isPowTwo(size) == 1);

    // setup arrays
    if(rank == 0){
        values = calloc(n, sizeof(double));
    }
    local_values = calloc(local_count, sizeof(double));

    compute();
    // TODO: check if it really is needed with a barrier
    MPI_Barrier(MPI_COMM_WORLD);
    gather();

    // print values on rank 0 and free memory
    if(rank == 0){
        print_values();
        free (values);
    }
    // free all local allocated memory
    free(local_values);

    // Finalize the MPI environment.
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}

void compute() {
    // shifting everything by 1 to include n, but not include 0
    int start = (local_count * rank) + 1;
    int end = start + local_count;
    int count = 0;
    // Compute values and put them in local values array
    for (double i=start; i<end; i++) {
        local_values[count++] = 1/(i*i);
    }
}

void gather() {
    MPI_Gather(local_values, local_count, MPI_DOUBLE, values, local_count, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
} 

void print_values() {
    double sum = 0;
    for (int i=0; i <n; i++) {
        sum += values[i];
    }
    double result = sqrt(sum*6);
    printf("Result: %.17f \n", result);
}
