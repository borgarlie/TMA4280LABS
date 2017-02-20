#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../common/common.h"
#include <mpi.h>
#include "mach1.h"
#include <unistd.h>


void vtestmach1 ();

int main ( int argc, char **argv ){
	vtestmach1();
    exit ( EXIT_SUCCESS );
}

void vtestmach1() {
	// initialise MPI
    int size, rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // assert that the number of processes are correct
    // COMMENT: I think this is an ugly way to handle this.
    assert(isPowTwo(size) == 1);

	FILE *f;
	if (rank == 0) {
		char filename[100];
		sprintf(filename, "stats/verification%d.txt", size); // puts string into buffer
        f = fopen(filename, "w");
    }

	for (int k=16; k<=24; k++) {
		int n = pow(2, k);
		double t1 = MPI_Wtime();
		double result = mpi_mach1(n, size, rank);
		double t2 = MPI_Wtime();
    	double time = t2-t1;
		if (rank == 0) {
			double diff = fabs(result - M_PI);
			fprintf(f, "Diff between pi and mach1(n=%d): %.50f. Time: %.10f\n", n, diff, time);
		}
	}

	fclose(f);

	// Finalize the MPI environment.
    MPI_Finalize();
}
