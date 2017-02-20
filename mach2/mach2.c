#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mach2.h"
#include <omp.h>

int main ( int argc, char **argv ){

    if(argc != 3){
        printf("Usage: %s <n_threads> <n_iterations>\n", argv[0]);
        exit(-1);
    }

    int nthreads = atoi(argv[1]);
    int n = atoi(argv[2]);

    double result = mach2(nthreads, n);

    printf("Result: %.50f \n", result);
        
    exit ( EXIT_SUCCESS );
}
