#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mach0.h"

double mach0 ( int n );
double arctan ( double x, int n);

int main ( int argc, char **argv ){

    if(argc != 2){
        printf("Usage: %s <n_iterations>\n", argv[0]);
        exit(-1);
    }

    int n = atoi(argv[1]);

    double result = mach0(n);

    printf("Result: %.17f \n", result);
        
    exit ( EXIT_SUCCESS );
}
