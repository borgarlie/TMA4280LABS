#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double zeta0 ( int n );

int main ( int argc, char **argv ){

    if(argc != 2){
        printf("Usage: %s <n_iterations>\n", argv[0]);
        exit(-1);
    }

    int n = atoi(argv[1]);

    double result = zeta0(n);

    printf("Result: %f \n", result);
        
    exit ( EXIT_SUCCESS );
}

double zeta0 (int n) {
    return n;
}
