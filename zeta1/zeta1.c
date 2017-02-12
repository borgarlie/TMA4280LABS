#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../common/common.h"

double zeta1 (int n);

int main ( int argc, char **argv ){

	if(argc != 2){
        printf("Usage: %s <n_iterations>\n", argv[0]);
        exit(-1);
    }

    int n = atoi(argv[1]);

    assert(isPowTwo(n) == 1);

    double result = zeta1(n);

    printf("Result: %.17f \n", result);
        
    exit ( EXIT_SUCCESS );
}

double zeta1 (int n) {

    double sum = 0;

    for (double i=1; i<=n; i++) {
        sum += 1/(i*i);
    }

    double result = sqrt(sum*6);

    return result;
}