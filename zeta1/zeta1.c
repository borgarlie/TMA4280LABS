#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

double zeta1 (int n);
int isPowTwo(int n);

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

int isPowTwo(int x) {
	while (((x % 2) == 0) && x > 1) /* While x is even and > 1 */
		x /= 2;
 	return (x == 1);
}

double zeta1 (int n) {

    double sum = 0;

    for (double i=1; i<=n; i++) {
        sum += 1/(i*i);
    }

    double result = sqrt(sum*6);

    return result;
}