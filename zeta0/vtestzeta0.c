#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "zeta0.h"

void vtestzeta0 ();

int main ( int argc, char **argv ){
	vtestzeta0();
    exit ( EXIT_SUCCESS );
}

void vtestzeta0() {

	FILE *f = fopen("verification.txt", "w");

	for (int k=1; k<=24; k++) {
		int n = pow(2, k);
		double result = zeta0(n);
		double diff = fabs(result - M_PI);
		fprintf(f, "Diff between pi and zeta0(n=%d): %.50f\n", n, diff);
	}

	fclose(f);
}
