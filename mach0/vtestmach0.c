#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mach0.h"

void vtestmach0 ();

int main ( int argc, char **argv ){
	vtestmach0();
    exit ( EXIT_SUCCESS );
}

void vtestmach0() {

	FILE *f = fopen("verification.txt", "w");

	for (int k=1; k<=24; k++) {
		int n = pow(2, k);
		double result = mach0(n);
		double diff = fabs(result - M_PI);
		fprintf(f, "Diff between pi and mach0(n=%d): %.50f\n", n, diff);
	}

	fclose(f);
}
