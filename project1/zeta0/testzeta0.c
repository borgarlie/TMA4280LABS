#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "zeta0.h"

void testzeta0 ();

int main ( int argc, char **argv ){
	testzeta0();
    exit ( EXIT_SUCCESS );
}

void testzeta0() {
	double result = zeta0(3);
	double expected = 2.85773;

	double epsilon = 0.001;
	double diff = fabs(result - expected);

	assert(diff <= epsilon);
}
