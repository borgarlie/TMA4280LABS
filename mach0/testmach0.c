#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "mach0.h"

void testmach0 ();

int main ( int argc, char **argv ){
	testmach0();
    exit ( EXIT_SUCCESS );
}

void testmach0() {
	double result = mach0(3);
	double expected = 3.141592;

	double epsilon = 0.001;
	double diff = fabs(result - expected);

	assert(diff <= epsilon);
}
