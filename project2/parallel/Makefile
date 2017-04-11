CC=mpicc
CFLAGS+=-std=c99 -O3
LDLIBS=-lm -fopenmp
TARGETS=parallel_poisson
NP=8
OMPP=4

all: ${TARGETS}

clean:
	-rm -f ${TARGETS}

run:
	mpirun -np ${NP} parallel_poisson ${OMPP} 128
