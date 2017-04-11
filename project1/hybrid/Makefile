CC=mpicc
CFLAGS+=-std=c99 -O3
LDLIBS=-lm -fopenmp
TARGETS=zeta1hybrid
NP=8
OMPP=4

all: ${TARGETS}

clean:
	-rm -f ${TARGETS}

run:
	mpirun -np ${NP} zeta1 {OMPP} 1024
