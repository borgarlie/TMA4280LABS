CC=mpicc
CFLAGS+=-std=c99 -O3
LDLIBS=-lm
TARGETS=zeta1
NP=8

all: ${TARGETS}

clean:
	-rm -f ${TARGETS}

run:
	mpirun -np ${NP} zeta1 1024
