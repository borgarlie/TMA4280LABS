CC=mpicc
CFLAGS+=-std=c99 -O3
LDLIBS=-lm
TARGETS=zeta1reduc
NP=8

all: ${TARGETS}

clean:
	-rm -f ${TARGETS}

run:
	mpirun -np ${NP} zeta1reduc 1024
