CC=mpicc
CFLAGS+=-std=c99 -O3
LDLIBS=-lm
TARGETS=mach1 vtestmach1
NP=8

all: ${TARGETS}

clean:
	-rm -f ${TARGETS}

run:
	mpirun -np ${NP} mach1 1024

vtest:
	for number in 1 2 4 8 16; do \
		mpirun -np $$number vtestmach1 ; \
	done