CFLAGS+=-std=c99 -O3
LDLIBS=-lm
TARGETS=zeta0 testzeta0 vtestzeta0

all: ${TARGETS}

clean:
	-rm -f ${TARGETS}

run:
	./zeta0 100

utest:
	./testzeta0

vtest:
	./vtestzeta0
