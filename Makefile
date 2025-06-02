CFLAGS=-std=c99 -O3 -march=native -Wall -g
LDLIBS=-lm

all: primes

clean:
	rm -f primes
