CFLAGS=-std=c99 -O3 -march=native -Wall -g -DNDEBUG
#CFLAGS=-std=c99 -O0 -march=native -Wall -g  # for debugging
LDLIBS=-lm

all: primes

primes: primes.c Makefile
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

clean:
	rm -f primes

.PHONY: all clean
