CPPFLAGS=-DNDEBUG
CFLAGS=-std=c99 -O3 -march=native -Wall -g
#CFLAGS=-std=c99 -O0 -march=native -Wall -g  # for debugging
LDLIBS=-lm

all: primes primes2 print-every

primes: primes.c Makefile
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

primes2: primes2.c Makefile
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

print-every: print-every.c Makefile
	$(CC) $(CFLAGS) $(CPPFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

clean:
	rm -f primes primes2 print-every

.PHONY: all clean
