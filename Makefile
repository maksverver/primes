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

# Runs some basic tests. Pretty slow. They are intended to catch regressions.
# To speed up the `primes` tests it's possible to run with cache enabled, e.g.
# PRIMES_CACHE=/tmp/primes.bin make check
check: check-primes check-primes2

check-primes: primes
	test "$$( ./primes -n 10000 | sha256sum | cut -c1-64 )" = 'de1b90e91ee8193f153cd9d6f79887a1ba05e2365a8ee230c6c6eb23c1ab5fe4'
	test "$$( ./primes | head -n 400000000 | sha256sum | cut -c1-64 )" = '8a539551865283276fb871cfd6c937bfc9d7b0b7f577d9400b39a43960111e5f'
	test "$$( ./primes 18446744073000000000 | sha256sum | cut -c1-64 )" = 'b541fb774c03b1fd6bc7944b4507865e79ba63d7358f6038756874b0c2a28734'

check-primes2: primes2
	test "$$( ./primes2 | head -n 400000000 | sha256sum | cut -c1-64 )" = '8a539551865283276fb871cfd6c937bfc9d7b0b7f577d9400b39a43960111e5f'

.PHONY: all clean check check-primes check-primes2
