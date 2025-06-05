/* Generates a long list of prime numbers.

Unlike primes.c, this program doesn't have any configurable options. It only
prints all the prime numbers, in decimal notation, one per line, in order, and
the code is optimized for being fast and concise. It's around 33% faster than
primes.c, though it does use about 4x as much memory (2.5 GiB vs 630 MiB).

Since this prints a lot of output, it's probably useful to filter it somehow,
for example:

% ./primes2 | ./print-every -n 1000000

...will print every 1 million-th prime number.

In theory, this program should print every prime number below 2**64 eventually
(without overflowing!) but since the total number of those primes is about 400
million billion (!) this will take too long in practice.

Conceptually, the algorithm divides the output range into 2**32 pages of 2**32
numbers each. It uses the Sieve of Eratosthenes to find the primes in the first
page (exactly 203,280,221 primes below 2**32), and uses those to sieve later
pages. It's not necessary to keep track of larger primes, because any number
below 2**64 that's not a prime has a divisor below 2**32.

An optimization of the segmented sieve algorithm (compared to primes.c) is that
this version keeps track of the last multiple of each prime used, which means
there is much less overhead per prime.

Also unlike primes.c this program does not cache the primes below 2**32 on disk,
so generating the initial small primes list may take a while (around 30 seconds).

The code here also includes some of the other tricks I came up with while
implementing primes.c, such as the optimization of printing consecutive decimal
numbers. */

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define PRIME_COUNT_UNDER_20 8
static const uint32_t primes_under_20[PRIME_COUNT_UNDER_20] = { 2, 3, 5, 7, 11, 13, 17, 19 };

#define PRIME_COUNT_32_BITS 203280221
#define PRIME_COUNT_SMALL (PRIME_COUNT_32_BITS - PRIME_COUNT_UNDER_20)
static uint32_t small_primes[PRIME_COUNT_SMALL];  /* ~775 MiB */
static uint64_t next_multiple[PRIME_COUNT_SMALL]; /* ~1.55 GiB */

static uint8_t init_pattern[3*5*7*11*13*17*19];   /* ~4.6 MiB */
static uint8_t odd_bitmap[UINT32_MAX / 16 + 1];   /*  256 MiB */

/* The following code is used to optimize printing of consecutive integers.

The average gap between primes is relatively small (approximately log(p_i)) so
copying the last number and adding the difference is significantly faster than
generating the decimal representation from scratch. */

/* pow10min1[i] = 10**i - 1, i.e., the maximum integer that fits into i digits. */
static const uint64_t pow10min1[21] = {
    0, 9, 99, 999, 9999, 99999, 999999, 9999999, 99999999, 999999999, 9999999999,
    99999999999, 999999999999, 9999999999999, 99999999999999, 999999999999999,
    9999999999999999, 99999999999999999, 999999999999999999, 9999999999999999999ULL,
    UINT64_MAX};
   
static uint8_t out_buf[32768];  /* 32 KiB, half of a pipe buffer */
static size_t out_pos;          /* write position in the buffer */
static size_t out_last_len;     /* length of last number written including '\n' */
static uint64_t out_last_val;   /* last number written */
static uint64_t out_next_max;   /* maximum number that fits in `out_last_len`. */

static void flush() {
    if (write(1, out_buf, out_pos) != out_pos) {
        perror("write() failed\n");
        exit(1);
    }
    out_pos = 0;
    out_last_val = 0;
    out_last_len = 0;
    out_next_max = 0;
}

static void output_prime(uint64_t p) {
    assert(p != 0 && p >= out_last_val);
    if (p <= out_next_max && out_last_len <= sizeof(out_buf) - out_pos) {
        uint64_t add = p - out_last_val;
        memcpy(&out_buf[out_pos], &out_buf[out_pos - out_last_len], out_last_len);
        out_pos += out_last_len;
        for (char *p = (char*)&out_buf[out_pos - 2]; add > 0; --p) {
            add += *p - '0';
            *p = add%10 + '0';
            add /= 10;
        }
        out_last_val = p;
        return;
    }

    if (out_pos >= sizeof(out_buf) - 100) flush();  /* ensure buffer space */
    out_last_val = p;
    out_last_len = snprintf(
            (char*)&out_buf[out_pos], sizeof(out_buf) - out_pos,
            "%" PRIu64 "\n", p);
    assert(1 <= out_last_len && out_last_len <= 21);
    out_pos += out_last_len;
    out_next_max = pow10min1[out_last_len - 1];
}

/* Next follows the code to calculate primes. */

/* Initializes odd_bitmap with a periodic pattern with multiples of the primes
under 20 already marked off, which speeds up the first few iterations of the
sieve, which are the slowest. */
static void initialize_odd_bitmap(uint64_t start) {
    memset(init_pattern, 0, sizeof(init_pattern));
    for (size_t i = 0; i < sizeof(init_pattern); ++i) {
        uint64_t p = start + 2*i + 1;
        if (p%3 && p%5 && p%7 && p%11 && p%13 && p%17 && p%19) {
            for (size_t k = 0; k < 8; ++k) {
                size_t j = i + k*sizeof(init_pattern);
                init_pattern[j >> 3] |= 1 << (j & 7);
            }
        }
    }
    for (size_t i = 0; ; ) {
        if (sizeof(init_pattern) <= sizeof(odd_bitmap) - i) {
            memcpy(odd_bitmap + i, init_pattern, sizeof(init_pattern));
            i += sizeof(init_pattern);
        } else {
            memcpy(odd_bitmap + i, init_pattern, sizeof(odd_bitmap) - i);
            break;
        }
    }
}

/* Generates and outputs all primes between 20 and UINT32_MAX, populating the
   small_primes[] and next_multiple[] tables along the way. */
static void generate_small_primes() {
    fprintf(stderr, "Calculating %d primes up to %u.\n", PRIME_COUNT_32_BITS, UINT32_MAX);
    initialize_odd_bitmap(0);
    size_t k = 0;
    for (size_t i = 1; i <= UINT32_MAX/2; ++i) {
        if (odd_bitmap[i >> 3] & (1 << (i & 7))) {
            size_t j = 2*i*i + 2*i;
            while (j <= UINT32_MAX/2) {
                odd_bitmap[j >> 3] &= ~(1 << (j & 7));
                j += 2*i + 1;
            }
            small_primes[k] = 2*i + 1;
            next_multiple[k] = j;  /* note odd index */
            ++k;
            output_prime(2*i + 1);
            if ((k & (k - 1)) == 0) flush();  /* flush output early at the start */
        }
    }
    if (k != PRIME_COUNT_SMALL) {  /* sanity check */
        fprintf(stderr, "Incorrect number of small primes!\n"
            "Expected: %zu\nFound: %zu\n", k, (size_t)PRIME_COUNT_SMALL);
        exit(1);
    }
    flush();
}

int main() {
    /* The smallest primes (under 20) are handled separately, so print them first. */
    for (size_t i = 0; i < PRIME_COUNT_UNDER_20; ++i) {
        output_prime(primes_under_20[i]);
    }
    flush();

    /* Generate 32-bit primes which are used to sieve larger numbers. */
    generate_small_primes();

    /* Sieve the larger numbers, one 32-bit page at a time. */
    size_t in_range = 0;
    for (uint64_t page = 1; page <= UINT32_MAX; ++page) {
        initialize_odd_bitmap(page << 32);
        uint64_t min_index = (page << 31);
        uint64_t max_index = (page << 31) + UINT32_MAX/2;
        while (in_range < PRIME_COUNT_SMALL && next_multiple[in_range] < max_index) {
            assert(next_multiple[in_range] >= min_index);
            next_multiple[in_range] -= min_index;
            ++in_range;
        }

        fprintf(stderr,
            "Sieving page %"PRIu64", from %"PRIu64" to %"PRIu64", using %zd primes.\n",
            page, (page << 32) + 1, (page << 32) + UINT32_MAX, in_range);

        /* This inner loop is doing most of the work. */
        for (size_t i = 0; i < in_range; ++i) {
            uint64_t m = next_multiple[i];
            uint64_t p = small_primes[i];
            do {
                odd_bitmap[m >> 3] &= ~(1 << (m & 7));
                m += p;
            } while (m <= UINT32_MAX/2);
            next_multiple[i] = m - (UINT32_MAX/2 + 1);
        }

        /* Output primes found on this page. */
        for (size_t j = 0; j <= UINT32_MAX/2; ++j) {
            if (odd_bitmap[j >> 3] & (1 << (j & 7))) {
                output_prime((min_index + j)*2 + 1);
            }
        }
        flush();
    }
}
