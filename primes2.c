/* Generates a long list of prime numbers.

Unlike primes.c, this program doesn't have any configurable options. It only
prints all the prime numbers, in decimal notation, one per line, in order, and
the code is optimized for being fast and concise. It's around 33% faster than
primes.c, though it does use around 3x as much memory (1.8 GiB vs 630 MiB).

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
this version keeps track of the next multiple of each prime (offset within the
page only, to save space).

Also unlike primes.c this program does not cache the primes below 2**32 on disk,
so generating the initial small primes list may take a while (around 20 seconds).

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

#define PRIME_COUNT_16_BITS      6542  /* number of 16-bit primes:       6,542 */
#define PRIME_COUNT_32_BITS 203280221  /* number of 32-bit primes: 203,280,221 */
#define PRIME_COUNT_TINY    (PRIME_COUNT_16_BITS - PRIME_COUNT_UNDER_20)
#define PRIME_COUNT_SMALL   (PRIME_COUNT_32_BITS - PRIME_COUNT_UNDER_20)
static uint32_t small_primes[PRIME_COUNT_SMALL];  /* ~775 MiB */
static uint32_t next_multiple[PRIME_COUNT_SMALL]; /* ~775 MiB */

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
    if (sizeof(out_buf) - out_pos < 24) {
        flush();  /* ensure buffer space */
    } else if (p <= out_next_max) {
        /* copy previous line */
        char tmp[24];
        memcpy(tmp, &out_buf[out_pos - out_last_len], 24);
        memcpy(&out_buf[out_pos], tmp, 24);
        out_pos += out_last_len;

        /* add difference to last value */
        uint64_t add = p - out_last_val;
        out_last_val = p;
        for (char *s = (char*)&out_buf[out_pos - 2]; add > 0; --s) {
            add += *s - '0';
            *s = add%10 + '0';
            add /= 10;
        }
        return;
    }

    out_last_val = p;
    out_last_len = snprintf(
            (char*)&out_buf[out_pos], sizeof(out_buf) - out_pos,
            "%" PRIu64 "\n", p);
    assert(1 <= out_last_len && out_last_len <= 21);
    out_pos += out_last_len;
    out_next_max = pow10min1[out_last_len - 1];
}

/* Next follows the code to calculate primes. */

static uint64_t square(uint32_t p) {
    return (uint64_t)p * p;
}

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
    initialize_odd_bitmap(0);
    size_t k = 0;
    for (size_t i = 1; i <= UINT16_MAX/2; ++i) {
        if (odd_bitmap[i >> 3] & (1 << (i & 7))) {
            size_t j = 2*i*i + 2*i;
            while (j <= UINT16_MAX/2) {
                odd_bitmap[j >> 3] &= ~(1 << (j & 7));
                j += 2*i + 1;
            }
            small_primes[k] = 2*i + 1;
            next_multiple[k] = j;  /* note odd index */
            ++k;
            output_prime(2*i + 1);
        }
    }
    flush();
    size_t n = k;  /* number of primes between 20 and 2**16 */
    assert(n == PRIME_COUNT_TINY);
    for (uint32_t page = 1; page <= UINT16_MAX; ++page) {
        for (size_t i = 0; i < n; ++i) {
            uint32_t m = next_multiple[i];
            uint32_t p = small_primes[i];
            while (m <= (page << 15) + UINT16_MAX/2) {
                odd_bitmap[(page << 12) | (m >> 3)] &= ~(1 << (m & 7));
                m += p;
            }
            next_multiple[i] = m;
        }
        /* Output primes found on this page. */
        for (size_t i = (page << 15); i <= (page << 15) + UINT16_MAX/2; ++i) {
            if (odd_bitmap[i >> 3] & (1 << (i & 7))) {
                small_primes[k] = 2*i + 1;
                next_multiple[k] = 2*i*i + 2*i;  /* note odd index; possible overflow */
                ++k;
                output_prime(2*i + 1);
            }
        }
        flush();
    }
    if (k != PRIME_COUNT_SMALL) {  /* sanity check */
        fprintf(stderr, "Incorrect number of small primes!\n"
            "Expected: %zu\nFound: %zu\n", k, (size_t)PRIME_COUNT_SMALL);
        exit(1);
    }
    for (size_t i = 0; i < k; ++i) next_multiple[i] &= UINT32_MAX/2;
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
        while ( in_range < PRIME_COUNT_SMALL &&
                (square(small_primes[in_range]) >> 32) <= page ) ++in_range;

        uint64_t start = (uint64_t) page << 32;
        initialize_odd_bitmap(start);

        fprintf(stderr,
            "Sieving page %"PRIu64", from %"PRIu64" to %"PRIu64", using %zd primes.\n",
            page, start + 1, start + UINT32_MAX, in_range);

        /* This inner loop is doing most of the work. */
        for (size_t i = 0; i < in_range; ++i) {
            uint32_t m = next_multiple[i];
            uint32_t p = small_primes[i];
            do {
                odd_bitmap[m >> 3] &= ~(1 << (m & 7));
                m += p;
            } while (m <= UINT32_MAX/2);
            next_multiple[i] = m & (UINT32_MAX/2);
        }

        /* Output primes found on this page. */
        for (size_t j = 0; j <= UINT32_MAX/2; ++j) {
            if (odd_bitmap[j >> 3] & (1 << (j & 7))) {
                output_prime(start + 2*j + 1);
            }
        }
        flush();
    }
}
