/* Program to print prime numbers below 2**64.

Uses a segmented Sieve of Eratosthenes, where the first 2**32 primes are
precalculated and stored on disk compressed to 120 MB using a prime wheel that
stores the primarily of every 2×3×5×7 = 210 integers in 48 bits = 6 bytes.

Compilation (example):

% cc -std=c99 -O3 -march=native -Wall -g primes.c  -lm -o primes

% ./primes -h

will print help on how to run.


TODO's:

    - maybe: option to check primes using trial division? this should be faster
        than sieving when max - min or count is small.

    - maybe: option to output composites instead of primes?

    - maybe:
        - store wheel data in shared memory backed by shm_open()?
        - mmap() prime wheel if possible?

    - maybe: add a checksum to the wheel data to verify it is intact?
      this would be a bit tricky when loading partial data. (though I could
      hardcode checksums of ranges in the source code: checksum of first N primes
      is Y etc.)

    - clean up/improve estimate_max_ub() logic?

    - support printing of interesting number with a minimum
      (removes the ability to mark primes interesting based on index)

    - maybe: add an option to use the cache readonly (load if it exists, but do
      not overwrite).

Not planned:

    - Option to print the last N primes instead of the first N.
      This can be achieved with `tail -N` instead (though it is a bit slower).

    - Option to skip the first N primes. This can be emulated with tail, for
      example to skip 90 primes and generate the next 10, you can do:

        ./primes -n 100 | tail -n +91
*/

#include <assert.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Assume we're on a little-endian architecture. (This only matters for
// load_wheel_bitmap() and save_wheel_bitmap() at this point. On a big-endian
// architecture these files would simply use a big-endian representation
// instead.)
#ifndef htole64
#define htole64(x) (x)
#endif
#ifndef le64toh
#define le64toh(x) (x)
#endif

const uint64_t wheel_magic_id = 3472612198888075856;  // "Prime210"

// Prime wheel of factors {2, 3, 5, 7} with size 2×3×5×7 = 210.
//
// 48 of every 210 integers are indivisible by 2, 3, 5 and 7.
// The table below lists those integers modulo 210.
//
// This is used to encode whether integers are prime using only 48 bits
// (6 bytes) per 210 numbers.
static const uint8_t wheel210_offset[48] = {
//   0    1    2    3    4    5    6    7    8    9
     1,  11,  13,  17,  19,  23,  29,  31,  37,  41,  // 10
    43,  47,  53,  59,  61,  67,  71,  73,  79,  83,  // 20
    89,  97, 101, 103, 107, 109, 113, 121, 127, 131,  // 30
   137, 139, 143, 149, 151, 157, 163, 167, 169, 173,  // 40
   179, 181, 187, 191, 193, 197, 199, 209
};

#ifndef NDEBUG
// Inverse of wheel210_offset. wheel210_index[i] = -1 if i is divisible by
// 2, 3, 5, or 7, otherwise it lists an index into wheel210_offset such that
// wheel210_offset[wheel_index[i]] = i.
static const int8_t wheel210_index[210] = {
//  0   1   2   3   4   5   6   7   8   9
   -1,  0, -1, -1, -1, -1, -1, -1, -1, -1,   //  10
   -1,  1, -1,  2, -1, -1, -1,  3, -1,  4,   //  20
   -1, -1, -1,  5, -1, -1, -1, -1, -1,  6,   //  30
   -1,  7, -1, -1, -1, -1, -1,  8, -1, -1,   //  40
   -1,  9, -1, 10, -1, -1, -1, 11, -1, -1,   //  50
   -1, -1, -1, 12, -1, -1, -1, -1, -1, 13,   //  60
   -1, 14, -1, -1, -1, -1, -1, 15, -1, -1,   //  70
   -1, 16, -1, 17, -1, -1, -1, -1, -1, 18,   //  80
   -1, -1, -1, 19, -1, -1, -1, -1, -1, 20,   //  90
   -1, -1, -1, -1, -1, -1, -1, 21, -1, -1,   // 100
   -1, 22, -1, 23, -1, -1, -1, 24, -1, 25,   // 110
   -1, -1, -1, 26, -1, -1, -1, -1, -1, -1,   // 120
   -1, 27, -1, -1, -1, -1, -1, 28, -1, -1,   // 130
   -1, 29, -1, -1, -1, -1, -1, 30, -1, 31,   // 140
   -1, -1, -1, 32, -1, -1, -1, -1, -1, 33,   // 150
   -1, 34, -1, -1, -1, -1, -1, 35, -1, -1,   // 160
   -1, -1, -1, 36, -1, -1, -1, 37, -1, 38,   // 170
   -1, -1, -1, 39, -1, -1, -1, -1, -1, 40,   // 180
   -1, 41, -1, -1, -1, -1, -1, 42, -1, -1,   // 190
   -1, 43, -1, 44, -1, -1, -1, 45, -1, 46,   // 200
   -1, -1, -1, -1, -1, -1, -1, -1, -1, 47,   // 210
};
#endif

// Same as wheel210_index, but only elements at odd indices (since all elements
// at even indices are -1).
static const int8_t wheel210_index_half[105] = {
//  1   3   5   7   9   11  13  15  17  19
    0, -1, -1, -1, -1,   1,  2, -1,  3,  4,   //  20
   -1,  5, -1, -1,  6,   7, -1, -1,  8, -1,   //  40
    9, 10, -1, 11, -1,  -1, 12, -1, -1, 13,   //  60
   14, -1, -1, 15, -1,  16, 17, -1, -1, 18,   //  80
   -1, 19, -1, -1, 20,  -1, -1, -1, 21, -1,   // 100
   22, 23, -1, 24, 25,  -1, 26, -1, -1, -1,   // 120
   27, -1, -1, 28, -1,  29, -1, -1, 30, 31,   // 140
   -1, 32, -1, -1, 33,  34, -1, -1, 35, -1,   // 160
   -1, 36, -1, 37, 38,  -1, 39, -1, -1, 40,   // 180
   41, -1, -1, 42, -1,  43, 44, -1, 45, 46,   // 200
   -1, -1, -1, -1, 47,
};

// Macros to calculate the minimum/maximum of two arguments.
//
// Note that these evaluate their arguments more than once, so should only
// be used with expressions that are cheap to compute and have no side-effects!
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

// Size in bytes of an array storing a bitmap of odd integers up to `max`.
#define ODD_BITMAP_SIZE(max)   ((max) / 16 + 1)

// Size in bytes of an array storing a bitmap of integers up to `max` where
// only 48/210 integers are stored (see the definition of wheel-210 above).
#define WHEEL_BITMAP_SIZE(max) (((max) / 210 + 1)*6)

#ifndef NDEBUG
static void static_assertions() {
    for (int i = 0; i < 48; ++i) {
        assert(wheel210_index[(int) wheel210_offset[i]] == i);
    }
    for (int i = 0; i < 210; ++i) {
        int j = wheel210_index[i];
        assert(j == -1 ||wheel210_offset[j] == i);
    }
    for (int i = 0; i < 105; ++i) {
        assert(wheel210_index[2*i + 0] == -1);
        assert(wheel210_index[2*i + 1] == wheel210_index_half[i]);
    }
}
#endif

/* Clears bit as indices `start`, `start + step`, `start + 2*step`, etc. until
   `start + i*step >= len`. */
static void clear_bits(uint8_t *bitmap, uint64_t start, uint64_t step, uint64_t len) {
    for (uint64_t i = start; i < len; i += step) {
        bitmap[i >> 3] &= ~(1 << (i & 7));
    }
}

/* Finds prime numbers no greater than `max` using the Sieve of Eratosthenes.

Returns the result as a bitmap of odd integers, where the i-th bit is set
iff. 2i + 1 is prime.

The result may have more bits set than expected because this actually rounds
`max` up to a multiple of 16 to simplify the implementation. */
static uint8_t *sieve_up_to_max(uint64_t max) {
    if (ODD_BITMAP_SIZE(max) > SIZE_MAX) {
        errno = ENOMEM;
        return NULL;
    }

    size_t bitmap_size = ODD_BITMAP_SIZE(max);
    uint8_t *bitmap = malloc(bitmap_size);
    if (bitmap == NULL) return NULL;  /* failed to allocate */
    memset(bitmap, 0xff, bitmap_size);
    bitmap[0] = 0xfe; /* 1 is not prime */

    uint64_t len = (uint64_t) bitmap_size * 8;  /* bitmap size in bits */
    for (uint64_t i = 1; ; ++i) {
        if (bitmap[i >> 3] & (1 << (i & 7))) {
            uint64_t j = 2*i*i + 2*i;
            if (j > len) return bitmap;
            clear_bits(bitmap, j, 2*i + 1, len);
        }
    }
}


/* Calculates a prime wheel bitmap.

The wheel bitmap uses 48 bits per 210 numbers (22.857% compression) which is
less than half of the odd numbers bitmap with 1 bit per 2 numbers (50%).

The wheel bitmap is calculated by first sieving into an odd bitmap, then
converting the odd bitmap to the wheel bitmap.
*/
static uint8_t *generate_prime_wheel(uint64_t max) {
    /* Check only the odd-bitmap size, which is greater than or equal to the
       wheel bitmap size. */
    if (ODD_BITMAP_SIZE(max) > SIZE_MAX) {
        errno = ENOMEM;
        return NULL;
    }

    /* Allocate memory in advance, so we don't waste our time sieving if there
       isn't enough memory for the wheel bitmap. */
    size_t wheel_size = WHEEL_BITMAP_SIZE(max);
    uint8_t *wheel_bitmap = malloc(wheel_size);
    if (wheel_bitmap == NULL) return NULL;

    /* Sieve to generate the odd bitmap. */
    uint8_t *odd_bitmap = sieve_up_to_max(max);
    size_t odd_size = ODD_BITMAP_SIZE(max);
    if (odd_bitmap == NULL) return NULL;

    /* Now we'll convert from the odd bitmap to the wheel bitmap. */
    memset(wheel_bitmap, 0, wheel_size);
    for (uint64_t i = 4; (i >> 3) < odd_size; ++i) {
        if (odd_bitmap[i >> 3] & (1 << (i & 7))) {
            assert(wheel210_index_half[i % 105] >= 0);
            uint64_t j = (i / 105) * 48 + wheel210_index_half[i % 105];
            assert((j >> 3) < wheel_size);
            wheel_bitmap[j >> 3] |= 1 << (j & 7);
        }
    }
    return wheel_bitmap;
}

/* Saves the wheel bitmap to a file.

The file format is very simple:

     - 8 bytes: `max` in little endian format
     - WHEEL_BITMAP_SIZE(max) bytes: the wheel bitmap
*/
static int save_wheel_bitmap(
        const char *cache_filename,
        const uint8_t *wheel,
        uint64_t max)
{
    FILE *fp = fopen(cache_filename, "wb");
    if (fp == NULL) return -1;
    int res = 0;
    uint64_t stored_max   = htole64(max);
    uint64_t stored_magic = htole64(wheel_magic_id);
    size_t byte_size = WHEEL_BITMAP_SIZE(max);
    if ( fwrite(&stored_magic, sizeof(stored_magic), 1, fp) != 1 ||
         fwrite(&stored_max,   sizeof(stored_max),   1, fp) != 1 ||
         fwrite(wheel, 1, byte_size, fp) != byte_size ) {
        /* Write failed. */
        errno = EIO;
        res = -2;
    }
    if (fclose(fp) != 0 && res == 0) res = -3;
    if (res != 0 && remove(cache_filename) != 0) {
        perror("failed to remove incomplete file");
    }
    return res;
}

static uint8_t *load_wheel_bitmap(const char *filename, uint64_t max) {
    if (WHEEL_BITMAP_SIZE(max) > SIZE_MAX) {
        errno = ENOMEM;
        return NULL;
    }
    size_t size = WHEEL_BITMAP_SIZE(max);
    uint8_t *data = NULL;  /* allocation comes later */

    FILE *fp = fopen(filename, "rb");
    if (fp == NULL) return NULL;

    /* Check the stored maximum. We can use the existing data if it's larger,
       but not if it's smaller. (Technically we could use the smaller
       bitmap to generate a larger one, but the complexity isn't worth it.) */
    uint64_t stored_magic = 0;
    uint64_t stored_max = 0;
    if (fread(&stored_magic, sizeof(stored_magic), 1, fp) != 1) goto eof;
    if (le64toh(stored_magic) != wheel_magic_id) {
        errno = EINVAL;
        fprintf(stderr, "cached prime wheel bitmap has invalid signature\n");
        goto fail;
    }

    if (fread(&stored_max, sizeof(stored_max), 1, fp) != 1) goto eof;
    if (le64toh(stored_max) < max) {
        fprintf(stderr, "note: not using cached prime wheel bitmap because it is too small\n");
        errno = ENOENT;  /* suppress scary warnings at the call site */
        goto fail;
    }

    /* Data should be good! Read the prefix we need. */
    if ((data = malloc(size)) == NULL) goto fail;
    if (fread(data, 1, size, fp) != size) goto eof;
    goto done;

eof:
    errno = EIO;
fail:
    free(data);
    data = NULL;
done:
    (void) fclose(fp);  /* don't care about failure */
    return data;
}


/* Callback function that receives a context pointer and the next prime.
Should return 0 to keep enumeration, or any other value to stop, which will
be returned from enumerate_primes(). */
typedef int (*enumerate_prime_callback_t)(void *context, uint64_t prime);

static int enumerate_small_primes_odd(
        uint8_t *odd_bitmap,
        uint32_t min,
        uint32_t max,
        enumerate_prime_callback_t callback,
        void *context)
{
    int res = 0;

    /* Handle 2 first. */
    if (max < 2 || (min <= 2 && (res = callback(context, 2)) != 0)) return res;

    /* Handle odd integers. Here max >= 2 so (max - 1) won't underflow. */
    for (uint32_t i = min / 2; i <= (max - 1)/2; ++i) {
        if (odd_bitmap[i >> 3] & (1 << (i & 7))) {
            if ((res = callback(context, 2*i + 1)) != 0) return res;
        }
    }

    return res;
}

static int enumerate_small_primes_sieved(
        uint32_t min,
        uint32_t max,
        enumerate_prime_callback_t callback,
        void *context)
{
    uint8_t *bitmap = sieve_up_to_max(max);
    if (bitmap == NULL) {
        perror("couldn't generate odd prime bitmap");
        return -1;
    }
    int res = enumerate_small_primes_odd(bitmap, min, max, callback, context);
    free(bitmap);
    return res;
}

static int enumerate_small_primes_wheel(
        const uint8_t *wheel_bitmap,
        uint32_t min,
        uint32_t max,
        enumerate_prime_callback_t callback,
        void *context)
{
    int res = 0;
    if (max < 2 || (min <= 2 && (res = callback(context, 2)) != 0)) return res;
    if (max < 3 || (min <= 3 && (res = callback(context, 3)) != 0)) return res;
    if (max < 5 || (min <= 5 && (res = callback(context, 5)) != 0)) return res;
    if (max < 7 || (min <= 7 && (res = callback(context, 7)) != 0)) return res;

    size_t i = min / 210;
    /* First chunk of 210 integers may contain primes p < min if min % 210 != 0 */
    if (min % 210 > 0) {
        for (size_t j = 0; j < 48; ++j) {
            if (wheel_bitmap[6*i + (j >> 3)] & (1 << (j & 7))) {
                uint64_t p = (uint64_t) 210*i + wheel210_offset[j];
                if (p < min) continue;
                if ((res = callback(context, p)) != 0) return res;
            }
        }
        ++i;
    }
    /* Middle chunks of integers: all primes are between min and max. */
    while (i < max / 210) {
        for (size_t j = 0; j < 48; ++j) {
            if (wheel_bitmap[6*i + (j >> 3)] & (1 << (j & 7))) {
                uint64_t p = (uint64_t) 210*i + wheel210_offset[j];
                if ((res = callback(context, p)) != 0) return res;
            }
        }
        ++i;
    }
    /* Final chunk of 210 integers may contain primes p > max */
    if (max % 210 > 0) {
        for (size_t j = 0; j < 48; ++j) {
            if (wheel_bitmap[6*i + (j >> 3)] & (1 << (j & 7))) {
                uint64_t p = (uint64_t) 210*i + wheel210_offset[j];
                if (p > max) break;
                if ((res = callback(context, p)) != 0) return res;
            }
        }
    }
    return 0;
}

struct sieve_segment_context {
    uint8_t *bitmap;
    uint64_t len, min, max;
};

/* Removes the multiples of `p` from the bitmap of length `len` (in bits) that
represents the odd integers starting from `min`.

Assumes min <= max, p <= UINT32_MAX, and `min`, `max` and `p` are all odd. */
int sieve_segment(void *ctx_arg, uint64_t p) {
    struct sieve_segment_context *ctx = ctx_arg;
    uint64_t start = p * p;
    if (start > ctx->max) return 1;
    if (start < ctx->min) {
        /* Round `start` up to the first multiple of p at least `min`,
           but be careful to avoid 64-bit overflow: */
        start = (ctx->min - 1) / p * p;
        if (start % 2 == 0) {
            if (ctx->max - start < p) return 0; /* skip this one */
            start += p;
        } else {
            if (ctx->max - start < 2*p) return 0; /* skip this one */
            start += 2*p;
        }
    }
    assert(ctx->min % 2 == 1);
    assert(start % 2 == 1);
    assert(start >= ctx->min);
    assert(start <= ctx->max);
    clear_bits(ctx->bitmap, (start - ctx->min)/2, p, ctx->len);
    return 0;
}

/* Enumerates primes in the range from min to max (inclusive).

`wheel` must be a pointer to the wheel bitmap of all primes below UINT32_MAX
(2**32). For each prime found, callback(context) is called. If the callback
returns 0, enumeration continues.

If `supports_flush` is true, then this function will invoke callback with
argument 0 whenever a large batch of output is complete. This allows e.g.
flushing output buffered so far. Note that this is purely informational;
there is no guarantee this will be called at all!

This function returns the nonzero value returned by callback(), or 0 if the
callback always returned 0 (or was never called at all).

This function may allocate up to 256 MiB for additional sieving. If allocation
fails, this function returns -1 and errno is set to ENOMEM.
*/
int enumerate_primes_sieved(
        const uint8_t *wheel,
        uint64_t min, uint64_t max,
        enumerate_prime_callback_t callback, void *context,
        bool supports_flush)
{
    /* Normalize min and max so they are both odd and min <= max.
       This handles some edge cases and is required by sieve_segment below. */
    if (min % 2 == 0) ++min;
    if (max % 2 == 0 && max > 0) --max;
    if (min > max) return 0;  /* empty range */

    if (max <= UINT32_MAX) {
        /* Range is small. Enumerate from wheel bitmap and return. */
        return enumerate_small_primes_wheel(wheel, min, max, callback, context);
    }

    /* Allocate memory in advance, so we fail early if there isn't enough memory
       to finish the operation. */
    uint8_t *bitmap = malloc(MIN(max - min, UINT32_MAX) / 16 + 1);
    if (bitmap == NULL) return -1;

    /* Handle the primes up to UINT32_MAX separately: */
    int res = 0;
    if (min <= UINT32_MAX) {
        res = enumerate_small_primes_wheel(wheel, min, UINT32_MAX, callback, context);
        if (res != 0) goto finish;
        min = (uint64_t) UINT32_MAX + 2;  /* +2 to keep min odd! */
    }

    /* Now min >= 2**32 and we will sieve the remaining numbers one segment of
       up to 2**32 integers at a time. Since there are 203,280,221 primes below
       UINT32_MAX, this can get relatively slow for later ranges. */
    for (;;) {
        if (supports_flush && (res = callback(context, 0)) != 0) goto finish;

        assert(min % 2 == 1 && max % 2 == 1);

        uint64_t diff = MIN(max - min, UINT32_MAX);
        memset(bitmap, 0xff, diff/16 + 1);

        struct sieve_segment_context ssc = {
            .bitmap = bitmap,
            .len = diff/2 + 1,
            .min = min,
            .max = min + diff,
        };
        enumerate_small_primes_wheel(wheel, 3, UINT32_MAX, sieve_segment, &ssc);

        for (uint64_t i = 0; i <= diff/2; ++i) {
            if (bitmap[i >> 3] & (1 << (i & 7))) {
                uint64_t p = min + 2*i;  /* note: `min` is odd here */
                res = callback(context, p);
            }
        }

        if (min + diff >= max) break;  /* avoid overflow */
        min += diff + 1;
    }

finish:
    free(bitmap);
    return res;
}

enum interest {
    INT_NONE = 0,
    INT_IDX  = 1,
    INT_DEC  = 2,
    INT_HEX  = 4,
};

struct print_context {
    FILE *fp;
    uint64_t count;
    uint64_t max_prime;
    uint64_t max_count;

    union {
        /* Used to optimize printing of consecutive integers.

        When buf_len > 0, buf stores a string of length len of the form
        "<buf_val>\n" followed by a terminating NUL character. Note that 64-bit
        numbers have most 20 digits, so we need at most 22 bytes. */
        struct {
            uint64_t buf_val;
            size_t buf_len;
            char buf[24];
        };

        /* for printing interesting numbers */
        struct {
            uint64_t prev_prime;
            enum interest prev_interest;
            uint64_t last_printed_index;
            uint64_t next_power_of_2;
            uint64_t next_power_of_10;
            uint64_t next_index_power_of_10;
        };
    };
};

static void update_print_buf(struct print_context *ctx, uint64_t new_val) {
    size_t old_val = ctx->buf_val;
    ctx->buf_val = new_val;
    if (new_val >= old_val && ctx->buf_len >= 2) {
        uint64_t add = new_val - old_val;
        for (size_t i = ctx->buf_len - 1; i--; ) {
            add += ctx->buf[i] - '0';
            ctx->buf[i] = add%10 + '0';
            add /= 10;
            if (add == 0) return;  /* updated */
        }
        /* if we get here, we didn't have enough space */
    }
    ctx->buf_len = snprintf(ctx->buf, sizeof(ctx->buf), "%" PRIu64 "\n", new_val);
}

static int print_callback(void *ctx_arg, uint64_t prime) {
    struct print_context *ctx = ctx_arg;
    if (prime == 0) {
        if (ctx->fp != NULL) fflush(ctx->fp);
        return 0;
    }
    if (prime > ctx->max_prime) return 1;
    if (ctx->count >= ctx->max_count) return 2;
    if (ctx->fp != NULL) {
        update_print_buf(ctx, prime);
        fwrite(ctx->buf, 1, ctx->buf_len, ctx->fp);
    }
    ++ctx->count;
    return 0;
}

static void print_interesting(struct print_context *ctx) {
    if (ctx->prev_prime > 0 && ctx->prev_interest != INT_NONE && ctx->fp != NULL) {
        if (ctx->count > ctx->last_printed_index + 1) fputs("\n", ctx->fp);
        fprintf(ctx->fp, "%20" PRIu64 ".%c %20" PRIu64 " (dec)%c %20" PRIx64 " (hex)%c\n",
                ctx->count,      (ctx->prev_interest & INT_IDX) ? '*' : ' ',
                ctx->prev_prime, (ctx->prev_interest & INT_DEC) ? '*' : ' ',
                ctx->prev_prime, (ctx->prev_interest & INT_HEX) ? '*' : ' ');
        fflush(ctx->fp);
        ctx->last_printed_index = ctx->count;
    }
    ctx->prev_interest = INT_NONE;
}

/* Prints only the primes that are slightly interesting, where interesting is defined as one
of the following:

  1. The index of the prime is within 3 of a power of 10 (e.g. index 97, 98, 99, 100, 101, 102).
  2. The prime is the highest below or the lowest above a power of 10 (e.g. 9973, 10007).
  3. The prime is the highest below or the lowest above a power of 2 (e.g. 0xfff1, 0x10001).

This function labels the interesting column with an asterisk. */
static int print_interesting_callback(void *ctx_arg, uint64_t prime) {
    struct print_context *ctx = ctx_arg;
    if (prime == 0) {
        if (ctx->fp != NULL) fflush(ctx->fp);
        return 0;
    }
    if (prime > ctx->max_prime || ctx->count >= ctx->max_count) return 1;

    enum interest next_interest = INT_NONE;
    if (ctx->next_index_power_of_10 + 2 - ctx->count <= 6) {
        next_interest = INT_IDX;
        if (ctx->next_index_power_of_10 + 2 == ctx->count) {
            ctx->next_index_power_of_10 = ctx->next_index_power_of_10 == 0 ? 10 :
                10*ctx->next_index_power_of_10;
        }
    }
    if (prime > ctx->next_power_of_2) {
        ctx->prev_interest |= INT_HEX;
        next_interest |= INT_HEX;
        ctx->next_power_of_2 = ctx->next_power_of_2 == 0 ? 2 : ctx->next_power_of_2 * 2;
    }
    if (prime > ctx->next_power_of_10) {
        ctx->prev_interest |= INT_DEC;
        next_interest |= INT_DEC;
        ctx->next_power_of_10 = ctx->next_power_of_10 == 0 ? 10 : ctx->next_power_of_10 * 10;
    }
    print_interesting(ctx);
    ctx->prev_interest = next_interest;
    ctx->prev_prime    = prime;
    ctx->count++;
    return 0;
}

static uint64_t arg_min_prime = 0;
static uint64_t arg_max_prime = UINT64_MAX;
static uint64_t arg_max_count = UINT64_MAX;
static bool arg_count_only = false;
static bool arg_interesting_only = false;
static const char *cache_filename = NULL;

bool parse_uint64(const char *arg, uint64_t *val) {
    /* strtoull() incorrectly converts negative numbers to unsigned */
    if (*arg == '\0' || *arg == '-') return false;
    char *end;
    errno = 0;
    /* this assumes unsigned long long is at least as large as uint64_t */
    *val = strtoull(arg, &end, 10);
    return errno == 0 && *end == '\0' && *val <= UINT64_MAX;
}

static const char *usage =
"Usage:\n"
"  primes               prints all primes below 2^64\n"
"  primes -h            prints this help, then exits\n"
"  primes <min>         prints primes starting from min\n"
"  primes <min> <max>   prints primes between min and max\n"
"\n"
"Options:\n"
"  -c                   count only (don't print)\n"
"  -C <filename>        set cache filename (overrides environment)\n"
"  -I                   print interesting primes only\n"
"  -n <n>               generate up to n primes\n"
"\n"
"Environmental variables:\n"
"  PRIMES_CACHE         cache filename (e.g. /tmp/primes.bin)\n"
;

bool parse_args(int argc, char *argv[]) {
    int have_args = 0;
    bool have_n = false;
    for (int i = 1; i < argc; ) {
        const char *arg = argv[i++];
        if (strcmp(arg, "-h") == 0) {
            fputs(usage, stdout);
            exit(0);
        } else if (strncmp(arg, "-C", 2) == 0) {
            if (cache_filename != NULL) {
                fprintf(stderr, "Duplicate option -C\n");
                return false;
            }
            if (arg[2] != '\0') {
                cache_filename = arg + 2;  /* -n123 */
            } else if (i < argc) {
                cache_filename = argv[i++];  /* -n 123 */
            } else {
                fprintf(stderr, "Missing argument for option -C\n");
                return false;
            }
        } else if (strcmp(arg, "-c") == 0) {
            if (arg_count_only) {
                fprintf(stderr, "Duplicate option -c\n");
                return false;
            }
            arg_count_only = true;
        } else if (strcmp(arg, "-I") == 0) {
            if (arg_interesting_only) {
                fprintf(stderr, "Duplicate option -I\n");
                return false;
            }
            arg_interesting_only = true;
        } else if (strncmp(arg, "-n", 2) == 0) {
            if (have_n) {
                fprintf(stderr, "Duplicate option -n\n");
                return false;
            }
            if (arg[2] != '\0') {
                arg = arg + 2;  /* -n123 */
            } else if (i < argc) {
                arg = argv[i++];  /* -n 123 */
            } else {
                fprintf(stderr, "Missing argument for option -n\n");
                return false;
            }
            if (!parse_uint64(arg, &arg_max_count)) {
                fprintf(stderr, "Invalid argument for option -n: %s\n", arg);
                return false;
            }
            have_n = true;
        } else if (have_args == 0) {
            if (!parse_uint64(arg, &arg_min_prime)) {
                fprintf(stderr, "Invalid min argument: %s\n", arg);
                return false;
            }
            ++have_args;
        } else if (have_args == 1) {
            if (!parse_uint64(arg, &arg_max_prime)) {
                fprintf(stderr, "Invalid max argument: %s\n", arg);
                return false;
            }
            ++have_args;
        } else {
            fprintf(stderr, "Too many arguments.\n");
            return false;
        }
    }
    if (arg_interesting_only) {
        if (arg_count_only) {
            fprintf(stderr, "Cannot specify both -I and -c\n");
            return false;
        }
        if (arg_min_prime > 2) {
            fprintf(stderr, "Cannot raise the minimum when specifying -I\n");
            return false;
        }
    }

    if (arg_min_prime > arg_max_prime) {
        fprintf(stderr, "Warning: min > max (no primes will match)\n");
    }
    if (arg_max_count < 1) {
        fprintf(stderr, "Warning: count = 0 (no primes will match)\n");
    }
    return true;
}

/* Generates the prime number wheel bitmap for primes up to `max`.

If cache_filename is NULL, this will just call generate_small_prime_wheel(max).

If cache_filename is not NULL, this function first tries to load an existing
bitmap from disk. If the cache file does not exist, or is not big enough, then
the wheel is generated anew, and the function tries to create/overwrite the
cache file with a new/bigger bitmap.

Neither failure to read nor failure to write are considered fatal errors,
so the only reason this function could return NULL is when there is not enough
memory to generate the wheel bitmap.
*/
static uint8_t *generate_cached_prime_wheel(
        const char *cache_filename,
        uint64_t max
) {
    if (cache_filename == NULL && max > 10000) {
        fprintf(stderr,
                "WARNING: generating large numbers of primes without cache is slow!\n"
                "Consider setting PRIMES_CACHE or use the -C argument to specify a cache\n"
                "location to speed up future invocations.\n");
    }
    uint8_t *wheel = NULL;
    if (cache_filename != NULL) {
        /* Try to load existing version from disk */
        wheel = load_wheel_bitmap(cache_filename, max);
        if (wheel != NULL) return wheel;  /* done! */
        if (errno == ENOENT) {
            /* file did not exist; don't print a warning, because this is
               expected to happen */
        } else {
            perror("failed to load small prime wheel bitmap (not fatal)");
        }
    }

    wheel = generate_prime_wheel(max);
    if (wheel == NULL) return NULL;  /* out of memory */

    if (cache_filename != NULL) {
        /* Try to save generated version to disk */
        if (save_wheel_bitmap(cache_filename, wheel, max) != 0) {
            perror("failed to save small prime wheel bitmap (not fatal)");
        }
    }
    return wheel;
}

/* Returns an upper bound on the maximum prime number generated, assuming that
we will generate at most the first `count` primes between `min` and `max`
(inclusive). This is useful to determine up to which number we need to sieve,
when `max` is very large but `min` and `count` are smaller.

Uses the following bounds to approximate:

    p(n) > n(ln n + ln ln n - 1) for n ≥ 6 (Dusart)
    p(n) < n(ln n + ln ln n    ) for n ≥ 6 (Rosser)

    pi(x) > x / (log(x) - 1  )  for x >=  5393 (Dusart)
    pi(x) < x / (log(x) - 1.1)  for x >= 60184 (Dusart)

https://en.wikipedia.org/wiki/Prime-counting_function#Inequalities
*/
uint64_t estimate_max_ub(uint64_t min, uint64_t max, uint64_t count) {
    if (max < 20) return max;  /* small enough already */

    // max_count is an upper bound on the number of primes below 2**64;
    // if count is greater than that, we will hit the `max` limit first.
    // Calculated as: ceil((2**64) / (log(2**64) - 1.1))
    const uint64_t max_count = 425418361821749888;
    if (count > max_count) return max;

    // upper bound on the number of primes below `min`.
    uint64_t skipped = min/2;
    if (min >= 4000) {
        skipped = ceil(min / (log(min) - 1.1));
    } else if (min >= 10) {
        skipped = ceil(min / log(min) * 1.25506);
    }
    if (skipped >= UINT64_MAX - count) return max;  /* prevent overflow */

    // i = index of the largest prime we will generate.
    uint64_t i = count + skipped;
    if (i <= 8) return 19;  /* p(8) = 19 */
    double lc = log(i);
    double ub = (lc + log(lc)) * i + 1;
    if (ub < max) return ub;
    return max;
}

int main(int argc, char *argv[]) {
    if (!parse_args(argc, argv)) {
        putc('\n', stderr);
        fputs(usage, stderr);
        exit(1);
    }

    if (cache_filename == NULL) {
        /* Cache filename was not specified on the command line.
           Try getting it from the environment instead. */
        cache_filename = getenv("PRIMES_CACHE");
    }

#ifndef NDEBUG
    static_assertions();
#endif
    if (arg_max_prime < arg_min_prime) return 0;

    struct print_context context = {
        .fp         = arg_count_only ? NULL : stdout,
        .count      = 0,
        .max_prime  = arg_max_prime,
        .max_count  = arg_max_count,
    };
    enumerate_prime_callback_t callback = print_callback;
    if (arg_interesting_only) callback = print_interesting_callback;

    uint64_t max_ub = estimate_max_ub(arg_min_prime, arg_max_prime, arg_max_count);
    if (max_ub < 2000000) {
        /* Only primes below 2 million. Sieving in memory is fast enough. */
        enumerate_small_primes_sieved(arg_min_prime, max_ub, callback, &context);
    } else {
        /* Use wheel bitmap if we can. Don't generate larger than 32-bit integers.
           Round up to 2^k - 1 to avoid repeated overwrites with tiny increases. */
        uint32_t wheel_max = UINT32_MAX;
        while ((wheel_max >> 1) >= max_ub) wheel_max >>= 1;
        const uint8_t *wheel = generate_cached_prime_wheel(cache_filename, wheel_max);
        if (wheel == NULL) {
            perror("failed to generate wheel bitmap");
            exit(1);
        } else {
            enumerate_primes_sieved(
                    wheel, arg_min_prime, max_ub,
                    callback, &context, true);
            free((void*) wheel);
        }
    }

    if (arg_interesting_only) {
        /* Print the final interesting number. This is necessary because
           print_interesting_callback lags behind by one prime, so for example
           `primes -I 1 20` would not print 19 without calling
           print_interesting() here.

           Note that it's still possible to miss the final line. For example, 13
           is considered interesting because it's the last four-bit prime (1101
           in binary), but that's only discovered when 17 is enumerated, so e.g.
           `primes -I 1 15` won't show 13 as being interesting. I don't have a
           good solution for this, but I don't think it matters much in practice.
        */
        print_interesting(&context);
    }

    if (arg_count_only) {
        printf("%" PRIu64 "\n", context.count);
    }

    return 0;
}
