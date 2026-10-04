#include "leechkem_ref.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bench_random.h"
#include "fips202.h"
#include "leech.h"
#include "leech_utils.h"

#define LEECH_DIM 24U
#define SEED_A_BYTES 16U
#define MAX_SS_BYTES 32U

static const uint16_t cdf_l1[] = {
    4643, 13363, 20579, 25843, 29227, 31145, 32103,
    32525, 32689, 32745, 32762, 32766, 32767
};
static const uint16_t cdf_l2[] = {
    5638, 15915, 23689, 28571, 31116, 32217,
    32613, 32731, 32760, 32766, 32767
};
static const uint16_t cdf_l3[] = {
    9142, 23462, 30338, 32361, 32725, 32765, 32767
};

static const lk_params parameter_sets[] = {
    {"LeechKEM-1", 880,  1, 18, 13, 5, 132, 256, 128, 128},
    {"LeechKEM-2", 1136, 2, 17, 13, 4, 216, 384, 192, 256},
    {"LeechKEM-3", 1640, 2, 17, 12, 5, 264, 512, 256, 256}
};

struct lk_ctx {
    const lk_params *p;
    uint32_t q;
    uint32_t qmask;
    uint32_t beta;
    size_t mu_bytes;
    size_t salt_bytes;
    size_t ss_bytes;
    size_t seedse_bytes;
    size_t pk_b_bytes;
    size_t pk_bytes;
    size_t ct_u_bytes;
    size_t ct_c_bytes;
    size_t ct0_bytes;
    size_t ct_bytes;

    uint8_t *pk;
    uint8_t *ct;
    uint8_t *ct_bad;
    uint8_t *ct_check;
    uint8_t *z;
    uint8_t *mu;
    uint8_t *mu_prime;
    uint8_t *seedse;
    uint8_t ss_enc[MAX_SS_BYTES];
    uint8_t ss_dec[MAX_SS_BYTES];

    int16_t *S;
    int16_t *E;
    int16_t *Sp;
    int16_t *Ep;
    int16_t *Epp;

    uint32_t *B;
    uint32_t *U;
    uint32_t *C;
    uint32_t *V;
    uint32_t *O;
    uint32_t *N;
    uint32_t *Ucheck;
    uint32_t *Ccheck;
    uint32_t *Arow;
    int64_t *acc;

    uint8_t *rowbuf;
    size_t rowbuf_bytes;
    uint8_t *noisebuf;
    size_t noisebuf_bytes;
    uint8_t *hashbuf;
    size_t hashbuf_bytes;
};

static size_t ceil_div(size_t a, size_t b)
{
    return (a + b - 1U) / b;
}

size_t lk_parameter_count(void)
{
    return sizeof(parameter_sets) / sizeof(parameter_sets[0]);
}

const lk_params *lk_parameter_at(size_t index)
{
    return index < lk_parameter_count() ? &parameter_sets[index] : NULL;
}

static const uint16_t *cdf_for(const lk_params *p, size_t *length)
{
    if (p == &parameter_sets[0] || p->n == 880) {
        *length = sizeof(cdf_l1) / sizeof(cdf_l1[0]);
        return cdf_l1;
    }
    if (p == &parameter_sets[1] || p->n == 1136) {
        *length = sizeof(cdf_l2) / sizeof(cdf_l2[0]);
        return cdf_l2;
    }
    *length = sizeof(cdf_l3) / sizeof(cdf_l3[0]);
    return cdf_l3;
}

static void xof(const lk_ctx *ctx, uint8_t *out, size_t outlen,
                const uint8_t *in, size_t inlen)
{
    if (ctx->p->shake_bits == 128) {
        shake128(out, outlen, in, inlen);
    } else {
        shake256(out, outlen, in, inlen);
    }
}

static void xof_public_matrix(uint8_t *out, size_t outlen,
                              const uint8_t *in, size_t inlen)
{
    /*
     * seed_A is 128 bits for every parameter set.  Following the FrodoKEM
     * SHAKE variant, public-matrix expansion therefore uses SHAKE128 at all
     * three levels; level-dependent SHAKE is retained for secrets and hashes.
     */
    shake128(out, outlen, in, inlen);
}

static void clear_unused_bits(uint8_t *bytes, unsigned bits)
{
    unsigned remainder = bits & 7U;
    if (remainder != 0) {
        bytes[(bits - 1U) / 8U] &= (uint8_t)((1U << remainder) - 1U);
    }
}

static unsigned get_bits(const uint8_t *bytes, size_t bit_offset, unsigned count)
{
    unsigned value = 0;
    unsigned i;
    for (i = 0; i < count; ++i) {
        value |= (unsigned)((bytes[(bit_offset + i) >> 3] >>
                  ((bit_offset + i) & 7U)) & 1U) << i;
    }
    return value;
}

static void put_bits(uint8_t *bytes, size_t bit_offset, unsigned count,
                     unsigned value)
{
    unsigned i;
    for (i = 0; i < count; ++i) {
        size_t position = bit_offset + i;
        uint8_t mask = (uint8_t)(1U << (position & 7U));
        uint8_t bit = (uint8_t)((value >> i) & 1U);
        bytes[position >> 3] =
            (uint8_t)((bytes[position >> 3] & (uint8_t)~mask) |
                      (uint8_t)((0U - bit) & mask));
    }
}

static void pack_values(uint8_t *out, size_t outlen, const uint32_t *values,
                        size_t count, unsigned bits)
{
    uint64_t reservoir = 0;
    unsigned available = 0;
    size_t i = 0, j = 0;
    uint32_t mask = (UINT32_C(1) << bits) - 1U;
    memset(out, 0, outlen);
    while (i < count) {
        reservoir |= ((uint64_t)(values[i++] & mask)) << available;
        available += bits;
        while (available >= 8 && j < outlen) {
            out[j++] = (uint8_t)reservoir;
            reservoir >>= 8;
            available -= 8;
        }
    }
    if (available != 0 && j < outlen) {
        out[j] = (uint8_t)reservoir;
    }
}

static void unpack_values(uint32_t *values, size_t count, const uint8_t *in,
                          size_t inlen, unsigned bits)
{
    uint64_t reservoir = 0;
    unsigned available = 0;
    size_t i = 0, j = 0;
    uint32_t mask = (UINT32_C(1) << bits) - 1U;
    while (i < count) {
        while (available < bits && j < inlen) {
            reservoir |= ((uint64_t)in[j++]) << available;
            available += 8;
        }
        values[i++] = (uint32_t)reservoir & mask;
        reservoir >>= bits;
        available -= bits;
    }
}

static int allocate_all(lk_ctx *ctx)
{
    const size_t n24 = ctx->p->n * LEECH_DIM;
    const size_t mn = ctx->p->mbar * ctx->p->n;
    const size_t m24 = ctx->p->mbar * LEECH_DIM;
    const size_t max_noise_count = n24 > (2U * mn + m24) ?
                                   n24 : (2U * mn + m24);

#define ALLOC_FIELD(field, count, type) \
    do { \
        ctx->field = (type *)calloc((count), sizeof(type)); \
        if (ctx->field == NULL) return 0; \
    } while (0)

    ALLOC_FIELD(pk, ctx->pk_bytes, uint8_t);
    ALLOC_FIELD(ct, ctx->ct_bytes, uint8_t);
    ALLOC_FIELD(ct_bad, ctx->ct_bytes, uint8_t);
    ALLOC_FIELD(ct_check, ctx->ct0_bytes, uint8_t);
    ALLOC_FIELD(z, ctx->mu_bytes, uint8_t);
    ALLOC_FIELD(mu, ctx->mu_bytes, uint8_t);
    ALLOC_FIELD(mu_prime, ctx->mu_bytes, uint8_t);
    ALLOC_FIELD(seedse, ctx->seedse_bytes, uint8_t);

    ALLOC_FIELD(S, n24, int16_t);
    ALLOC_FIELD(E, n24, int16_t);
    ALLOC_FIELD(Sp, mn, int16_t);
    ALLOC_FIELD(Ep, mn, int16_t);
    ALLOC_FIELD(Epp, m24, int16_t);

    ALLOC_FIELD(B, n24, uint32_t);
    ALLOC_FIELD(U, mn, uint32_t);
    ALLOC_FIELD(C, m24, uint32_t);
    ALLOC_FIELD(V, m24, uint32_t);
    ALLOC_FIELD(O, m24, uint32_t);
    ALLOC_FIELD(N, m24, uint32_t);
    ALLOC_FIELD(Ucheck, mn, uint32_t);
    ALLOC_FIELD(Ccheck, m24, uint32_t);
    ALLOC_FIELD(Arow, ctx->p->n, uint32_t);
    ALLOC_FIELD(acc, mn, int64_t);

    ctx->rowbuf_bytes = ceil_div(ctx->p->n * ctx->p->logq, 8U);
    ALLOC_FIELD(rowbuf, ctx->rowbuf_bytes, uint8_t);
    ctx->noisebuf_bytes = 2U * max_noise_count;
    ALLOC_FIELD(noisebuf, ctx->noisebuf_bytes, uint8_t);
    ctx->hashbuf_bytes = 1U + ctx->pk_bytes + ctx->mu_bytes + ctx->ct_bytes;
    ALLOC_FIELD(hashbuf, ctx->hashbuf_bytes, uint8_t);

#undef ALLOC_FIELD
    return 1;
}

int lk_ctx_create(lk_ctx **out, const lk_params *params)
{
    lk_ctx *ctx;
    if (out == NULL || params == NULL) return 0;
    ctx = (lk_ctx *)calloc(1, sizeof(*ctx));
    if (ctx == NULL) return 0;
    ctx->p = params;
    ctx->q = UINT32_C(1) << params->logq;
    ctx->qmask = ctx->q - 1U;
    ctx->beta = UINT32_C(1) << params->logbeta;
    ctx->mu_bytes = ceil_div(params->mu_bits, 8U);
    ctx->salt_bytes = params->salt_bits / 8U;
    ctx->ss_bytes = params->ss_bits / 8U;
    ctx->seedse_bytes = ctx->ss_bytes;
    ctx->pk_b_bytes = ceil_div(params->n * LEECH_DIM * params->logq, 8U);
    ctx->pk_bytes = SEED_A_BYTES + ctx->pk_b_bytes;
    ctx->ct_u_bytes = ceil_div(params->mbar * params->n * params->logq, 8U);
    ctx->ct_c_bytes = ceil_div(params->mbar * LEECH_DIM * params->logq, 8U);
    ctx->ct0_bytes = ctx->ct_u_bytes + ctx->ct_c_bytes;
    ctx->ct_bytes = ctx->ct0_bytes + ctx->salt_bytes;

    if (!allocate_all(ctx)) {
        lk_ctx_destroy(ctx);
        return 0;
    }
    if (!leech_set_modulus(8U * ctx->beta)) {
        lk_ctx_destroy(ctx);
        return 0;
    }
    *out = ctx;
    return 1;
}

void lk_ctx_destroy(lk_ctx *ctx)
{
    if (ctx == NULL) return;
#define FREE_FIELD(field) do { free(ctx->field); ctx->field = NULL; } while (0)
    FREE_FIELD(pk); FREE_FIELD(ct); FREE_FIELD(ct_bad); FREE_FIELD(ct_check);
    FREE_FIELD(z); FREE_FIELD(mu); FREE_FIELD(mu_prime); FREE_FIELD(seedse);
    FREE_FIELD(S); FREE_FIELD(E); FREE_FIELD(Sp); FREE_FIELD(Ep);
    FREE_FIELD(Epp); FREE_FIELD(B); FREE_FIELD(U); FREE_FIELD(C);
    FREE_FIELD(V); FREE_FIELD(O); FREE_FIELD(N); FREE_FIELD(Ucheck);
    FREE_FIELD(Ccheck); FREE_FIELD(Arow); FREE_FIELD(acc); FREE_FIELD(rowbuf);
    FREE_FIELD(noisebuf); FREE_FIELD(hashbuf);
#undef FREE_FIELD
    free(ctx);
}

static void generate_a_row(lk_ctx *ctx, const uint8_t seed_a[SEED_A_BYTES],
                           uint32_t row)
{
    uint8_t input[1 + 4 + SEED_A_BYTES];
    input[0] = 0xA0;
    input[1] = (uint8_t)row;
    input[2] = (uint8_t)(row >> 8);
    input[3] = (uint8_t)(row >> 16);
    input[4] = (uint8_t)(row >> 24);
    memcpy(input + 5, seed_a, SEED_A_BYTES);
    xof_public_matrix(ctx->rowbuf, ctx->rowbuf_bytes, input, sizeof(input));
    unpack_values(ctx->Arow, ctx->p->n, ctx->rowbuf,
                  ctx->rowbuf_bytes, ctx->p->logq);
}

static void sample_noise(lk_ctx *ctx, int16_t *out, size_t count,
                         uint8_t domain, const uint8_t *seed, size_t seedlen)
{
    uint8_t input[1 + MAX_SS_BYTES];
    const uint16_t *cdf;
    size_t cdf_len, i, j;
    cdf = cdf_for(ctx->p, &cdf_len);
    input[0] = domain;
    memcpy(input + 1, seed, seedlen);
    xof(ctx, ctx->noisebuf, 2U * count, input, 1U + seedlen);
    for (i = 0; i < count; ++i) {
        uint16_t word = (uint16_t)ctx->noisebuf[2U * i] |
                        (uint16_t)((uint16_t)ctx->noisebuf[2U * i + 1U] << 8);
        uint16_t magnitude = 0;
        uint16_t threshold = word >> 1;
        uint16_t sign = word & 1U;
        for (j = 0; j + 1U < cdf_len; ++j) {
            magnitude = (uint16_t)(magnitude +
                        ((uint16_t)(cdf[j] - threshold) >> 15));
        }
        out[i] = sign ? -(int16_t)magnitude : (int16_t)magnitude;
    }
}

static void matrix_as_plus_e(lk_ctx *ctx)
{
    size_t i, j, k;
    const size_t n = ctx->p->n;
    for (i = 0; i < n; ++i) {
        generate_a_row(ctx, ctx->pk, (uint32_t)i);
        for (j = 0; j < LEECH_DIM; ++j) {
            int64_t sum = ctx->E[i * LEECH_DIM + j];
            for (k = 0; k < n; ++k) {
                sum += (int64_t)ctx->Arow[k] *
                       ctx->S[k * LEECH_DIM + j];
            }
            ctx->B[i * LEECH_DIM + j] = (uint32_t)sum & ctx->qmask;
        }
    }
}

static void matrix_sa_plus_e(lk_ctx *ctx, uint32_t *out)
{
    size_t r, j, k;
    const size_t n = ctx->p->n;
    const size_t mbar = ctx->p->mbar;
    for (r = 0; r < mbar; ++r) {
        for (j = 0; j < n; ++j) {
            ctx->acc[r * n + j] = ctx->Ep[r * n + j];
        }
    }
    for (k = 0; k < n; ++k) {
        generate_a_row(ctx, ctx->pk, (uint32_t)k);
        for (r = 0; r < mbar; ++r) {
            int16_t scalar = ctx->Sp[r * n + k];
            for (j = 0; j < n; ++j) {
                ctx->acc[r * n + j] += (int64_t)scalar * ctx->Arow[j];
            }
        }
    }
    for (r = 0; r < mbar; ++r) {
        for (j = 0; j < n; ++j) {
            out[r * n + j] = (uint32_t)ctx->acc[r * n + j] & ctx->qmask;
        }
    }
}

static void matrix_sb_plus_e(lk_ctx *ctx, uint32_t *out)
{
    size_t r, j, k;
    const size_t n = ctx->p->n;
    for (r = 0; r < ctx->p->mbar; ++r) {
        for (j = 0; j < LEECH_DIM; ++j) {
            int64_t sum = ctx->Epp[r * LEECH_DIM + j];
            for (k = 0; k < n; ++k) {
                sum += (int64_t)ctx->Sp[r * n + k] *
                       ctx->B[k * LEECH_DIM + j];
            }
            out[r * LEECH_DIM + j] = (uint32_t)sum & ctx->qmask;
        }
    }
}

static void matrix_us(lk_ctx *ctx)
{
    size_t r, j, k;
    const size_t n = ctx->p->n;
    for (r = 0; r < ctx->p->mbar; ++r) {
        for (j = 0; j < LEECH_DIM; ++j) {
            int64_t sum = 0;
            for (k = 0; k < n; ++k) {
                sum += (int64_t)ctx->U[r * n + k] *
                       ctx->S[k * LEECH_DIM + j];
            }
            ctx->N[r * LEECH_DIM + j] =
                (ctx->C[r * LEECH_DIM + j] - (uint32_t)sum) & ctx->qmask;
        }
    }
}

/* Exact coefficient allocation from Table 1 in the manuscript.  Index i is
 * the zero-based coefficient a_{i+1}; a_1 is fixed to zero.  The remaining
 * widths consume b_1,...,b_36 consecutively. */
static const uint8_t f_coefficient_width[LEECH_DIM] = {
    0, 1,1,1,1,1,1,2,1,1,1,2,1,2,2,2,1,2,2,2,2,2,2,3
};

static int mog_generator_entry(size_t row, size_t column)
{
    /* The upstream table uses uint8_t storage; its sole negative entry is -3. */
    if (row == 23U && column == 0U) return -3;
    return (int)leech[row][column];
}

static void f_message_to_coefficients(const uint8_t *message,
                                      size_t bit_offset,
                                      unsigned coefficients[LEECH_DIM])
{
    size_t i;
    coefficients[0] = 0;
    for (i = 1; i < LEECH_DIM; ++i) {
        unsigned width = f_coefficient_width[i];
        coefficients[i] = get_bits(message, bit_offset, width);
        bit_offset += width;
    }
}

static int positive_mod8(int64_t value)
{
    int result = (int)(value % 8);
    return result < 0 ? result + 8 : result;
}

/* The paper orders the g_tau bits by bit plane: b_{g,24r+j} is the
 * r-th binary digit assigned to coordinate j. */
static unsigned g_tau_digit_from_message(const uint8_t *message,
                                         size_t block_offset,
                                         unsigned tau, unsigned coordinate)
{
    unsigned digit = 0;
    unsigned r;
    for (r = 0; r < tau - 1U; ++r) {
        digit |= get_bits(message,
                          block_offset + 36U + LEECH_DIM * r + coordinate,
                          1U) << r;
    }
    return digit;
}

static void g_tau_digit_to_message(uint8_t *message, size_t block_offset,
                                   unsigned tau, unsigned coordinate,
                                   unsigned digit)
{
    unsigned r;
    for (r = 0; r < tau - 1U; ++r) {
        put_bits(message,
                 block_offset + 36U + LEECH_DIM * r + coordinate,
                 1U, (digit >> r) & 1U);
    }
}

/* Coordinate isometry documented in the manuscript: swap every adjacent
 * pair.  The corrected standard MOG matrix in leech_utils.c makes all 24
 * generator rows map to distance-zero Vardy--Be'ery decoder points. */
static void mog_to_decoder_coordinates(const uint32_t source[LEECH_DIM],
                                       uint32_t destination[LEECH_DIM])
{
    unsigned j;
    for (j = 0; j < LEECH_DIM / 2U; ++j) {
        destination[2U * j] = source[2U * j + 1U];
        destination[2U * j + 1U] = source[2U * j];
    }
}

static void decoder_to_mog_coordinates(const uint32_t source[LEECH_DIM],
                                       uint32_t destination[LEECH_DIM])
{
    unsigned j;
    for (j = 0; j < LEECH_DIM / 2U; ++j) {
        destination[2U * j] = source[2U * j + 1U];
        destination[2U * j + 1U] = source[2U * j];
    }
}

/* Invert a*G_MOG^int modulo 8 by back substitution.  The generator is lower
 * triangular and the coefficient widths give one permitted solution for each
 * diagonal congruence. */
static int f_point_to_message(uint8_t *message, size_t bit_offset,
                              const uint32_t point_scaled[LEECH_DIM],
                              uint32_t beta)
{
    unsigned coefficients[LEECH_DIM] = {0};
    unsigned point_mod8[LEECH_DIM];
    int i;
    size_t column, row;

    for (column = 0; column < LEECH_DIM; ++column) {
        if (point_scaled[column] % beta != 0U) return 0;
        point_mod8[column] = (point_scaled[column] / beta) & 7U;
    }

    for (i = (int)LEECH_DIM - 1; i >= 1; --i) {
        int64_t remainder = point_mod8[(size_t)i];
        unsigned diagonal;
        for (row = (size_t)i + 1U; row < LEECH_DIM; ++row) {
            remainder -= (int64_t)coefficients[row] *
                         mog_generator_entry(row, (size_t)i);
        }
        diagonal = (unsigned)mog_generator_entry((size_t)i, (size_t)i);
        remainder = positive_mod8(remainder);
        if (diagonal == 0U || (unsigned)remainder % diagonal != 0U) return 0;
        coefficients[(size_t)i] = (unsigned)remainder / diagonal;
        if (coefficients[(size_t)i] >=
            (1U << f_coefficient_width[(size_t)i])) return 0;
    }

    /* a_1 is fixed to zero, so coordinate 1 provides a consistency check. */
    {
        int64_t reconstructed = 0;
        for (row = 1; row < LEECH_DIM; ++row) {
            reconstructed += (int64_t)coefficients[row] *
                             mog_generator_entry(row, 0);
        }
        if ((unsigned)positive_mod8(reconstructed) != point_mod8[0]) return 0;
    }

    for (row = 1; row < LEECH_DIM; ++row) {
        unsigned width = f_coefficient_width[row];
        put_bits(message, bit_offset, width, coefficients[row]);
        bit_offset += width;
    }
    return 1;
}

static void lattice_encode(const lk_ctx *ctx, const uint8_t *message,
                           uint32_t *encoded)
{
    size_t block;
    const unsigned block_bits = 36U + LEECH_DIM * (ctx->p->tau - 1U);
    for (block = 0; block < ctx->p->mbar; ++block) {
        unsigned coefficients[LEECH_DIM];
        unsigned j, i;
        f_message_to_coefficients(message, block * block_bits, coefficients);
        for (j = 0; j < LEECH_DIM; ++j) {
            int64_t mog_integer_coordinate = 0;
            unsigned residual = g_tau_digit_from_message(
                message, block * block_bits, ctx->p->tau, j);
            for (i = 1; i < LEECH_DIM; ++i) {
                mog_integer_coordinate += (int64_t)coefficients[i] *
                                          mog_generator_entry(i, j);
            }
            /* f is defined modulo 2 Z^24, so reduce a*G_MOG^int modulo 8
             * before selecting its canonical representative and adding the
             * g_tau digit.  This matters for the -3 entry in generator row
             * 24 and for larger coefficient combinations. */
            mog_integer_coordinate = positive_mod8(mog_integer_coordinate);
            /* beta*(a*G_MOG^int/4 + 2*residual), reduced modulo q. */
            mog_integer_coordinate =
                (int64_t)(ctx->beta / 4U) * mog_integer_coordinate +
                (int64_t)(2U * ctx->beta) * residual;
            encoded[block * LEECH_DIM + j] =
                (uint32_t)mog_integer_coordinate & ctx->qmask;
        }
    }
}

static int lattice_decode(lk_ctx *ctx, const uint32_t *received,
                          uint8_t *message)
{
    size_t block;
    const unsigned block_bits = 36U + LEECH_DIM * (ctx->p->tau - 1U);
    const uint32_t decoder_q = 8U * ctx->beta;
    const uint64_t residual_modulus = UINT64_C(4) * ctx->q;
    const uint64_t quant_step = UINT64_C(8) * ctx->beta;
    const uint64_t half_step = UINT64_C(4) * ctx->beta;
    const unsigned residual_mask = (1U << (ctx->p->tau - 1U)) - 1U;
    memset(message, 0, ctx->mu_bytes);
    if (!leech_set_modulus(decoder_q)) return 0;

    for (block = 0; block < ctx->p->mbar; ++block) {
        uint32_t target_mog[LEECH_DIM], target[LEECH_DIM];
        uint32_t point[LEECH_DIM], point_original[LEECH_DIM];
        uint64_t cv, distance;
        unsigned j;
        const uint32_t *row = received + block * LEECH_DIM;

        for (j = 0; j < LEECH_DIM; ++j) {
            target_mog[j] = (uint32_t)
                ((UINT64_C(4) * row[j]) % decoder_q);
        }
        mog_to_decoder_coordinates(target_mog, target);
        decoder_L24(target, &cv, &distance);
        decode_pt(point, cv, 0);
        (void)distance;
        decoder_to_mog_coordinates(point, point_original);
        if (!f_point_to_message(message, block * block_bits,
                                point_original, ctx->beta)) return 0;
        for (j = 0; j < LEECH_DIM; ++j) {
            uint64_t scaled_residual =
                (UINT64_C(4) * row[j] + residual_modulus -
                 point_original[j]) % residual_modulus;
            unsigned digit = (unsigned)((scaled_residual + half_step) /
                             quant_step) & residual_mask;
            g_tau_digit_to_message(message, block * block_bits,
                                   ctx->p->tau, j, digit);
        }
    }
    clear_unused_bits(message, ctx->p->mu_bits);
    return 1;
}

static void hash_g(lk_ctx *ctx, const uint8_t *message, const uint8_t *salt,
                   uint8_t *out)
{
    size_t offset = 0;
    ctx->hashbuf[offset++] = 0x47;
    memcpy(ctx->hashbuf + offset, ctx->pk, ctx->pk_bytes);
    offset += ctx->pk_bytes;
    memcpy(ctx->hashbuf + offset, message, ctx->mu_bytes);
    offset += ctx->mu_bytes;
    memcpy(ctx->hashbuf + offset, salt, ctx->salt_bytes);
    offset += ctx->salt_bytes;
    xof(ctx, out, ctx->seedse_bytes, ctx->hashbuf, offset);
}

static void hash_h(lk_ctx *ctx, const uint8_t *message, const uint8_t *ct,
                   uint8_t *out)
{
    size_t offset = 0;
    ctx->hashbuf[offset++] = 0x48;
    memcpy(ctx->hashbuf + offset, ctx->pk, ctx->pk_bytes);
    offset += ctx->pk_bytes;
    memcpy(ctx->hashbuf + offset, message, ctx->mu_bytes);
    offset += ctx->mu_bytes;
    memcpy(ctx->hashbuf + offset, ct, ctx->ct_bytes);
    offset += ctx->ct_bytes;
    xof(ctx, out, ctx->ss_bytes, ctx->hashbuf, offset);
}

static void expand_encryption_noise(lk_ctx *ctx, const uint8_t *seed)
{
    size_t mn = ctx->p->mbar * ctx->p->n;
    size_t m24 = ctx->p->mbar * LEECH_DIM;
    sample_noise(ctx, ctx->Sp, mn, 0x51, seed, ctx->seedse_bytes);
    sample_noise(ctx, ctx->Ep, mn, 0x52, seed, ctx->seedse_bytes);
    sample_noise(ctx, ctx->Epp, m24, 0x53, seed, ctx->seedse_bytes);
}

static void unpack_public_key(lk_ctx *ctx)
{
    unpack_values(ctx->B, ctx->p->n * LEECH_DIM,
                  ctx->pk + SEED_A_BYTES, ctx->pk_b_bytes, ctx->p->logq);
}

static void pke_encrypt_core(lk_ctx *ctx, const uint8_t *message,
                             const uint8_t *seed, uint32_t *u_out,
                             uint32_t *c_out, uint8_t *packed_ct0)
{
    size_t i;
    expand_encryption_noise(ctx, seed);
    unpack_public_key(ctx);
    matrix_sa_plus_e(ctx, u_out);
    matrix_sb_plus_e(ctx, ctx->V);
    lattice_encode(ctx, message, ctx->O);
    for (i = 0; i < ctx->p->mbar * LEECH_DIM; ++i) {
        c_out[i] = (ctx->V[i] + ctx->O[i]) & ctx->qmask;
    }
    pack_values(packed_ct0, ctx->ct_u_bytes, u_out,
                ctx->p->mbar * ctx->p->n, ctx->p->logq);
    pack_values(packed_ct0 + ctx->ct_u_bytes, ctx->ct_c_bytes, c_out,
                ctx->p->mbar * LEECH_DIM, ctx->p->logq);
}

int lk_keygen(lk_ctx *ctx)
{
    uint8_t seed_s[MAX_SS_BYTES], seed_e[MAX_SS_BYTES];
    if (ctx == NULL) return 0;
    if (randombytes(ctx->pk, SEED_A_BYTES) != 0 ||
        randombytes(seed_s, ctx->seedse_bytes) != 0 ||
        randombytes(seed_e, ctx->seedse_bytes) != 0 ||
        randombytes(ctx->z, ctx->mu_bytes) != 0) {
        return 0;
    }
    clear_unused_bits(ctx->z, ctx->p->mu_bits);
    sample_noise(ctx, ctx->S, ctx->p->n * LEECH_DIM,
                 0x41, seed_s, ctx->seedse_bytes);
    sample_noise(ctx, ctx->E, ctx->p->n * LEECH_DIM,
                 0x42, seed_e, ctx->seedse_bytes);
    matrix_as_plus_e(ctx);
    pack_values(ctx->pk + SEED_A_BYTES, ctx->pk_b_bytes, ctx->B,
                ctx->p->n * LEECH_DIM, ctx->p->logq);
    return 1;
}

int lk_encaps(lk_ctx *ctx)
{
    uint8_t *salt;
    if (ctx == NULL) return 0;
    salt = ctx->ct + ctx->ct0_bytes;
    if (randombytes(ctx->mu, ctx->mu_bytes) != 0 ||
        randombytes(salt, ctx->salt_bytes) != 0) {
        return 0;
    }
    clear_unused_bits(ctx->mu, ctx->p->mu_bits);
    hash_g(ctx, ctx->mu, salt, ctx->seedse);
    pke_encrypt_core(ctx, ctx->mu, ctx->seedse, ctx->U, ctx->C, ctx->ct);
    hash_h(ctx, ctx->mu, ctx->ct, ctx->ss_enc);
    return 1;
}

static uint8_t constant_time_is_zero(const uint8_t *a, const uint8_t *b,
                                     size_t length)
{
    uint32_t difference = 0;
    size_t i;
    for (i = 0; i < length; ++i) difference |= (uint32_t)(a[i] ^ b[i]);
    return (uint8_t)(((difference | (0U - difference)) >> 31) ^ 1U);
}

static int decaps_core(lk_ctx *ctx, const uint8_t *ct)
{
    uint8_t k0[MAX_SS_BYTES], k1[MAX_SS_BYTES];
    uint8_t valid, mask;
    size_t i;
    const uint8_t *salt = ct + ctx->ct0_bytes;

    unpack_values(ctx->U, ctx->p->mbar * ctx->p->n,
                  ct, ctx->ct_u_bytes, ctx->p->logq);
    unpack_values(ctx->C, ctx->p->mbar * LEECH_DIM,
                  ct + ctx->ct_u_bytes, ctx->ct_c_bytes, ctx->p->logq);
    matrix_us(ctx);
    if (!lattice_decode(ctx, ctx->N, ctx->mu_prime)) return 0;
    hash_g(ctx, ctx->mu_prime, salt, ctx->seedse);
    pke_encrypt_core(ctx, ctx->mu_prime, ctx->seedse,
                     ctx->Ucheck, ctx->Ccheck, ctx->ct_check);
    hash_h(ctx, ctx->mu_prime, ct, k0);
    hash_h(ctx, ctx->z, ct, k1);
    valid = constant_time_is_zero(ctx->ct_check, ct, ctx->ct0_bytes);
    mask = (uint8_t)(0U - valid);
    for (i = 0; i < ctx->ss_bytes; ++i) {
        ctx->ss_dec[i] = (uint8_t)((k0[i] & mask) | (k1[i] & (uint8_t)~mask));
    }
    return 1;
}

int lk_decaps_valid(lk_ctx *ctx)
{
    return ctx != NULL && decaps_core(ctx, ctx->ct);
}

int lk_decaps_invalid(lk_ctx *ctx)
{
    return ctx != NULL && decaps_core(ctx, ctx->ct_bad);
}

int lk_prepare_invalid_ciphertext(lk_ctx *ctx)
{
    if (ctx == NULL) return 0;
    memcpy(ctx->ct_bad, ctx->ct, ctx->ct_bytes);
    ctx->ct_bad[0] ^= 1U;
    return 1;
}

int lk_selftest(lk_ctx *ctx)
{
    uint8_t *test_message;
    unsigned trial;
    size_t bit;
    int ok = 0;
    if (ctx == NULL) return 0;
    test_message = (uint8_t *)calloc(ctx->mu_bytes, 1);
    if (test_message == NULL) return 0;

    /* Check the paper's g_tau indexing directly, independently of decoding.
     * Setting b_{g,24r+j}=1 must add 2^(r+1) beta in coordinate j. */
    {
        unsigned r, j, k;
        for (r = 0; r < ctx->p->tau - 1U; ++r) {
            for (j = 0; j < LEECH_DIM; ++j) {
                uint32_t expected = (uint32_t)
                    (((uint64_t)(2U * ctx->beta) << r) & ctx->qmask);
                memset(ctx->mu, 0, ctx->mu_bytes);
                put_bits(ctx->mu, 36U + LEECH_DIM * r + j, 1U, 1U);
                lattice_encode(ctx, ctx->mu, ctx->O);
                for (k = 0; k < ctx->p->mbar * LEECH_DIM; ++k) {
                    uint32_t coordinate_expected = (k == j) ? expected : 0U;
                    if (ctx->O[k] != coordinate_expected) {
                        fprintf(stderr,
                                "self-test: g_tau formula mismatch at "
                                "r=%u, j=%u, coordinate=%u\n",
                                r, j, k);
                        goto cleanup;
                    }
                }
            }
        }
    }

    /* Zero, every basis bit, and random messages exercise exact f/g inversion. */
    memset(ctx->mu, 0, ctx->mu_bytes);
    lattice_encode(ctx, ctx->mu, ctx->O);
    if (!lattice_decode(ctx, ctx->O, test_message) ||
        memcmp(ctx->mu, test_message, ctx->mu_bytes) != 0) {
        fprintf(stderr, "self-test: zero-message lattice round trip failed\n");
        goto cleanup;
    }
    for (bit = 0; bit < ctx->p->mu_bits; ++bit) {
        memset(ctx->mu, 0, ctx->mu_bytes);
        ctx->mu[bit >> 3] = (uint8_t)(1U << (bit & 7U));
        lattice_encode(ctx, ctx->mu, ctx->O);
        if (!lattice_decode(ctx, ctx->O, test_message) ||
            memcmp(ctx->mu, test_message, ctx->mu_bytes) != 0) {
            size_t debug_bit;
            fprintf(stderr,
                    "self-test: basis-bit lattice round trip failed at %zu\n",
                    bit);
            fprintf(stderr, "self-test: decoded set bits:");
            for (debug_bit = 0; debug_bit < ctx->p->mu_bits; ++debug_bit) {
                if ((test_message[debug_bit >> 3] >> (debug_bit & 7U)) & 1U) {
                    fprintf(stderr, " %zu", debug_bit);
                }
            }
            fprintf(stderr, "\n");
            goto cleanup;
        }
    }
    for (trial = 0; trial < 32U; ++trial) {
        if (randombytes(ctx->mu, ctx->mu_bytes) != 0) {
            fprintf(stderr, "self-test: random message generation failed\n");
            goto cleanup;
        }
        clear_unused_bits(ctx->mu, ctx->p->mu_bits);
        lattice_encode(ctx, ctx->mu, ctx->O);
        if (!lattice_decode(ctx, ctx->O, test_message) ||
            memcmp(ctx->mu, test_message, ctx->mu_bytes) != 0) {
            fprintf(stderr,
                    "self-test: random lattice round trip failed at %u\n",
                    trial);
            goto cleanup;
        }
    }

    if (!lk_keygen(ctx) || !lk_encaps(ctx) || !lk_decaps_valid(ctx)) {
        fprintf(stderr, "self-test: full valid-path execution failed\n");
        goto cleanup;
    }
    if (memcmp(ctx->ss_enc, ctx->ss_dec, ctx->ss_bytes) != 0) {
        fprintf(stderr, "self-test: valid encaps/decaps shared-key mismatch\n");
        goto cleanup;
    }
    if (!lk_prepare_invalid_ciphertext(ctx) || !lk_decaps_invalid(ctx)) {
        fprintf(stderr, "self-test: invalid-path execution failed\n");
        goto cleanup;
    }
    if (memcmp(ctx->ss_enc, ctx->ss_dec, ctx->ss_bytes) == 0) {
        fprintf(stderr, "self-test: invalid ciphertext was not rejected\n");
        goto cleanup;
    }
    ok = 1;

cleanup:
    free(test_message);
    return ok;
}

int lk_shared_secret_matches(const lk_ctx *ctx)
{
    return ctx != NULL &&
           memcmp(ctx->ss_enc, ctx->ss_dec, ctx->ss_bytes) == 0;
}

size_t lk_public_key_bytes(const lk_ctx *ctx)
{
    return ctx != NULL ? ctx->pk_bytes : 0;
}

size_t lk_ciphertext_bytes(const lk_ctx *ctx)
{
    return ctx != NULL ? ctx->ct_bytes : 0;
}

size_t lk_secret_key_bytes_estimate(const lk_ctx *ctx)
{
    if (ctx == NULL) return 0;
    return ctx->pk_bytes + 2U * ctx->p->n * LEECH_DIM + ctx->mu_bytes;
}

uint64_t lk_checksum(const lk_ctx *ctx)
{
    uint64_t value = 0;
    size_t i;
    if (ctx == NULL) return 0;
    for (i = 0; i < ctx->ss_bytes; ++i) {
        value = (value << 5) ^ (value >> 2) ^ ctx->ss_dec[i] ^ ctx->ss_enc[i];
    }
    value ^= ctx->pk[0];
    value ^= (uint64_t)ctx->ct[0] << 32;
    return value;
}
