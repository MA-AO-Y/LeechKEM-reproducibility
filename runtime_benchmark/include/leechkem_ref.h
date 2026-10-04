#ifndef LEECHKEM_REF_H
#define LEECHKEM_REF_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    const char *name;
    size_t n;
    size_t mbar;
    unsigned logq;
    unsigned logbeta;
    unsigned tau;
    unsigned mu_bits;
    unsigned salt_bits;
    unsigned ss_bits;
    unsigned shake_bits;
} lk_params;

typedef struct lk_ctx lk_ctx;

size_t lk_parameter_count(void);
const lk_params *lk_parameter_at(size_t index);

int lk_ctx_create(lk_ctx **out, const lk_params *params);
void lk_ctx_destroy(lk_ctx *ctx);

int lk_keygen(lk_ctx *ctx);
int lk_encaps(lk_ctx *ctx);
int lk_decaps_valid(lk_ctx *ctx);
int lk_decaps_invalid(lk_ctx *ctx);
int lk_prepare_invalid_ciphertext(lk_ctx *ctx);
int lk_selftest(lk_ctx *ctx);
int lk_shared_secret_matches(const lk_ctx *ctx);

size_t lk_public_key_bytes(const lk_ctx *ctx);
size_t lk_ciphertext_bytes(const lk_ctx *ctx);
size_t lk_secret_key_bytes_estimate(const lk_ctx *ctx);
uint64_t lk_checksum(const lk_ctx *ctx);

#endif
