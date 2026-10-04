#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>

#include "leech.h"
#include "leech_utils.h"

static int entry(unsigned row, unsigned column)
{
    if (row == 23U && column == 0U) return -3;
    return (int)leech[row][column];
}

int main(void)
{
    unsigned row, j, failures = 0;
    if (!leech_set_modulus(8U)) return 1;
    for (row = 0; row < 24U; ++row) {
        uint32_t mog[24], vb[24], decoded[24];
        uint64_t cv = 0, distance = 0;
        int same = 1;
        for (j = 0; j < 24U; ++j) {
            int value = entry(row, j) % 8;
            mog[j] = (uint32_t)(value < 0 ? value + 8 : value);
        }
        for (j = 0; j < 12U; ++j) {
            vb[2U * j] = mog[2U * j + 1U];
            vb[2U * j + 1U] = mog[2U * j];
        }
        decoder_L24(vb, &cv, &distance);
        decode_pt(decoded, cv, 0);
        for (j = 0; j < 24U; ++j) {
            if (decoded[j] != vb[j]) same = 0;
        }
        printf("row %2u: exact=%s distance=%" PRIu64 "\n",
               row + 1U, same ? "yes" : "NO", distance);
        if (!same) {
            ++failures;
            printf("  target :");
            for (j = 0; j < 24U; ++j) printf(" %u", vb[j]);
            printf("\n  decoded:");
            for (j = 0; j < 24U; ++j) printf(" %u", decoded[j]);
            printf("\n");
        }
    }
    if (failures == 0U) printf("all 24 generator rows verified\n");
    return failures == 0U ? 0 : 1;
}
