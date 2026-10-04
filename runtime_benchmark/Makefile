CC ?= gcc
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -Wpedantic
FRODO_CFLAGS ?= -O3 -std=gnu11 -Wall -Wextra -Wpedantic

BUILD_DIR := build
LEECH_INC := -Iinclude -Ithird_party/avanpo_leech \
	-Ithird_party/PQCrypto-LWEKE/common/sha3
FRODO_INC := -Iinclude -Ithird_party/PQCrypto-LWEKE/FrodoKEM/src \
	-Ithird_party/PQCrypto-LWEKE/common/sha3 \
	-Ithird_party/PQCrypto-LWEKE/common/random

LEECH_SRC := src/leechkem_bench.c src/leechkem_ref.c src/bench_random.c \
	third_party/avanpo_leech/leech.c \
	third_party/avanpo_leech/leech_utils.c \
	third_party/PQCrypto-LWEKE/common/sha3/fips202.c

FRODO_COMMON := src/frodo_bench.c src/bench_random.c \
	third_party/PQCrypto-LWEKE/FrodoKEM/src/util.c \
	third_party/PQCrypto-LWEKE/common/sha3/fips202.c

.PHONY: all clean

all: $(BUILD_DIR)/leechkem_bench $(BUILD_DIR)/frodo640_bench \
	$(BUILD_DIR)/frodo976_bench $(BUILD_DIR)/frodo1344_bench \
	$(BUILD_DIR)/verify_mog_mapping

$(BUILD_DIR):
	mkdir -p $@

$(BUILD_DIR)/leechkem_bench: $(LEECH_SRC) | $(BUILD_DIR)
	$(CC) $(CFLAGS) $(LEECH_INC) $(LEECH_SRC) -o $@

$(BUILD_DIR)/verify_mog_mapping: tools/verify_mog_mapping.c \
	third_party/avanpo_leech/leech.c \
	third_party/avanpo_leech/leech_utils.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -Ithird_party/avanpo_leech $^ -o $@

$(BUILD_DIR)/frodo640_bench: $(FRODO_COMMON) \
	third_party/PQCrypto-LWEKE/FrodoKEM/src/frodo640.c | $(BUILD_DIR)
	$(CC) $(FRODO_CFLAGS) -DNIX -D_AMD64_ -D_REFERENCE_ \
		-D_SHAKE128_FOR_A_ -DFRODO_LEVEL=640 $(FRODO_INC) \
		$(FRODO_COMMON) third_party/PQCrypto-LWEKE/FrodoKEM/src/frodo640.c -o $@

$(BUILD_DIR)/frodo976_bench: $(FRODO_COMMON) \
	third_party/PQCrypto-LWEKE/FrodoKEM/src/frodo976.c | $(BUILD_DIR)
	$(CC) $(FRODO_CFLAGS) -DNIX -D_AMD64_ -D_REFERENCE_ \
		-D_SHAKE128_FOR_A_ -DFRODO_LEVEL=976 $(FRODO_INC) \
		$(FRODO_COMMON) third_party/PQCrypto-LWEKE/FrodoKEM/src/frodo976.c -o $@

$(BUILD_DIR)/frodo1344_bench: $(FRODO_COMMON) \
	third_party/PQCrypto-LWEKE/FrodoKEM/src/frodo1344.c | $(BUILD_DIR)
	$(CC) $(FRODO_CFLAGS) -DNIX -D_AMD64_ -D_REFERENCE_ \
		-D_SHAKE128_FOR_A_ -DFRODO_LEVEL=1344 $(FRODO_INC) \
		$(FRODO_COMMON) third_party/PQCrypto-LWEKE/FrodoKEM/src/frodo1344.c -o $@

clean:
	rm -f $(BUILD_DIR)/leechkem_bench $(BUILD_DIR)/frodo640_bench \
		$(BUILD_DIR)/frodo976_bench $(BUILD_DIR)/frodo1344_bench \
		$(BUILD_DIR)/verify_mog_mapping
