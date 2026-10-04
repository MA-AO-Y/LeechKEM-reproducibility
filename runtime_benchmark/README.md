# LeechKEM full runtime benchmark (paper-exact v4)

This corrected release follows the manuscript definition
`[g_tau(b_g)]_j = sum_r 2^(r+1) b_{g,24r+j}`.  The self-test checks this
coordinate assignment directly in addition to encode/decode round trips.

This package measures complete KEM operations for LeechKEM-1/2/3 and compares
them on the same machine with FrodoKEM-640/976/1344-SHAKE.  Unlike the earlier
single-input harness, it rotates single-operation measurements through an
independently generated case pool and also reports a fresh end-to-end KEM cycle.

## Run in VS Code + WSL Ubuntu

Open the extracted folder in VS Code, choose **Terminal > New Terminal**, and
run:

```bash
chmod +x run_full_benchmark.sh
./run_full_benchmark.sh
```

Requirements are GCC/Make and Python 3 (the standard Ubuntu/WSL packages are
`build-essential` and `python3`).

The default paper run uses 31 independently generated cases, 31 timing samples,
and a minimum 20 ms batch.  A short installation check is:

```bash
./run_full_benchmark.sh 5 1
```

The default run may take several minutes.  Keep the laptop plugged in, select
the normal/high-performance power mode, and close heavy applications.  For
the paper, run the default command at least twice and use a run made under a
quiet, stable system state.

After the run, send back:

- `results/runtime_results.csv`
- `results/environment.txt`

Ready-to-paste tables are generated as `results/runtime_results.tex` and
`results/runtime_results.md`.

## What is included in each timing

| Operation | Included work |
|---|---|
| KeyGen | Seed generation, SHAKE expansion, error sampling, public-matrix generation, matrix multiplication, and public-key serialization |
| Encaps | Message/salt generation, `G`, error expansion, PKE encryption, Leech encoding, ciphertext serialization, and final `H` |
| Decaps (valid) | PKE decryption, Leech decoding, `G`, full deterministic re-encryption, two `H` calls, ciphertext comparison, and constant-time key selection |
| Decaps (invalid) | The same fixed-work decapsulation path using a one-bit-corrupted ciphertext, ending in the implicit-rejection key |
| Full KEM cycle | Fresh KeyGen, fresh Encaps, and valid Decaps in one timed invocation |

Thus this is the complete runtime requested for an implementation evaluation;
it is not merely the roughly 50--100 microsecond Leech decoder microbenchmark.

## Experimental choices

- LeechKEM is a portable, scalar proof-of-concept reference implementation,
  compiled with `-O3`.  It is intended for reproducible algorithm-level
  timing, not deployment or side-channel certification.
- The public matrix is generated with SHAKE128 from a 128-bit seed.  Secret
  expansion and `G`/`H` use SHAKE128 for LeechKEM-1 and SHAKE256 for
  LeechKEM-2/3.  Domain-separation bytes are used.
- The complete serialized public key is hashed, matching the current
  manuscript pseudocode.  Public keys and ciphertexts are bit-packed.
- The 36-bit Leech quotient representative uses the manuscript's exact
  coefficient-allocation table and normalized MOG generator matrix.  Its
  inverse is recovered from the Vardy--Be'ery decoder output by the documented
  coordinate isometry followed by triangular back substitution.
  The residual `g_tau` digits are then recovered coordinate by coordinate.
  Their bits are stored by bit plane, so coordinate `j` uses message positions
  `24*r+j`, exactly as in the manuscript formula.
- FrodoKEM uses Microsoft's official portable scalar reference code, the
  SHAKE public-matrix variant, and the same `-O3` optimization level.  AVX2 and
  AES-specific paths are disabled so the comparison is between scalar
  reference implementations.
- A deterministic benchmark-only RNG is used for both schemes so OS entropy
  acquisition is not included and both programs are reproducible.  This RNG
  must never be used in production cryptography.
- Before single-operation timing, the program prepares 31 independent
  `(pk, sk, ct, ss)` cases for each scheme (or the sample count supplied on the
  command line).  Valid and invalid Decaps rotate through this pool.  Case
  preparation and correctness checks are outside the Decaps timer; each
  invalid ciphertext is derived from its corresponding valid ciphertext by a
  one-bit change.
- `Full KEM cycle` is separate from the single-operation measurements.  Every
  timed invocation generates a fresh key pair and ciphertext and then performs
  valid decapsulation.  It is therefore a direct KeyGen--Encaps--Decaps total,
  rather than a sum of separately reported medians.
- The programs run functional self-tests before timing: a direct check of the
  manuscript's `g_tau` formula, zero, every basis bit, and 32 random noiseless
  Leech encode/decode round trips; valid KEM shared-key agreement; and
  invalid-ciphertext rejection.
- `tools/verify_mog_mapping.c` independently checks that all 24 rows of the
  MOG generator map to distance-zero points in the decoder lattice before a
  benchmark run starts.

## Corrected upstream utility-table entry

The decoder algorithm is vendored from `avanpo/leech-decoding` unchanged.
Its separate debugging utility table contained a transcription error in MOG
generator row 19.  That one row is corrected here to the standard row printed
in Figure 3.1 of the decoder author's thesis and in the manuscript:

```text
2 0 0 2 | 2 2 0 0 | 2 0 2 0 | 0 0 0 0 | 2 0 2 0 | 0 0 0 0
```

With that correction, the manuscript's adjacent-pair coordinate permutation
maps all 24 MOG generator rows to distance-zero decoder points.  The executable
certificate `build/verify_mog_mapping` performs this check before every run.

## Files

- `README_zh.md`: Chinese running and interpretation guide
- `src/leechkem_ref.c`: complete LeechKEM reference path
- `src/leechkem_bench.c`: LeechKEM timing harness and tests
- `src/frodo_bench.c`: official FrodoKEM wrapper and timing harness
- `tools/verify_mog_mapping.c`: 24-row MOG/decoder-coordinate certificate
- `paper/benchmark_method_en.tex`: English experimental-method paragraph
- `paper/benchmark_method_zh.tex`: Chinese counterpart
- `tools/merge_results.py`: CSV merger and Markdown/LaTeX table generator
- `THIRD_PARTY.md`: pinned source revisions and licenses

## Important interpretation limit

These measurements support a statement about the supplied scalar reference
implementations on the recorded platform.  They do not by themselves show
that an optimized LeechKEM implementation is faster or slower than every
optimized FrodoKEM implementation.  In the paper, describe the comparison as
a proof-of-concept/reference benchmark and report the compiler, CPU,
WSL/kernel, case-pool size, sample count, median, and p95 from
`results/environment.txt`.
