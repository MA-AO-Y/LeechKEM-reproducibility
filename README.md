# LeechKEM Reproducibility Package

This repository contains the fixed-parameter numerical reproduction package
for the paper *LeechKEM: A Plain-LWE Key Encapsulation Mechanism with
Leech-Lattice Quotient Decoding*. It reproduces the conditional single-user,
multi-challenge classical-ROM estimates and comparison figures for the three
fixed LeechKEM parameter sets. It does not search over the LWE dimension `n`
or provide a production implementation of the KEM. The fixed dimensions are
880, 1136, and 1640.

## Quick reproduction

Python 3.10 or later is sufficient for the quick path. Sage is not required.
From this directory, run:

```bash
python -m pip install -r requirements.txt
python run.py
python tests.py
```

Using the existing Sage conda environment is also fine:

```bash
conda activate sage
python -m pip install -r requirements.txt
python run.py
python tests.py
```

The expected design-point output is:

```text
L1: n=880,  N=128, 131.034770695 bits, margin=3.034770695
L2: n=1136, N=128, 194.495437386 bits, margin=2.495437386
L3: n=1640, N=64,  260.654644874 bits, margin=4.654644874
```

`tests.py` must finish with `OK`.

The generated files are written to `results/`:

- `L1_comparison`, `L2_comparison`, and `L3_comparison` in PNG and vector PDF;
- one CSV per parameter set, containing every integer `N` from 1 through 128;
- one audit JSON per set, containing the input hashes and every bound term;
- `summary.json`, containing the three design-point results.

Every integer challenge count is evaluated. Plot markers are shown only at
selected values to keep the figures readable. The LeechKEM-1 and LeechKEM-2
design points are `N=128`; the LeechKEM-3 design point is `N=64`.

## Package layout

- `security.py`: fixed parameters and log-domain reduction calculations;
- `run.py`: input validation, evaluation, CSV/JSON output, and plotting;
- `data/raw/`: the nine archived IND-CPA, `d_P`, and `d_Q` JSON files used for
  the paper;
- `recompute_inputs.py`: optional fixed-parameter input recomputation;
- `generators/`: the correctness-bound implementations and lattice-estimator
  wrapper;
- `third_party/frodo_estimates/`: the fixed Frodo cost model and complete
  lattice-estimator source needed by the IND-CPA recomputation;
- `tests.py`: design-point regression and an independent ordinary-domain check
  of the Frodo formula.

See `REPRODUCIBILITY.md` for the two reproduction levels and `THIRD_PARTY.md`
for source commits and licenses.

## Security convention and formulas

All curves use the normalized distinguishing advantage

```text
Adv = 2 * |Pr[win] - 1/2|.
```

The outputs are conditional concrete reduction estimates in the classical ROM
under the stated core-SVP/C-LSF assumptions. They are not measured attack-cost
curves. The calculation fixes

```text
t = 2^32,        qG = qH = t / 2^18 = 2^14.
```

The IND-CPA input is interpreted as a normalized per-gate advantage bound.
For `N` challenges, its contribution is bounded by
`N*t*2^(-bCPA)`. The GHS salt-multiplicity parameter is `T=1`, and the salt
collision bad event `N*(N-1)/(2*2^salt)` is included separately.

The LeechKEM message lengths are 132, 216, and 264 bits. The final shared-secret
lengths are 128, 192, and 256 bits. These are distinct parameters. The package
uses the fixed direct-H transform:

```text
H(pk || mu || ct0 || salt)
```

and replaces `mu` by the implicit-rejection value `z` on decapsulation failure.

Let `M=2^mu`, `S=2^salt`, `p=1-1/alpha`, and
`R=s_N*D_alpha(P||Q)/ln(2)`. The LeechKEM calculation is

```text
A = qG*dQ + N/M + 2*qG/M + 2*N*t*2^(-bCPA)
B = min(1, 2*((A*2^R)^p + dP + qH/M
              + N*(N-1)/(M*S) + N*(N-1)/(2*S)))
bits = log2(t) - log2(B)
s_N = 2*n*24 + N*(2*rows*n + rows*24)
```

The outer factor two converts the GHS unnormalized IND advantage to the
normalized convention used by Frodo. Likewise, the `2*AdvCPA` term above is
the normalized form of the `4*AdvCPA` term in GHS Theorem 7.

For Frodo, the code evaluates Theorem 1, Equation (3), directly at the same
`t` and query budget. The inner terms are

```text
6*q/M + 2^(-512) + N/M + q*delta + 2*N*t*2^(-bCPA),
```

and the outer terms are

```text
delta + q/M + N*(N-1)/(M*S).
```

The total bound is multiplied by two, and Frodo uses `nbar=mbar=8` in its
sample count. Its `bCPA`, `delta`, `alpha`, and Renyi-divergence values are
literature inputs; this package does not re-prove them.

References:

- GHS multi-target analysis: <https://cic.iacr.org/p/3/1/20/pdf>
- FrodoKEM specification and analysis: <https://cic.iacr.org/p/2/3/25/pdf>

## Independent numerical verification

The final design-point results were independently recomputed from the nine raw
JSON files with 70-digit decimal arithmetic, without calling `security.py`:

```text
L1  131.0347706951662980513884273086928405265882
L2  194.4954373857996936708442221488567837089557
L3  260.6546448735815944668066877081217364037309
```

They agree with `run.py` to approximately `1e-14` bit. Each audit JSON records
SHA-256 hashes of the raw inputs, and `run.py` verifies the dimensions, modulus,
message length, beta, and C-LSF hybrid loss before evaluation.

The final computation uses stable binary64 log-domain arithmetic. The original
`d_P` and `d_Q` JSON files retain their Arb enclosures, but the whole final
calculation is not an end-to-end interval-arithmetic certificate.

## Optional recomputation of all nine inputs

The quick path reads the supplied raw JSON files. To recompute the correctness
inputs for the same fixed parameter sets, use:

```bash
python -m pip install -r requirements-recompute.txt
python recompute_inputs.py --stage dp
python recompute_inputs.py --stage dq
```

The IND-CPA lattice stage must be run with a SageMath-enabled Python. The fixed
Frodo cost model and lattice-estimator source are already included in this
archive, so no personal path or separate download is needed:

```bash
conda activate sage
python recompute_inputs.py --stage lattice
```

After all three stages have populated `data/regenerated/`, evaluate them with:

```bash
python run.py \
  --inputs data/regenerated \
  --output results_regenerated
```

Use `--set L1`, `--set L2`, or `--set L3` to recompute only one set. None of
these commands searches over `n`. The correctness stages, especially `d_Q`,
can take much longer than `run.py`; this is expected. The bundled source
versions are recorded in `THIRD_PARTY.md`.

## Citation

A version-specific citation and DOI will be added here when the first Zenodo
release is published.
