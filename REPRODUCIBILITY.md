# Reproducibility levels

This archive supports two checks.

## 1. Recreate all paper tables and comparison figures

This fast path reads the nine archived input JSON files and recreates every
CSV, audit JSON, PNG, and PDF in `results/`:

```bash
python -m pip install -r requirements.txt
python run.py
python tests.py
```

This is sufficient to reproduce the numerical results and figures reported in
the paper from the archived security inputs.

## 2. Recompute the nine security inputs

The `d_P` and `d_Q` calculations use ordinary Python plus `python-flint`. The
IND-CPA calculation additionally requires SageMath. Its exact Frodo cost model
and lattice-estimator source are bundled, so no personal path or Git checkout
is required.

Run the following commands from the archive root in a SageMath-enabled Python
environment:

```bash
python -m pip install -r requirements-recompute.txt
python recompute_inputs.py --stage dp
python recompute_inputs.py --stage dq
python recompute_inputs.py --stage lattice
python run.py --inputs data/regenerated --output results_regenerated
```

The three recomputation stages intentionally write to `data/regenerated/`.
They do not overwrite `data/raw/`, which is the archived reference used for
the paper. Runtime fields in estimator output can differ across machines; the
security values and final summaries are the quantities to compare.

To use a different compatible estimator checkout, pass
`--estimates-dir PATH`. The default is the bundled directory
`third_party/frodo_estimates/`.
