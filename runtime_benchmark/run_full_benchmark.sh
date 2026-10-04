#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$project_dir"

samples="${1:-31}"
minimum_ms="${2:-20}"

if ! [[ "$samples" =~ ^[0-9]+$ ]] || (( samples < 5 || samples % 2 == 0 )); then
    echo "Error: samples must be an odd integer greater than or equal to 5." >&2
    exit 2
fi
if ! [[ "$minimum_ms" =~ ^[0-9]+([.][0-9]+)?$ ]] || [[ "$minimum_ms" == "0" ]]; then
    echo "Error: min-ms must be a positive number." >&2
    exit 2
fi

mkdir -p build results

echo "[1/8] Compiling scalar reference implementations..."
make all

echo "[2/8] Verifying the 24 MOG generator rows in decoder coordinates..."
./build/verify_mog_mapping

echo "[3/8] Running LeechKEM-1/2/3 with independent case pools..."
./build/leechkem_bench --samples "$samples" --min-ms "$minimum_ms" \
    --output results/leech_full.csv

echo "[4/8] Running FrodoKEM-640-SHAKE with an independent case pool..."
./build/frodo640_bench --samples "$samples" --min-ms "$minimum_ms" \
    --output results/frodo640_full.csv

echo "[5/8] Running FrodoKEM-976-SHAKE with an independent case pool..."
./build/frodo976_bench --samples "$samples" --min-ms "$minimum_ms" \
    --output results/frodo976_full.csv

echo "[6/8] Running FrodoKEM-1344-SHAKE with an independent case pool..."
./build/frodo1344_bench --samples "$samples" --min-ms "$minimum_ms" \
    --output results/frodo1344_full.csv

echo "[7/8] Recording the experiment environment..."
{
    echo "timestamp_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "kernel=$(uname -a)"
    echo "cpu=$(awk -F: '/model name/{gsub(/^[[:space:]]+/, "", $2); print $2; exit}' /proc/cpuinfo)"
    echo "logical_cpus=$(getconf _NPROCESSORS_ONLN)"
    echo "compiler=$(${CC:-gcc} --version | head -n 1)"
    echo "leech_cflags=-O3 -std=c11 -Wall -Wextra -Wpedantic"
    echo "frodo_cflags=-O3 -std=gnu11 -Wall -Wextra -Wpedantic -D_REFERENCE_ -D_SHAKE128_FOR_A_"
    echo "samples=$samples"
    echo "independent_case_pool_size=$samples"
    echo "minimum_batch_ms=$minimum_ms"
    echo "full_cycle=KeyGen + Encaps + valid Decaps (fresh randomness per iteration)"
    echo "rng=deterministic benchmark stream (OS entropy acquisition excluded)"
} > results/environment.txt

echo "[8/8] Producing CSV, Markdown, and LaTeX tables..."
python3 tools/merge_results.py

echo
echo "Done. Please send back these two files:"
echo "  results/runtime_results.csv"
echo "  results/environment.txt"
echo "A ready-to-paste table is in results/runtime_results.tex"
