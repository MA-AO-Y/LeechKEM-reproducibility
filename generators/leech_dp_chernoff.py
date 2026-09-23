#!/usr/bin/env python3
"""Certified Chernoff upper bounds for LeechPKE correctness under P=chi.

The GHS 2026 modular proof needs two different correctness quantities:

    d_P : correctness error of the implemented finite Frodo distribution P=chi;
    d_Q : correctness error of the ideal rounded-Gaussian reference Q.

This program computes d_P.  Its companion ``leech_dq_chernoff.py`` computes
d_Q.  For one oriented Leech Voronoi facet represented by ``a=4v``, write

    T_a = Z_a + sum_{j=1}^{2n} W_{a,j},
    Z_a = <X,a>,                 W_a = R Z_a,

where the coordinates of X and R are independent samples from P.  For every
theta>0,

    Pr[T_a >= tau_a]
      <= exp(-theta*tau_a) M_Za(theta) M_Wa(theta)^(2n),
    tau_a = beta*||a||^2/8.

The program selects a good theta with binary64/SciPy, then re-evaluates that
*fixed* theta using Arb ball arithmetic.  The optimization itself is not part
of validity: every positive fixed theta gives a valid Chernoff upper bound.
The seven oriented Leech Voronoi classes are evaluated separately, multiplied
by their exact multiplicities, summed, and multiplied by the number of message
rows.

The finite distributions are reconstructed exactly from FrodoKEM's 15-bit
magnitude CDF tables.  The final presets are:

    L1: n=824,  beta=2^13, chi_Frodo-640,  rows=1
    L2: n=1072, beta=2^13, chi_Frodo-976,  rows=2
    L3: n=1504, beta=2^12, chi_Frodo-1344, rows=2

Dependencies:
    python -m pip install numpy scipy python-flint

Examples:
    python leech_dp_chernoff.py --preset L1
    python leech_dp_chernoff.py --preset all --dps 100
"""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable
import argparse
import json
import math
import sys
import time

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

try:
    from flint import arb, ctx
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "python-flint is required for the Arb verification step. Install with:\n"
        "  python -m pip install python-flint\n"
    ) from exc


@dataclass(frozen=True)
class NoiseSpec:
    name: str
    nominal_sigma: float
    support: int
    cdf_table: tuple[int, ...]


NOISES = {
    "frodo640": NoiseSpec(
        "chi_Frodo-640",
        2.8,
        12,
        (4643, 13363, 20579, 25843, 29227, 31145, 32103,
         32525, 32689, 32745, 32762, 32766, 32767),
    ),
    "frodo976": NoiseSpec(
        "chi_Frodo-976",
        2.3,
        10,
        (5638, 15915, 23689, 28571, 31116, 32217,
         32613, 32731, 32760, 32766, 32767),
    ),
    "frodo1344": NoiseSpec(
        "chi_Frodo-1344",
        1.4,
        6,
        (9142, 23462, 30338, 32361, 32725, 32765, 32767),
    ),
}


@dataclass(frozen=True)
class Preset:
    key: str
    name: str
    n: int
    beta: int
    noise_key: str
    rows: int


PRESETS = {
    "L1": Preset("L1", "LeechKEM-1", 880, 1 << 13, "frodo640", 1),
    "L2": Preset("L2", "LeechKEM-2", 1136, 1 << 13, "frodo976", 2),
    "L3": Preset("L3", "LeechKEM-3", 1640, 1 << 12, "frodo1344", 2),
}


@dataclass(frozen=True)
class VoronoiClass:
    name: str
    representative: str
    coefficients: tuple[int, ...]
    multiplicity: int

    @property
    def counts(self) -> Counter[int]:
        return Counter(abs(x) for x in self.coefficients if x)

    @property
    def norm2(self) -> int:
        return sum(x * x for x in self.coefficients)

    @property
    def l1(self) -> int:
        return sum(abs(x) for x in self.coefficients)


VORONOI_CLASSES = (
    VoronoiClass("V1", "(2^8,0^16)", (2,) * 8 + (0,) * 16, 97_152),
    VoronoiClass("V2", "(3,1^23)", (3,) + (1,) * 23, 98_304),
    VoronoiClass("V3", "(4^2,0^22)", (4,) * 2 + (0,) * 22, 1_104),
    VoronoiClass("V4", "(2^12,0^12)", (2,) * 12 + (0,) * 12, 5_275_648),
    VoronoiClass("V5", "(3^3,1^21)", (3,) * 3 + (1,) * 21, 8_290_304),
    VoronoiClass(
        "V6", "(4,2^8,0^15)",
        (4,) + (2,) * 8 + (0,) * 15, 3_108_864,
    ),
    VoronoiClass("V7", "(5,1^23)", (5,) + (1,) * 23, 98_304),
)

assert sum(v.multiplicity for v in VORONOI_CLASSES) == 16_969_680


def magnitude_counts(noise: NoiseSpec) -> tuple[int, ...]:
    """Return exact counts for |X|=0,...,support out of 2^15."""
    cdf = noise.cdf_table
    if len(cdf) != noise.support + 1:
        raise ValueError(f"{noise.name}: invalid CDF length")
    counts = [cdf[0] + 1]
    counts.extend(cdf[i] - cdf[i - 1] for i in range(1, noise.support))
    counts.append(32767 - cdf[noise.support - 1])
    if sum(counts) != 32768 or min(counts) <= 0:
        raise ValueError(f"{noise.name}: invalid magnitude counts")
    return tuple(counts)


def signed_float_pmf(noise: NoiseSpec) -> tuple[np.ndarray, np.ndarray]:
    counts = magnitude_counts(noise)
    values = np.arange(-noise.support, noise.support + 1, dtype=np.float64)
    probabilities = np.zeros(values.size, dtype=np.float64)
    probabilities[noise.support] = counts[0] / 32768.0
    for k in range(1, noise.support + 1):
        probability_one_sign = counts[k] / 65536.0
        probabilities[noise.support - k] = probability_one_sign
        probabilities[noise.support + k] = probability_one_sign
    if not np.isclose(probabilities.sum(), 1.0, atol=1e-15):
        raise ArithmeticError(f"{noise.name}: PMF is not normalized")
    return values, probabilities


class FastFiniteDistribution:
    """Binary64 evaluator used only to select candidate theta values."""

    def __init__(self, noise: NoiseSpec):
        self.noise = noise
        self.values, self.probabilities = signed_float_pmf(noise)
        self.log_probabilities = np.log(self.probabilities)

    def log_m(self, u: float) -> float:
        return float(logsumexp(self.log_probabilities + u * self.values))

    def log_mz(self, theta_r: float, vc: VoronoiClass) -> float:
        return sum(
            count * self.log_m(coefficient * theta_r)
            for coefficient, count in vc.counts.items()
        )

    def log_mw(self, theta: float, vc: VoronoiClass) -> float:
        return float(
            logsumexp(
                self.log_probabilities
                + np.array(
                    [self.log_mz(theta * r, vc) for r in self.values],
                    dtype=np.float64,
                )
            )
        )

    def objective(self, theta: float, n: int, beta: int, vc: VoronoiClass) -> float:
        threshold = beta * vc.norm2 / 8.0
        return (
            -theta * threshold
            + self.log_mz(theta, vc)
            + 2.0 * n * self.log_mw(theta, vc)
        )


def finite_support_maximum(noise: NoiseSpec, n: int, vc: VoronoiClass) -> int:
    return noise.support * vc.l1 * (1 + 2 * n * noise.support)


def select_theta(
    model: FastFiniteDistribution,
    n: int,
    beta: int,
    vc: VoronoiClass,
) -> tuple[float, dict[str, float]]:
    """Select a good positive theta; validity never depends on optimality."""
    high = 0.25
    previous = model.objective(high / 2.0, n, beta, vc)
    current = model.objective(high, n, beta, vc)
    while current < previous and high < 8.0:
        high *= 2.0
        previous, current = current, model.objective(high, n, beta, vc)
    result = minimize_scalar(
        lambda theta: model.objective(theta, n, beta, vc),
        bounds=(1.0e-12, high),
        method="bounded",
        options={"xatol": 1.0e-13, "maxiter": 1000},
    )
    if not result.success or result.x <= 0:
        raise RuntimeError(f"theta optimization failed for {vc.name}: {result}")
    return float(result.x), {
        "search_upper": high,
        "fast_objective": float(result.fun),
        "fast_log2_single_facet": float(result.fun / math.log(2.0)),
    }


def _arb(value: int | float | str) -> arb:
    if isinstance(value, float):
        return arb(repr(value))
    return arb(value)


def _upper_float(value: arb) -> float:
    return float(value.upper())


def _ball_json(value: arb, digits: int = 30) -> dict[str, Any]:
    return {
        "arb": value.str(digits, radius=True),
        "lower_float": float(value.lower()),
        "upper_float": float(value.upper()),
    }


class ArbFiniteDistribution:
    """Exact-rational Frodo PMF with Arb-enclosed exponential arithmetic."""

    def __init__(self, noise: NoiseSpec, dps: int):
        ctx.dps = int(dps)
        self.noise = noise
        counts = magnitude_counts(noise)
        self.p0 = _arb(counts[0]) / 32768
        self.p_positive = {
            k: _arb(counts[k]) / 65536
            for k in range(1, noise.support + 1)
        }
        total = self.p0 + 2 * sum(self.p_positive.values(), arb(0))
        if not total.contains(1):
            raise ArithmeticError(f"{noise.name}: Arb PMF does not contain one")

    def m(self, u: arb) -> arb:
        if u < 0:
            u = -u
        value = self.p0
        for k, probability in self.p_positive.items():
            value += 2 * probability * (_arb(k) * u).cosh()
        return value

    def mz(self, theta_r: arb, vc: VoronoiClass) -> arb:
        value = arb(1)
        for coefficient, count in sorted(vc.counts.items()):
            value *= self.m(abs(theta_r) * coefficient) ** count
        return value

    def mw(self, theta: arb, vc: VoronoiClass) -> arb:
        value = self.p0  # M_Z(0)=1
        for r, probability in self.p_positive.items():
            value += 2 * probability * self.mz(theta * r, vc)
        return value

    def certify_class(
        self,
        theta_decimal: str,
        n: int,
        beta: int,
        vc: VoronoiClass,
    ) -> dict[str, Any]:
        theta = arb(theta_decimal)
        if theta <= 0:
            raise ValueError("theta must be positive")
        threshold = _arb(beta * vc.norm2) / 8
        mz = self.mz(theta, vc)
        mw = self.mw(theta, vc)
        log_single = -theta * threshold + mz.log() + 2 * n * mw.log()
        single = log_single.exp()
        contribution = vc.multiplicity * single
        log2 = _arb(2).log()
        return {
            "name": vc.name,
            "representative_a": vc.representative,
            "multiplicity": vc.multiplicity,
            "norm2": vc.norm2,
            "l1": vc.l1,
            "threshold": beta * vc.norm2 // 8,
            "theta": theta_decimal,
            "log2_single_facet_upper": _ball_json(log_single / log2),
            "log2_class_contribution_upper": _ball_json(contribution.log() / log2),
            "mz_upper": _ball_json(mz),
            "mw_upper": _ball_json(mw),
            "single_facet_probability_upper": single,
            "class_contribution_upper": contribution,
        }


def certify_best_near_candidate(
    verifier: ArbFiniteDistribution,
    candidate: float,
    n: int,
    beta: int,
    vc: VoronoiClass,
) -> dict[str, Any]:
    trials: list[dict[str, Any]] = []
    for relative_offset in (-2e-4, -5e-5, 0.0, 5e-5, 2e-4):
        theta = candidate * (1.0 + relative_offset)
        if theta <= 0:
            continue
        row = verifier.certify_class(format(theta, ".17g"), n, beta, vc)
        trials.append(row)
    if not trials:
        raise RuntimeError(f"no positive theta trial for {vc.name}")
    best = min(
        trials,
        key=lambda row: _upper_float(row["class_contribution_upper"]),
    )
    best["theta_trials"] = [
        {
            "theta": row["theta"],
            "log2_class_contribution_upper": row[
                "log2_class_contribution_upper"
            ],
        }
        for row in trials
    ]
    return best


def compute_preset(preset: Preset, dps: int) -> dict[str, Any]:
    started = time.time()
    noise = NOISES[preset.noise_key]
    fast = FastFiniteDistribution(noise)
    verifier = ArbFiniteDistribution(noise, dps)
    rows: list[dict[str, Any]] = []
    print(
        f"{preset.name}: n={preset.n}, beta={preset.beta}, "
        f"P={noise.name}, rows={preset.rows}",
        flush=True,
    )
    for vc in VORONOI_CLASSES:
        threshold = preset.beta * vc.norm2 // 8
        support_max = finite_support_maximum(noise, preset.n, vc)
        if support_max < threshold:
            rows.append(
                {
                    "name": vc.name,
                    "representative_a": vc.representative,
                    "multiplicity": vc.multiplicity,
                    "norm2": vc.norm2,
                    "l1": vc.l1,
                    "threshold": threshold,
                    "support_max": support_max,
                    "method": "finite-support-zero",
                    "log2_single_facet_upper": None,
                    "log2_class_contribution_upper": None,
                    "class_contribution_upper": arb(0),
                }
            )
            print(f"  {vc.name}: probability is exactly zero by support", flush=True)
            continue

        candidate, fast_info = select_theta(
            fast, preset.n, preset.beta, vc
        )
        print(
            f"  {vc.name}: candidate theta={candidate:.12g}, "
            f"fast single-facet log2={fast_info['fast_log2_single_facet']:.6f}",
            flush=True,
        )
        row = certify_best_near_candidate(
            verifier, candidate, preset.n, preset.beta, vc
        )
        row["support_max"] = support_max
        row["method"] = "Chernoff-Arb-fixed-theta"
        row["fast_optimization"] = fast_info
        print(
            f"      certified class log2 upper <= "
            f"{row['log2_class_contribution_upper']['upper_float']:.6f}",
            flush=True,
        )
        rows.append(row)

    total = arb(0)
    for row in rows:
        total += row.pop("class_contribution_upper")
        row.pop("single_facet_probability_upper", None)
    total *= preset.rows
    if total > 1:
        total = arb(1)
    total_log2 = total.log() / _arb(2).log() if total > 0 else None

    values, probabilities = signed_float_pmf(noise)
    realized_sigma = math.sqrt(float(np.sum(probabilities * values * values)))
    result = {
        "schema": "leech-dp-finite-chernoff-v1",
        "definition": {
            "P": noise.name,
            "quantity": "d_P = delta_P(1)",
            "method": "seven-class Chernoff union bound with fixed-theta Arb verification",
            "formal_scope": (
                "The Frodo CDF probabilities are exact rationals. Arb encloses "
                "the fixed-theta exponential arithmetic; candidate optimization "
                "is not required for validity."
            ),
        },
        "parameters": asdict(preset),
        "noise": {
            **asdict(noise),
            "magnitude_counts_out_of_2^15": magnitude_counts(noise),
            "realized_sigma_binary64": realized_sigma,
        },
        "verification": {"arb_decimal_digits": dps},
        "classes": rows,
        "log2_dP_upper": (
            {"arb": "-inf", "lower_float": -math.inf, "upper_float": -math.inf}
            if total_log2 is None
            else _ball_json(total_log2)
        ),
        "dP_upper": _ball_json(total),
        "elapsed_seconds": time.time() - started,
    }
    print(
        f"  total log2(d_P) upper <= "
        f"{result['log2_dP_upper']['upper_float']:.9f}",
        flush=True,
    )
    return result


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--preset", choices=("L1", "L2", "L3", "all"), default="L1"
    )
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument(
        "--output", type=Path,
        help="JSON path; with --preset all, this is an output directory",
    )
    args = parser.parse_args(argv)
    if args.dps < 40:
        parser.error("--dps must be at least 40")
    return args


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    if args.preset == "all":
        presets = list(PRESETS.values())
        output_dir = args.output or Path("leech_dp_results")
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        presets = [PRESETS[args.preset]]
        output_dir = None

    for preset in presets:
        result = compute_preset(preset, args.dps)
        if output_dir is not None:
            path = output_dir / f"{preset.name.lower().replace('-', '_')}_dp.json"
        elif args.output is not None:
            path = args.output
        else:
            slug = preset.name.lower().replace("-", "_")
            path = Path(f"{slug}_dp_chernoff.json")
        path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
        print(f"Wrote {path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
