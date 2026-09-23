#!/usr/bin/env python3
"""Certified Chernoff upper bound for LeechPKE correctness under Frodo's Q.

This program computes an upper bound on

    d_Q = delta_Q(1),

the one-user correctness error of LeechPKE when every scalar secret/error
sample is drawn from Frodo's *rounded Gaussian* reference distribution

    Q = Psi_{sigma * sqrt(2*pi)} = Round(N(0, sigma^2)).

For one oriented Leech Voronoi facet represented by a = 4v,

    T_a = Z_a + sum_{j=1}^{2n} W_{a,j},
    Z_a = <X,a>,       W_a = R Z_a,

where all coordinates of X and the scalar R are independent Q samples.  For
every admissible theta > 0,

    Pr[T_a >= tau_a]
      <= exp(-theta*tau_a) M_Za(theta) M_Wa(theta)^(2n),
    tau_a = beta*||a||^2/8.

The seven Leech Voronoi classes are evaluated separately, multiplied by their
oriented multiplicities, summed, and finally multiplied by the number c of
message rows.

Numerical strategy
------------------
1. A fast finite-sum evaluation and bounded scalar minimization find a good
   candidate theta for every class.  This optimization is not trusted for the
   validity of the bound.
2. The fixed candidate is re-evaluated with Arb ball arithmetic
   (python-flint).  Both infinite sums, M_Q and M_W, receive analytic Gaussian
   tail upper bounds.  Any admissible fixed theta gives a valid Chernoff bound;
   failing to find the exact minimizer can only make the result looser.

The JSON output includes Arb enclosures, per-class contributions, truncation
diagnostics, and the final log2(d_Q) enclosure.

Dependencies
------------
    python -m pip install numpy scipy python-flint

Examples
--------
    python leech_dq_chernoff.py --preset L1
    python leech_dq_chernoff.py --preset all --dps 100 --tail-bits 120
    python leech_dq_chernoff.py --n 808 --beta 8192 --sigma 2.8 --rows 1
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
from scipy.special import log_ndtr, logsumexp, ndtr

try:
    from flint import arb, ctx
except ImportError as exc:  # pragma: no cover - exercised only without dependency
    raise SystemExit(
        "python-flint is required for the Arb verification step. Install with:\n"
        "  python -m pip install python-flint\n"
    ) from exc


@dataclass(frozen=True)
class Preset:
    name: str
    n: int
    beta: int
    sigma: float
    rows: int


PRESETS = {
    "L1": Preset("LeechKEM-1", n=880, beta=1 << 13, sigma=2.8, rows=1),
    "L2": Preset("LeechKEM-2", n=1136, beta=1 << 13, sigma=2.3, rows=2),
    "L3": Preset("LeechKEM-3", n=1640, beta=1 << 12, sigma=1.4, rows=2),
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


# a = 4v.  Multiplicities already contain all signs and coordinate permutations.
VORONOI_CLASSES = (
    VoronoiClass("V1", "(2^8,0^16)", (2,) * 8 + (0,) * 16, 97_152),
    VoronoiClass("V2", "(3,1^23)", (3,) + (1,) * 23, 98_304),
    VoronoiClass("V3", "(4^2,0^22)", (4,) * 2 + (0,) * 22, 1_104),
    VoronoiClass("V4", "(2^12,0^12)", (2,) * 12 + (0,) * 12, 5_275_648),
    VoronoiClass("V5", "(3^3,1^21)", (3,) * 3 + (1,) * 21, 8_290_304),
    VoronoiClass("V6", "(4,2^8,0^15)", (4,)  + (2,) * 8 + (0,) * 15, 3_108_864),
    VoronoiClass("V7", "(5,1^23)", (5,) + (1,) * 23, 98_304),
)

assert sum(v.multiplicity for v in VORONOI_CLASSES) == 16_969_680


def _logdiffexp(a: np.ndarray | float, b: np.ndarray | float):
    """Return log(exp(a)-exp(b)) for a >= b, stably."""
    a_arr = np.asarray(a)
    b_arr = np.asarray(b)
    return a_arr + np.log1p(-np.exp(b_arr - a_arr))


class FastRoundedGaussian:
    """Fast finite-grid evaluator used only to select candidate theta values."""

    def __init__(self, sigma: float, radius: int = 160):
        self.sigma = float(sigma)
        self.radius = int(radius)
        positive = np.arange(1, self.radius + 1, dtype=np.float64)
        z_lo = -(positive - 0.5) / self.sigma
        z_hi = -(positive + 0.5) / self.sigma
        log_positive = _logdiffexp(log_ndtr(z_lo), log_ndtr(z_hi))
        q0 = float(ndtr(0.5 / self.sigma) - ndtr(-0.5 / self.sigma))
        self.k = np.arange(-self.radius, self.radius + 1, dtype=np.float64)
        self.log_q = np.concatenate(
            [log_positive[::-1], np.array([math.log(q0)]), log_positive]
        )

    def log_mq(self, u: float) -> float:
        return float(logsumexp(self.log_q + abs(float(u)) * self.k))

    def log_mz(self, theta_r: float, vc: VoronoiClass) -> float:
        return sum(
            count * self.log_mq(coefficient * theta_r)
            for coefficient, count in vc.counts.items()
        )

    def log_mw(self, theta: float, vc: VoronoiClass) -> float:
        values = np.empty_like(self.k)
        for index, r in enumerate(self.k):
            values[index] = self.log_q[index] + self.log_mz(theta * r, vc)
        return float(logsumexp(values))

    def objective(
        self, theta: float, n: int, beta: int, vc: VoronoiClass
    ) -> float:
        threshold = beta * vc.norm2 / 8.0
        return (
            -theta * threshold
            + self.log_mz(theta, vc)
            + 2.0 * n * self.log_mw(theta, vc)
        )


def select_theta(
    model: FastRoundedGaussian,
    n: int,
    beta: int,
    sigma: float,
    vc: VoronoiClass,
) -> tuple[float, dict[str, float]]:
    """Find a numerical candidate in the proven convergence interval."""
    theta_limit = 1.0 / (sigma * sigma * math.sqrt(vc.norm2))
    lower = theta_limit * 1.0e-8
    upper = theta_limit * 0.97
    result = minimize_scalar(
        lambda theta: model.objective(theta, n, beta, vc),
        bounds=(lower, upper),
        method="bounded",
        options={"xatol": theta_limit * 1.0e-11, "maxiter": 300},
    )
    if not result.success or not (lower < result.x < upper):
        raise RuntimeError(f"theta optimization failed for {vc.name}: {result}")
    return float(result.x), {
        "theta_limit": theta_limit,
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
    """Serializable Arb enclosure plus convenient non-authoritative floats."""
    lower = float(value.lower())
    upper = float(value.upper())
    return {
        "arb": value.str(digits, radius=True),
        "lower_float": lower,
        "upper_float": upper,
    }


class ArbRoundedGaussian:
    """Arb evaluator with analytic upper bounds for both infinite tails."""

    def __init__(
        self,
        sigma: float,
        dps: int,
        tail_bits: int,
        max_inner_radius: int = 10000,
        max_outer_radius: int = 10000,
    ):
        ctx.dps = int(dps)
        self.sigma = _arb(sigma)
        self.sigma_float = float(sigma)
        self.sqrt2 = _arb(2).sqrt()
        self.half = _arb(1) / 2
        self.tail_ratio = _arb(2) ** (-int(tail_bits))
        self.tail_bits = int(tail_bits)
        self.max_inner_radius = int(max_inner_radius)
        self.max_outer_radius = int(max_outer_radius)
        self._q_cache: dict[int, arb] = {}
        self.max_inner_used = 0
        self.max_outer_used = 0
        self.max_inner_tail_ratio = arb(0)
        self.max_outer_tail_ratio = arb(0)

    def q(self, k: int) -> arb:
        """Exact rounded-Gaussian bin probability, enclosed by Arb."""
        k = abs(int(k))
        if k not in self._q_cache:
            if k == 0:
                x = self.half / (self.sigma * self.sqrt2)
                value = x.erf()
            else:
                lo = (_arb(k) - self.half) / (self.sigma * self.sqrt2)
                hi = (_arb(k) + self.half) / (self.sigma * self.sqrt2)
                value = (lo.erfc() - hi.erfc()) / 2
            if not value.is_finite() or value <= 0:
                raise ArithmeticError(f"failed to enclose q_{k}: {value}")
            self._q_cache[k] = value
        return self._q_cache[k]

    def normal_survival(self, z: arb) -> arb:
        return (z / self.sqrt2).erfc() / 2

    def mq_tail_upper(self, u: arb, radius: int) -> arb:
        """Upper-bound sum_{|k|>radius} q_k exp(u*k), for u >= 0."""
        b = _arb(radius) + self.half
        sigma2 = self.sigma * self.sigma
        shifted = sigma2 * u
        gaussian_mgf = (sigma2 * u * u / 2).exp()
        positive = self.normal_survival((b - shifted) / self.sigma)
        negative = self.normal_survival((b + shifted) / self.sigma)
        return (u / 2).exp() * gaussian_mgf * (positive + negative)

    def mq_upper(self, u: arb) -> tuple[arb, dict[str, Any]]:
        """Upper enclosure for M_Q(u), including the omitted infinite tail."""
        if u < 0:
            u = -u
        u_float = abs(float(u.mid()))
        # Center of the exponentially tilted continuous Gaussian plus a wide
        # initial safety window.  The loop, not this heuristic, is decisive.
        radius = max(
            8,
            int(
                math.ceil(
                    self.sigma_float**2 * u_float
                    + self.sigma_float
                    * math.sqrt(2.0 * (self.tail_bits + 8) * math.log(2.0))
                    + 1.0
                )
            ),
        )
        if radius > self.max_inner_radius:
            raise RuntimeError("initial M_Q radius exceeds configured maximum")

        while True:
            partial = self.q(0)
            for k in range(1, radius + 1):
                partial += 2 * self.q(k) * (_arb(k) * u).cosh()
            tail = self.mq_tail_upper(u, radius)
            ratio = tail / partial
            if ratio.upper() <= self.tail_ratio.lower():
                break
            radius += max(4, radius // 5)
            if radius > self.max_inner_radius:
                raise RuntimeError(
                    f"M_Q tail did not reach 2^-{self.tail_bits} by radius "
                    f"{self.max_inner_radius}"
                )

        self.max_inner_used = max(self.max_inner_used, radius)
        if ratio.upper() > self.max_inner_tail_ratio.upper():
            self.max_inner_tail_ratio = ratio
        return partial + tail, {
            "radius": radius,
            "tail_ratio": ratio,
        }

    def mz_upper(
        self, theta_r: arb, vc: VoronoiClass
    ) -> tuple[arb, list[dict[str, Any]]]:
        value = arb(1)
        diagnostics: list[dict[str, Any]] = []
        for coefficient, count in sorted(vc.counts.items()):
            mq, diag = self.mq_upper(abs(theta_r) * coefficient)
            value *= mq**count
            diagnostics.append(
                {
                    "coefficient": coefficient,
                    "count": count,
                    "inner_radius": diag["radius"],
                    "inner_tail_ratio_arb": diag["tail_ratio"].str(16, radius=True),
                }
            )
        return value, diagnostics

    def mw_tail_upper(
        self, theta: arb, vc: VoronoiClass, radius: int
    ) -> arb:
        """Bound E[exp(BR^2+A|R|); |R|>radius] for R=Round(G)."""
        sigma2 = self.sigma * self.sigma
        a_linear = theta * vc.l1 / 2
        b_quad = sigma2 * theta * theta * vc.norm2 / 2
        c_quad = 1 / (2 * sigma2) - b_quad
        if c_quad <= 0:
            raise ValueError("theta is outside the M_W convergence interval")
        d_linear = a_linear + b_quad
        constant = b_quad / 4 + a_linear / 2
        cutoff = _arb(radius) + self.half
        shifted_cutoff = c_quad.sqrt() * (
            cutoff - d_linear / (2 * c_quad)
        )
        prefactor = (
            constant + d_linear * d_linear / (4 * c_quad)
        ).exp() / (self.sigma * (2 * c_quad).sqrt())
        return prefactor * shifted_cutoff.erfc()

    def mw_upper(
        self, theta: arb, vc: VoronoiClass
    ) -> tuple[arb, dict[str, Any]]:
        sigma2_float = self.sigma_float**2
        theta_float = float(theta.mid())
        b_quad = sigma2_float * theta_float**2 * vc.norm2 / 2.0
        c_quad = 1.0 / (2.0 * sigma2_float) - b_quad
        a_linear = theta_float * vc.l1 / 2.0
        d_linear = a_linear + b_quad
        if c_quad <= 0:
            raise ValueError("theta is outside the M_W convergence interval")
        mode = max(0.0, d_linear / (2.0 * c_quad))
        radius = max(
            8,
            int(
                math.ceil(
                    mode
                    + math.sqrt((self.tail_bits + 8) * math.log(2.0) / c_quad)
                    + 1.0
                )
            ),
        )
        if radius > self.max_outer_radius:
            raise RuntimeError("initial M_W radius exceeds configured maximum")

        last_inner: list[dict[str, Any]] = []
        while True:
            partial = self.q(0)  # M_Z(0) = 1
            for r in range(1, radius + 1):
                mz, last_inner = self.mz_upper(theta * r, vc)
                partial += 2 * self.q(r) * mz
            tail = self.mw_tail_upper(theta, vc, radius)
            ratio = tail / partial
            if ratio.upper() <= self.tail_ratio.lower():
                break
            radius += max(4, radius // 5)
            if radius > self.max_outer_radius:
                raise RuntimeError(
                    f"M_W tail did not reach 2^-{self.tail_bits} by radius "
                    f"{self.max_outer_radius}"
                )

        self.max_outer_used = max(self.max_outer_used, radius)
        if ratio.upper() > self.max_outer_tail_ratio.upper():
            self.max_outer_tail_ratio = ratio
        return partial + tail, {
            "radius": radius,
            "tail_ratio": ratio,
            "last_inner": last_inner,
        }

    def certify_class(
        self,
        theta_decimal: str,
        n: int,
        beta: int,
        vc: VoronoiClass,
    ) -> dict[str, Any]:
        theta = arb(theta_decimal)
        theta_limit = 1 / (self.sigma * self.sigma * _arb(vc.norm2).sqrt())
        if theta <= 0 or theta >= theta_limit:
            raise ValueError(
                f"theta={theta_decimal} is outside (0, {theta_limit}) for {vc.name}"
            )
        threshold = _arb(beta * vc.norm2) / 8
        mz, mz_diag = self.mz_upper(theta, vc)
        mw, mw_diag = self.mw_upper(theta, vc)
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
            "theta_limit": _ball_json(theta_limit),
            "log2_single_facet_upper": _ball_json(log_single / log2),
            "log2_class_contribution_upper": _ball_json(
                contribution.log() / log2
            ),
            "single_facet_probability_upper": single,
            "class_contribution_upper": contribution,
            "mz_upper": _ball_json(mz),
            "mw_upper": _ball_json(mw),
            "mz_diagnostics": mz_diag,
            "mw_radius": mw_diag["radius"],
            "mw_tail_ratio_arb": mw_diag["tail_ratio"].str(20, radius=True),
        }


def certify_best_near_candidate(
    verifier: ArbRoundedGaussian,
    candidate: float,
    n: int,
    beta: int,
    sigma: float,
    vc: VoronoiClass,
) -> dict[str, Any]:
    """Certify a small deterministic neighborhood and retain the best bound."""
    theta_limit = 1.0 / (sigma * sigma * math.sqrt(vc.norm2))
    relative_offsets = (-2e-4, -5e-5, 0.0, 5e-5, 2e-4)
    trials: list[dict[str, Any]] = []
    for offset in relative_offsets:
        theta = candidate * (1.0 + offset)
        if not (0.0 < theta < theta_limit):
            continue
        # Seventeen significant digits preserve the selected binary64 value as
        # an exact decimal input to Arb.
        theta_decimal = format(theta, ".17g")
        row = verifier.certify_class(theta_decimal, n, beta, vc)
        trials.append(row)
    if not trials:
        raise RuntimeError(f"no admissible theta trial for {vc.name}")
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


def compute_preset(
    preset: Preset,
    dps: int,
    tail_bits: int,
    fast_radius: int,
) -> dict[str, Any]:
    started = time.time()
    fast = FastRoundedGaussian(preset.sigma, radius=fast_radius)
    verifier = ArbRoundedGaussian(
        preset.sigma,
        dps=dps,
        tail_bits=tail_bits,
    )
    rows: list[dict[str, Any]] = []
    print(
        f"{preset.name}: n={preset.n}, beta={preset.beta}, "
        f"sigma={preset.sigma}, rows={preset.rows}",
        flush=True,
    )
    for vc in VORONOI_CLASSES:
        candidate, fast_info = select_theta(
            fast, preset.n, preset.beta, preset.sigma, vc
        )
        print(
            f"  {vc.name}: candidate theta={candidate:.12g}, "
            f"fast single-facet log2={fast_info['fast_log2_single_facet']:.6f}",
            flush=True,
        )
        row = certify_best_near_candidate(
            verifier,
            candidate,
            preset.n,
            preset.beta,
            preset.sigma,
            vc,
        )
        row["fast_optimization"] = fast_info
        print(
            f"      certified class log2 upper <= "
            f"{row['log2_class_contribution_upper']['upper_float']:.6f}, "
            f"outer radius={row['mw_radius']}",
            flush=True,
        )
        rows.append(row)

    total = arb(0)
    for row in rows:
        total += row.pop("class_contribution_upper")
        row.pop("single_facet_probability_upper")
    total *= preset.rows
    if total > 1:
        total = arb(1)
    total_log2 = total.log() / _arb(2).log()
    result = {
        "schema": "leech-dq-rounded-gaussian-chernoff-v1",
        "definition": {
            "Q": "Round(N(0,sigma^2)); Frodo Psi_{sigma*sqrt(2*pi)}",
            "quantity": "d_Q = delta_Q(1)",
            "method": "seven-class Chernoff union bound with Arb and analytic tails",
            "formal_scope": (
                "Arb encloses fixed-theta arithmetic and special functions; "
                "theta optimization only selects candidates and is not needed "
                "for validity"
            ),
        },
        "parameters": asdict(preset),
        "verification": {
            "arb_decimal_digits": dps,
            "requested_relative_tail_bits": tail_bits,
            "fast_grid_radius": fast_radius,
            "maximum_inner_radius_used": verifier.max_inner_used,
            "maximum_outer_radius_used": verifier.max_outer_used,
            "maximum_inner_tail_ratio_arb": verifier.max_inner_tail_ratio.str(
                20, radius=True
            ),
            "maximum_outer_tail_ratio_arb": verifier.max_outer_tail_ratio.str(
                20, radius=True
            ),
        },
        "classes": rows,
        "log2_dQ_upper": _ball_json(total_log2),
        "dQ_upper": _ball_json(total),
        "elapsed_seconds": time.time() - started,
    }
    print(
        f"  total log2(d_Q) upper <= "
        f"{result['log2_dQ_upper']['upper_float']:.9f}",
        flush=True,
    )
    return result


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--preset",
        choices=("L1", "L2", "L3", "all", "custom"),
        default="L1",
        help="parameter set; custom requires --n, --beta, --sigma and --rows",
    )
    parser.add_argument("--n", type=int)
    parser.add_argument("--beta", type=int)
    parser.add_argument("--sigma", type=float)
    parser.add_argument("--rows", type=int)
    parser.add_argument("--name", default="LeechKEM-custom")
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--tail-bits", type=int, default=120)
    parser.add_argument("--fast-radius", type=int, default=160)
    parser.add_argument(
        "--output",
        type=Path,
        help="JSON output path; for --preset all, this is treated as a directory",
    )
    args = parser.parse_args(argv)
    if args.dps < 40:
        parser.error("--dps must be at least 40")
    if args.tail_bits < 30:
        parser.error("--tail-bits must be at least 30")
    if args.fast_radius < 40:
        parser.error("--fast-radius must be at least 40")
    if args.preset == "custom":
        if any(value is None for value in (args.n, args.beta, args.sigma, args.rows)):
            parser.error("custom mode requires --n, --beta, --sigma and --rows")
        if args.n <= 0 or args.beta <= 0 or args.sigma <= 0 or args.rows <= 0:
            parser.error("custom parameters must be positive")
    elif any(value is not None for value in (args.n, args.beta, args.sigma, args.rows)):
        parser.error("--n/--beta/--sigma/--rows are only valid with --preset custom")
    return args


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)
    if args.preset == "all":
        presets = list(PRESETS.values())
        output_dir = args.output or Path("leech_dq_results")
        output_dir.mkdir(parents=True, exist_ok=True)
    elif args.preset == "custom":
        presets = [
            Preset(
                args.name,
                n=args.n,
                beta=args.beta,
                sigma=args.sigma,
                rows=args.rows,
            )
        ]
        output_dir = None
    else:
        presets = [PRESETS[args.preset]]
        output_dir = None

    for preset in presets:
        result = compute_preset(
            preset,
            dps=args.dps,
            tail_bits=args.tail_bits,
            fast_radius=args.fast_radius,
        )
        if output_dir is not None:
            path = output_dir / f"{preset.name.lower().replace('-', '_')}_dq.json"
        elif args.output is not None:
            path = args.output
        else:
            slug = preset.name.lower().replace("-", "_")
            path = Path(f"{slug}_dq_chernoff.json")
        path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
        print(f"Wrote {path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
