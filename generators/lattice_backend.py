from typing import Any
import time

MODEL_SPECS = (
    ("Q-RW-Sieve", "qrw", r"\qrandwalksieve"),
    ("Q-Grover-Sieve", "qgrover", r"\qgroversieve"),
    ("C-LSF-Sieve", "clsf", r"\clsfsieve"),
)

ATTACK_SPECS = (
    ("usvp", "Primal uSVP"),
    ("bdd", "Primal BDD"),
    ("dual", "Basic dual"),
    ("dual_hybrid", "Dual-sieve-FFT"),
)

DENY_LIST = ("arora-gb", "bkw", "bdd_hybrid", "bdd_mitm_hybrid")
_ESTIMATOR_COMPONENTS = None


def load_estimator_components():
    global _ESTIMATOR_COMPONENTS
    if _ESTIMATOR_COMPONENTS is not None:
        return _ESTIMATOR_COMPONENTS
    try:
        from sage.all import log as sage_log
        from lattice_estimator.estimator import LWE, Simulator
        from lattice_estimator.estimator.nd import DiscreteGaussian
        from lattice_estimator.estimator.lwe_parameters import LWEParameters
        from cost_models import reduction_algs, svp_models
    except Exception as exc:
        raise RuntimeError(
            "Sage/Frodo estimator imports failed. Run this stage with a "
            "SageMath-enabled Python. The fixed Frodo cost model and "
            "lattice-estimator source are bundled under third_party/."
        ) from exc
    _ESTIMATOR_COMPONENTS = (
        sage_log,
        LWE,
        Simulator,
        DiscreteGaussian,
        LWEParameters,
        reduction_algs,
        svp_models,
    )
    return _ESTIMATOR_COMPONENTS


def _as_int_or_none(value):
    try:
        return int(value)
    except Exception:
        return None


def _parse_attack(key: str, item, n: int, sage_log) -> dict[str, Any]:
    if item is None or isinstance(item, Exception):
        error = "missing" if item is None else repr(item)
        raise RuntimeError(f"estimator attack {key!r} failed: {error}")
    try:
        log2_cost = float(sage_log(item["rop"], 2))
    except Exception as exc:
        raise RuntimeError(f"cannot read estimator cost for {key!r}") from exc
    output: dict[str, Any] = {"log2_cost": log2_cost}
    for field in ("beta", "eta", "d", "beta_", "m", "zeta", "t"):
        if field in item:
            output[field] = _as_int_or_none(item[field])
    if key == "dual_hybrid":
        try:
            output["effective_dimension"] = int(
                item["m"] + n - item["zeta"] - item["t"]
            )
        except Exception:
            output["effective_dimension"] = None
    return output


def estimate_one_side(
    name: str,
    n: int,
    q: int,
    m: int,
    sigma: float,
) -> dict[str, Any]:
    (
        sage_log,
        LWE,
        Simulator,
        DiscreteGaussian,
        LWEParameters,
        reduction_algs,
        svp_models,
    ) = load_estimator_components()
    distribution = DiscreteGaussian(sigma, 0)
    parameters = LWEParameters(
        n=n,
        q=q,
        Xs=distribution,
        Xe=distribution,
        m=m,
        tag=name,
    )
    models: dict[str, Any] = {}
    for display_name, slug, svp_name in MODEL_SPECS:
        started = time.time()
        print(f"    {display_name} ...", flush=True)
        reduction_model = reduction_algs[r"\core"](
            svp_models[svp_name], svp_name
        )
        raw = LWE.estimate(
            parameters,
            red_cost_model=reduction_model,
            red_shape_model=Simulator.GSA,
            deny_list=DENY_LIST,
            catch_exceptions=True,
            quiet=True,
        )
        attacks = {
            key: _parse_attack(key, raw.get(key), n, sage_log)
            for key, _ in ATTACK_SPECS
        }
        minimum = min(item["log2_cost"] for item in attacks.values())
        models[slug] = {
            "display_name": display_name,
            "svp_name": svp_name,
            "attacks": attacks,
            "raw_minimum": minimum,
            "elapsed_seconds": time.time() - started,
        }
        print(
            f"      minimum={minimum:.3f} bits "
            f"({models[slug]['elapsed_seconds']:.1f} s)"
        )
    return {"n": n, "q": q, "m": m, "sigma": sigma, "models": models}
