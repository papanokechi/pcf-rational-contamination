#!/usr/bin/env python3
"""
Phase 1 — ln(2) universality pre-filter and PSLQ scan.

For degrees d=2..6, sample 200 random PCF families (same generator as T1),
pre-filter for ln(2) proximity, then run PSLQ with ln(2)-enriched basis.

Output: results/ln2_universality.jsonl
"""

import argparse
import json
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

from mpmath import mp, mpf, nstr, log, pslq, pi, zeta, catalan, fabs

RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)
OUTPUT_FILE = RESULTS_DIR / "ln2_universality.jsonl"
RATIONAL_FILE = RESULTS_DIR / "rational_limits.jsonl"

# ── PCF evaluation (same as siarc_t1t6_relay.py) ──

def eval_pcf_forward(a_func, b_func, nterms, dps):
    """Standard Wallis: A_{-1}=1, A_0=b(0), B_{-1}=0, B_0=1, returns A_n/B_n."""
    with mp.workdps(dps + 60):
        b0 = b_func(0)
        if nterms == 0:
            return +mpf(b0)
        Ap = mpf(1)          # A_{-1}
        Ac = mpf(b0)         # A_0 = b(0)
        Bp = mpf(0)          # B_{-1}
        Bc = mpf(1)          # B_0 = 1
        for n in range(1, nterms + 1):
            an = a_func(n)
            bn = b_func(n)
            An = bn * Ac + an * Ap
            Bn = bn * Bc + an * Bp
            Ap, Ac = Ac, An
            Bp, Bc = Bc, Bn
            if n % 100 == 0 and abs(Bc) > mpf('1e300'):
                scale = mpf(1) / abs(Bc)
                Ap *= scale
                Ac *= scale
                Bp *= scale
                Bc *= scale
        if abs(Bc) < mpf('1e-100'):
            return None
    mp.dps = dps
    return +(Ac / Bc) if abs(Bc) > 0 else None


def make_poly_func(coeffs):
    def f(n):
        return sum(mpf(c) * mpf(n) ** i for i, c in enumerate(coeffs))
    return f


def generate_family(rng, d):
    a_coeffs = [rng.randint(-4, 4) for _ in range(d + 1)]
    b_coeffs = [rng.randint(-4, 4) for _ in range(d + 1)]
    if a_coeffs[-1] == 0:
        a_coeffs[-1] = 1
    if b_coeffs[-1] == 0:
        b_coeffs[-1] = 1
    if b_coeffs[0] == 0:
        b_coeffs[0] = 1
    return a_coeffs, b_coeffs


def safe_pslq(basis, maxcoeff=1000000, maxsteps=5000):
    try:
        return pslq(basis, maxcoeff=maxcoeff, maxsteps=maxsteps)
    except Exception:
        return None


# ── Known iteration-6 ln(2) relation reference ──
# All three prior hits: -13 + K + 13*ln2 - K*ln2 = 0  =>  (K-13)(1-ln2) = 0
REF_RELATION_COEFF_SIGNS = (-1, +1, +1, -1)  # signs of (1, K, ln2, K*ln2) coeffs


def run_prefilter(emit_rational=False):
    t0 = time.time()
    results = []
    rational_records = []

    for d in [2, 3, 4, 5, 6]:
        dt0 = time.time()
        seed_val = 700 + d * 100  # base_seed + 700 + d*100 per spec
        rng = random.Random(seed_val)
        n_sampled = 0
        n_prefilter = 0
        n_pslq_hits = 0
        n_near_miss = 0
        n_rational = 0

        print(f"[d={d}] scanning 200 families...", file=sys.stderr, flush=True)

        for fam_idx in range(200):
            a_coeffs, b_coeffs = generate_family(rng, d)
            a_func = make_poly_func(list(a_coeffs))
            b_func = make_poly_func(list(b_coeffs))
            n_sampled += 1

            # ── Convergence check at dps=60 ──
            mp.dps = 60
            try:
                v1 = eval_pcf_forward(a_func, b_func, 1999, 55)
                v2 = eval_pcf_forward(a_func, b_func, 2000, 55)
                if v1 is None or v2 is None:
                    results.append({
                        "degree": d, "family_id": fam_idx,
                        "prefilter_pass": False, "pslq_residual": None,
                        "relation": None, "coefficient_parity_match": False,
                        "discriminant": None
                    })
                    continue
                if abs(v1 - v2) > mpf('1e-50'):
                    results.append({
                        "degree": d, "family_id": fam_idx,
                        "prefilter_pass": False, "pslq_residual": None,
                        "relation": None, "coefficient_parity_match": False,
                        "discriminant": None
                    })
                    continue
            except Exception:
                results.append({
                    "degree": d, "family_id": fam_idx,
                    "prefilter_pass": False, "pslq_residual": None,
                    "relation": None, "coefficient_parity_match": False,
                    "discriminant": None
                })
                continue

            # ── Step 0 (Iteration 9 patch): rationality pre-screen ──
            mp.dps = 55
            K_low = eval_pcf_forward(a_func, b_func, 2000, 50)
            if K_low is None:
                results.append({
                    "degree": d, "family_id": fam_idx,
                    "prefilter_pass": False, "pslq_residual": None,
                    "relation": None, "coefficient_parity_match": False,
                    "discriminant": None, "rational_K": None
                })
                continue

            mp.dps = 50
            rat_rel = safe_pslq([mpf(1), K_low], maxcoeff=10000, maxsteps=500)
            if rat_rel is not None and rat_rel[1] != 0:
                rat_resid = float(fabs(rat_rel[0] + rat_rel[1] * K_low))
                if rat_resid < 1e-20:
                    # K is rational — classify as Type R, skip extended basis
                    K_rat = Fraction(-rat_rel[0], rat_rel[1])
                    n_rational += 1
                    results.append({
                        "degree": d, "family_id": fam_idx,
                        "prefilter_pass": False, "pslq_residual": None,
                        "relation": None, "coefficient_parity_match": False,
                        "discriminant": None,
                        "rational_K": str(K_rat),
                        "rational_residual": rat_resid
                    })
                    if emit_rational:
                        rational_records.append({
                            "degree": d, "family_id": fam_idx,
                            "K_rational": str(K_rat),
                            "residual_1K": rat_resid,
                            "dps": 50
                        })
                    continue

            ln2_val = log(2)
            prefilter_pass = False
            for p in range(1, 9):
                for q in range(1, 9):
                    if abs(K_low - mpf(p) / mpf(q) * ln2_val) < mpf('1e-8'):
                        prefilter_pass = True
                        break
                    # Also check K = rational + rational*ln(2)
                    # i.e. |K - p/q| < 1e-8 * |ln2| (near rational that might relate to ln2)
                if prefilter_pass:
                    break

            # Also check: does PSLQ with small basis find anything at dps=50?
            # Quick screen: [1, K, ln2]
            if not prefilter_pass:
                mp.dps = 55
                quick_rel = safe_pslq([mpf(1), K_low, ln2_val], maxcoeff=100, maxsteps=500)
                if quick_rel is not None and quick_rel[1] != 0:
                    resid = abs(quick_rel[0] + quick_rel[1] * K_low + quick_rel[2] * ln2_val)
                    if float(resid) < 1e-8:
                        prefilter_pass = True

            if not prefilter_pass:
                results.append({
                    "degree": d, "family_id": fam_idx,
                    "prefilter_pass": False, "pslq_residual": None,
                    "relation": None, "coefficient_parity_match": False,
                    "discriminant": None
                })
                continue

            n_prefilter += 1

            # ── Full PSLQ at dps=150 with ln(2)-enriched basis ──
            mp.dps = 160
            K_full = eval_pcf_forward(a_func, b_func, 2000, 150)
            if K_full is None:
                results.append({
                    "degree": d, "family_id": fam_idx,
                    "prefilter_pass": True, "pslq_residual": None,
                    "relation": None, "coefficient_parity_match": False,
                    "discriminant": None
                })
                continue

            mp.dps = 150
            ln2_150 = log(2)
            basis = [mpf(1), K_full, ln2_150, K_full * ln2_150, ln2_150 ** 2]
            rel = safe_pslq(basis, maxcoeff=1000000, maxsteps=5000)

            pslq_residual = None
            relation = None
            coeff_parity_match = False
            discriminant = None

            if rel is not None and rel[1] != 0:
                resid = abs(sum(mpf(c) * b for c, b in zip(rel, basis)))
                pslq_residual = float(resid)
                relation = list(rel)

                # Discriminant: norm of relation vector
                discriminant = float(sum(c * c for c in rel) ** 0.5)

                if pslq_residual < 1e-40:
                    n_pslq_hits += 1
                    # Check coefficient parity match with reference
                    signs = tuple(
                        1 if c > 0 else (-1 if c < 0 else 0)
                        for c in [rel[0], rel[1], rel[2], rel[3]]
                    )
                    coeff_parity_match = (signs == REF_RELATION_COEFF_SIGNS)
                elif pslq_residual < 1e-20:
                    n_near_miss += 1

            results.append({
                "degree": d, "family_id": fam_idx,
                "prefilter_pass": True,
                "pslq_residual": pslq_residual,
                "relation": relation,
                "coefficient_parity_match": coeff_parity_match,
                "discriminant": discriminant
            })

        dt = time.time() - dt0
        print(
            f"[d={d}] done in {dt:.1f}s: "
            f"{n_sampled} sampled, {n_rational} rational, {n_prefilter} prefilter pass, "
            f"{n_pslq_hits} PSLQ hits, {n_near_miss} near-misses",
            file=sys.stderr, flush=True
        )

        # Checkpoint: if total elapsed > 2400s, warn
        if time.time() - t0 > 2400:
            print(f"WARNING: elapsed {time.time()-t0:.0f}s > 2400s checkpoint",
                  file=sys.stderr, flush=True)

    # Write output
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        for r in results:
            f.write(json.dumps(r, default=str) + "\n")

    total_time = time.time() - t0
    total_hits = sum(1 for r in results if r["pslq_residual"] is not None and r["pslq_residual"] < 1e-40)
    total_rational = sum(1 for r in results if r.get("rational_K") is not None)
    print(f"\nPhase 1 complete in {total_time:.1f}s. "
          f"Total PSLQ hits: {total_hits}. "
          f"Total rational: {total_rational}. "
          f"Output: {OUTPUT_FILE}", file=sys.stderr, flush=True)

    if emit_rational and rational_records:
        with open(RATIONAL_FILE, "w", encoding="utf-8") as f:
            for r in rational_records:
                f.write(json.dumps(r, default=str) + "\n")
        print(f"Rational limits: {len(rational_records)} written to {RATIONAL_FILE}",
              file=sys.stderr, flush=True)

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--emit-rational", action="store_true",
                        help="Emit rational K families to rational_limits.jsonl")
    args = parser.parse_args()
    run_prefilter(emit_rational=args.emit_rational)
