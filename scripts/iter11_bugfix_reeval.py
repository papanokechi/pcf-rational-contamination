#!/usr/bin/env python3
"""
Iteration 11 — Post-bug-fix complete re-evaluation of all 43 rational hits.

The bug: eval_pcf_forward initializes Pc = b(1) instead of Pc = b(0).
This adds exactly b(1) to the standard PCF value.
Corrected value = code_value - b(1) = standard Wallis value.

For each of the 43 families:
  1. Compute corrected K (= code_K - b(1), or equivalently standard Wallis)
  2. Classify:
     - a(1)=0: K_corrected = b(0) (trivially rational)
     - a(k)=0, k>=2: K_corrected = finite CF value (rational)
     - otherwise: check if still rational via PSLQ on [1, K_corrected]
  3. For a(1)=0 cases: verify K_corrected == b(0) exactly
  4. For finite CF cases: verify K_corrected is rational
"""

import json
import random
import time
from fractions import Fraction
from pathlib import Path

from mpmath import mp, mpf, pslq, fabs

RESULTS = Path("results")


def eval_poly(coeffs, n):
    return sum(c * n**i for i, c in enumerate(coeffs))


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


def reconstruct_family(d, fam_id):
    seed_val = 700 + d * 100
    rng = random.Random(seed_val)
    for _ in range(fam_id):
        generate_family(rng, d)
    return generate_family(rng, d)


def standard_wallis_exact(a_coeffs, b_coeffs, N):
    """Standard Wallis: A_{-1}=1, A_0=b(0), B_{-1}=0, B_0=1."""
    b0 = Fraction(eval_poly(b_coeffs, 0))
    Ap = Fraction(1)
    Ac = b0
    Bp = Fraction(0)
    Bc = Fraction(1)
    for n in range(1, N + 1):
        an = Fraction(eval_poly(a_coeffs, n))
        bn = Fraction(eval_poly(b_coeffs, n))
        An = bn * Ac + an * Ap
        Bn = bn * Bc + an * Bp
        Ap, Ac = Ac, An
        Bp, Bc = Bc, Bn
    return Ac / Bc if Bc != 0 else None


def direct_backward_eval(a_coeffs, b_coeffs, depth):
    """Backward recursion (ground truth)."""
    val = Fraction(eval_poly(b_coeffs, depth))
    for n in range(depth - 1, 0, -1):
        an1 = Fraction(eval_poly(a_coeffs, n + 1))
        bn = Fraction(eval_poly(b_coeffs, n))
        if val == 0:
            return None
        val = bn + an1 / val
    a1 = Fraction(eval_poly(a_coeffs, 1))
    b0 = Fraction(eval_poly(b_coeffs, 0))
    if val == 0:
        return None
    return b0 + a1 / val


def standard_wallis_mp(a_coeffs, b_coeffs, N, dps=150):
    """Standard Wallis with mpmath for high precision."""
    mp.dps = dps + 60
    b0 = mpf(eval_poly(b_coeffs, 0))
    Ap = mpf(1)
    Ac = b0
    Bp = mpf(0)
    Bc = mpf(1)
    for n in range(1, N + 1):
        an = mpf(eval_poly(a_coeffs, n))
        bn = mpf(eval_poly(b_coeffs, n))
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
    mp.dps = dps
    if abs(Bc) < mpf('1e-100'):
        return None
    return +(Ac / Bc)


def find_integer_roots(coeffs, lo=-10, hi=10):
    """Find integer roots of polynomial in [lo, hi]."""
    return [n for n in range(lo, hi + 1) if eval_poly(coeffs, n) == 0]


def main():
    t0 = time.time()
    print("=" * 70)
    print("ITERATION 11 — POST-BUG-FIX RE-EVALUATION")
    print("=" * 70)

    # Load all 43 hits from the original scan
    hits = [json.loads(l) for l in open(RESULTS / "ln2_universality.jsonl")
            if json.loads(l).get("relation") is not None]

    print(f"\n  Total hits to re-evaluate: {len(hits)}")

    categories = {
        "trivial_a1_zero": [],      # a(1)=0 → K = b(0) trivially
        "finite_cf": [],             # a(k)=0 for k>=2 → finite CF
        "still_rational": [],        # K ∈ Q but not via above mechanisms
        "no_longer_rational": [],    # Was "rational" only due to bug offset
        "divergent": [],             # Doesn't converge after fix
    }

    print(f"\n  {'d':>2} {'fam':>4} {'a(1)':>5} {'b(0)':>5} {'b(1)':>5} "
          f"{'old(code)':>12} {'corrected':>12} {'mechanism':>18} {'K_exact':>15}")
    print(f"  {'─'*2} {'─'*4} {'─'*5} {'─'*5} {'─'*5} "
          f"{'─'*12} {'─'*12} {'─'*18} {'─'*15}")

    for h in hits:
        d = h["degree"]
        fid = h["family_id"]
        a_coeffs, b_coeffs = reconstruct_family(d, fid)

        a1 = eval_poly(a_coeffs, 1)
        b0 = eval_poly(b_coeffs, 0)
        b1 = eval_poly(b_coeffs, 1)
        a_roots = find_integer_roots(a_coeffs, 0, 20)

        # Old code value
        old_K = float(h.get("K", h.get("value", 0)))

        # Corrected: standard Wallis at high precision
        mp.dps = 150
        K_corrected = standard_wallis_mp(a_coeffs, b_coeffs, 2000, dps=150)
        if K_corrected is None:
            categories["divergent"].append({
                "degree": d, "family_id": fid,
                "a_coeffs": a_coeffs, "b_coeffs": b_coeffs,
            })
            print(f"  {d:2d} {fid:4d} {a1:5d} {b0:5d} {b1:5d} "
                  f"{old_K:12.4f} {'DIVERGENT':>12} {'divergent':>18}")
            continue

        # Classify
        if a1 == 0:
            # Should give K = b(0) trivially
            K_exact = standard_wallis_exact(a_coeffs, b_coeffs, 20)
            if K_exact is not None and K_exact == Fraction(b0):
                mechanism = "trivial_a1_zero"
                k_str = str(K_exact)
            else:
                # Hmm, a(1)=0 but K != b(0)? Check deeper
                K_deep = standard_wallis_exact(a_coeffs, b_coeffs, 100)
                if K_deep == Fraction(b0):
                    mechanism = "trivial_a1_zero"
                    k_str = str(K_deep)
                else:
                    mechanism = "a1_zero_unexpected"
                    k_str = f"{float(K_corrected):.6f}"

            categories.get(mechanism, categories["still_rational"]).append({
                "degree": d, "family_id": fid,
                "a_coeffs": a_coeffs, "b_coeffs": b_coeffs,
                "K_corrected": float(K_corrected),
                "K_exact": str(K_exact) if K_exact else None,
                "b_0": b0,
            })
        elif any(r >= 2 for r in a_roots):
            # Finite CF
            k_root = min(r for r in a_roots if r >= 2)
            K_exact = standard_wallis_exact(a_coeffs, b_coeffs, k_root)
            mechanism = "finite_cf"
            k_str = str(K_exact) if K_exact else f"{float(K_corrected):.6f}"
            categories["finite_cf"].append({
                "degree": d, "family_id": fid,
                "a_coeffs": a_coeffs, "b_coeffs": b_coeffs,
                "K_corrected": float(K_corrected),
                "K_exact": str(K_exact) if K_exact else None,
                "termination_point": k_root,
                "b_0": b0, "b_1": b1,
            })
        else:
            # Neither a(1)=0 nor finite CF — check if still rational
            mp.dps = 150
            rel = pslq([mpf(1), K_corrected], maxcoeff=10000)
            if rel and rel[1] != 0:
                K_rat = Fraction(-rel[0], rel[1])
                # Verify
                residual = abs(K_corrected - float(K_rat))
                if residual < mpf('1e-100'):
                    mechanism = "still_rational"
                    k_str = str(K_rat)
                    categories["still_rational"].append({
                        "degree": d, "family_id": fid,
                        "a_coeffs": a_coeffs, "b_coeffs": b_coeffs,
                        "K_corrected": float(K_corrected),
                        "K_exact": str(K_rat),
                        "b_0": b0, "b_1": b1,
                        "a_roots": a_roots,
                    })
                else:
                    mechanism = "no_longer_rational"
                    k_str = f"{float(K_corrected):.6f}"
                    categories["no_longer_rational"].append({
                        "degree": d, "family_id": fid,
                        "K_corrected": float(K_corrected),
                    })
            else:
                mechanism = "no_longer_rational"
                k_str = f"{float(K_corrected):.6f}"
                categories["no_longer_rational"].append({
                    "degree": d, "family_id": fid,
                    "K_corrected": float(K_corrected),
                })

        print(f"  {d:2d} {fid:4d} {a1:5d} {b0:5d} {b1:5d} "
              f"{old_K:12.4f} {float(K_corrected):12.4f} "
              f"{mechanism:>18} {k_str:>15}")

    # ── Summary ──
    print(f"\n{'=' * 70}")
    print(f"CLASSIFICATION SUMMARY")
    print(f"{'=' * 70}")

    for cat, items in categories.items():
        print(f"  {cat}: {len(items)}")

    total = sum(len(v) for v in categories.values())
    print(f"  TOTAL: {total}")

    # Detail on non-trivial cases
    if categories["still_rational"]:
        print(f"\n  *** NON-TRIVIAL RATIONAL LIMITS (post-fix): ***")
        for f in categories["still_rational"]:
            print(f"    d={f['degree']}, fam={f['family_id']}: "
                  f"K = {f['K_exact']}, b(0) = {f['b_0']}, "
                  f"a_roots = {f.get('a_roots', [])}")

    if categories["no_longer_rational"]:
        print(f"\n  *** NO LONGER RATIONAL (bug-induced false positives): ***")
        for f in categories["no_longer_rational"]:
            print(f"    d={f['degree']}, fam={f['family_id']}: "
                  f"K_corrected = {f['K_corrected']:.10f}")

    # ── Verify a(1)=0 → K = b(0) ──
    if categories["trivial_a1_zero"]:
        all_trivial = all(
            f["K_exact"] == str(f["b_0"])
            for f in categories["trivial_a1_zero"]
            if f.get("K_exact")
        )
        print(f"\n  All a(1)=0 families give K = b(0): {all_trivial}")

    # ── Check finite CF values ──
    if categories["finite_cf"]:
        print(f"\n  Finite CF families (post-fix):")
        for f in categories["finite_cf"]:
            old_rational = float(Fraction(f["K_exact"])) if f["K_exact"] else None
            print(f"    d={f['degree']}, fam={f['family_id']}: "
                  f"K = {f['K_exact']} = {old_rational}, "
                  f"terminates at a({f['termination_point']})=0")

    elapsed = time.time() - t0

    # Save results
    output = {
        "bug_description": "eval_pcf_forward initializes Pc=b(1) instead of Pc=b(0), adding +b(1) offset",
        "correction": "Use standard Wallis: A_0=b(0), B_0=1",
        "categories": {k: len(v) for k, v in categories.items()},
        "details": {k: v for k, v in categories.items()},
        "elapsed_s": round(elapsed, 1),
    }
    with open(RESULTS / "iter11_bugfix_reeval.json", "w") as f:
        json.dump(output, f, indent=2, default=str)

    print(f"\n  Output: results/iter11_bugfix_reeval.json")
    print(f"  Completed in {elapsed:.1f}s")


if __name__ == "__main__":
    main()
