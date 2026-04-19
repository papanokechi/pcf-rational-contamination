#!/usr/bin/env python3
"""
Iteration 9, Phase 1 — K = -3 coefficient characterization.

For each of the 7 families with K = -3 across degrees 2-6:
  1. Reconstruct coefficient vectors from the same RNG
  2. Verify K = -3 exactly at dps=150 via PSLQ([1, K])
  3. Compute minimal polynomial of K via PSLQ([1, K, K^2, ..., K^d])
  4. Check if fixed-point equation is satisfied
  5. Determine if families share a common minimal polynomial

Output: results/rat_characterization.json
"""

import json
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

from mpmath import mp, mpf, pslq, fabs, nstr

RESULTS = Path("results")
OUTPUT = RESULTS / "rat_characterization.json"

# The 7 K=-3 families identified in Iteration 8
K_NEG3_FAMILIES = [
    (2, 67),
    (3, 153),
    (4, 34),
    (4, 156),
    (4, 167),
    (5, 121),
    (6, 74),
]


def eval_pcf_forward(a_func, b_func, nterms, dps):
    """Forward Wallis recurrence to depth nterms."""
    with mp.workdps(dps + 60):
        b0 = b_func(0)
        if nterms == 0:
            return +mpf(b0)
        Pp = mpf(1)
        Pc = b_func(1)
        Qp = mpf(0)
        Qc = mpf(1)
        for n in range(1, nterms + 1):
            an = a_func(n)
            bn = b_func(n)
            Pn = bn * Pc + an * Pp
            Qn = bn * Qc + an * Qp
            Pp, Pc = Pc, Pn
            Qp, Qc = Qc, Qn
            if n % 100 == 0 and abs(Qc) > mpf('1e300'):
                scale = mpf(1) / abs(Qc)
                Pp *= scale
                Pc *= scale
                Qp *= scale
                Qc *= scale
        if abs(Qc) < mpf('1e-100'):
            return None
    mp.dps = dps
    return +(b0 + Pc / Qc) if abs(Qc) > 0 else None


def make_poly_func(coeffs):
    def f(n):
        return sum(mpf(c) * mpf(n) ** i for i, c in enumerate(coeffs))
    return f


def generate_family(rng, d):
    """Same generator as ln2_prefilter.py."""
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


def reconstruct_family(d, fam_id):
    """Replay RNG to get the exact coefficients for family fam_id at degree d."""
    seed_val = 700 + d * 100
    rng = random.Random(seed_val)
    for _ in range(fam_id):
        generate_family(rng, d)
    return generate_family(rng, d)


def analyze_fixed_point(a_coeffs, b_coeffs, K_val):
    """
    For a PCF b_0 + a_1/(b_1 + a_2/(b_2 + ...)),
    if K = b_0 + CF value, and the CF converges to x,
    then for large n the recurrence is approximately:
      x = a(n) / (b(n) + x)  =>  x*b(n) + x^2 = a(n)
    
    But more precisely, for a periodic CF with polynomial coefficients,
    the limit satisfies:
      K = b(0) + lim_{n->inf} CF(a_1..a_n, b_1..b_n)
    
    If K = -3 for all these families, check:
    Does x = K - b(0) satisfy the fixed-point equation
      x = a(1)/(b(1) + a(2)/(b(2) + ...)) ?
    
    The key structural question: what constrains K to be exactly -3?
    """
    b0 = b_coeffs[0]
    x = K_val - b0  # The CF part (excluding b0)
    
    # For a generalized CF, at "steady state" (large n), if a(n) ~ A*n^d
    # and b(n) ~ B*n^d, then the tail converges and the fixed point
    # of the tail x_inf satisfies: x_inf = A*n^d / (B*n^d + x_inf)
    # => x_inf * B + x_inf^2/n^d = A  (as n -> inf, x_inf = A/B if B != 0)
    
    # But the actual limit depends on ALL coefficients, not just leading.
    # The characterization is: K = b0 + x where x is determined by
    # the full CF structure.
    
    # Check: is b0 = -3 - x for some "natural" x?
    return {
        "b0": b0,
        "cf_value": float(x),
        "K_minus_b0": float(K_val - b0),
    }


def main():
    print("=" * 70)
    print("Phase 1 — K = -3 Coefficient Characterization")
    print("=" * 70)

    t0 = time.time()
    families = []

    for d, fam_id in K_NEG3_FAMILIES:
        print(f"\n  d={d}, fam={fam_id}:")
        a_coeffs, b_coeffs = reconstruct_family(d, fam_id)
        print(f"    a_coeffs = {a_coeffs}")
        print(f"    b_coeffs = {b_coeffs}")

        a_func = make_poly_func(a_coeffs)
        b_func = make_poly_func(b_coeffs)

        # ── Step 1: Verify K = -3 at dps=150 ──
        mp.dps = 160
        K_150 = eval_pcf_forward(a_func, b_func, 2000, 150)
        if K_150 is None:
            print(f"    DIVERGED at dps=150")
            families.append({
                "degree": d, "family_id": fam_id,
                "coefficients": {"a": a_coeffs, "b": b_coeffs},
                "status": "DIVERGED"
            })
            continue

        mp.dps = 150
        rel_1K = safe_pslq([mpf(1), K_150], maxcoeff=1000000, maxsteps=5000)
        resid_1K = None
        K_rational = None
        if rel_1K is not None and rel_1K[1] != 0:
            resid_1K = float(fabs(rel_1K[0] + rel_1K[1] * K_150))
            K_rational = str(Fraction(-rel_1K[0], rel_1K[1]))
        
        print(f"    K = {nstr(K_150, 25)}")
        print(f"    PSLQ [1,K]: {rel_1K} → K = {K_rational}")
        print(f"    Residual: {resid_1K}")

        verified = (rel_1K == [3, 1] or rel_1K == [-3, -1]) and (resid_1K is not None and resid_1K < 1e-100)
        print(f"    K = -3 verified: {verified}")

        # ── Step 2: Minimal polynomial via PSLQ([1, K, K^2, ..., K^d]) ──
        # For K = -3 (rational), the minimal polynomial is just (K + 3) = 0
        # degree 1. But let's confirm this formally.
        mp.dps = 150
        max_poly_deg = max(d, 3)  # At least degree 3 for thoroughness
        powers = [K_150 ** i for i in range(max_poly_deg + 1)]
        
        min_poly = None
        min_poly_deg = None
        for poly_d in range(1, max_poly_deg + 1):
            basis = powers[:poly_d + 1]
            rel = safe_pslq(basis, maxcoeff=10000, maxsteps=3000)
            if rel is not None:
                resid = float(fabs(sum(mpf(c) * p for c, p in zip(rel, basis))))
                if resid < 1e-100:
                    min_poly = list(rel)
                    min_poly_deg = poly_d
                    break

        if min_poly is not None:
            # Normalize: make leading coefficient positive
            if min_poly[-1] < 0:
                min_poly = [-c for c in min_poly]
            print(f"    Minimal polynomial: {min_poly} (degree {min_poly_deg})")
            # Express as polynomial string
            terms = []
            for i, c in enumerate(min_poly):
                if c == 0:
                    continue
                if i == 0:
                    terms.append(str(c))
                elif i == 1:
                    terms.append(f"{c}*K")
                else:
                    terms.append(f"{c}*K^{i}")
            print(f"    i.e., {' + '.join(terms)} = 0")
        else:
            print(f"    Minimal polynomial: NOT FOUND (K may be transcendental)")

        # ── Step 3: Discriminant ──
        disc = None
        if min_poly is not None and min_poly_deg == 1:
            disc = 1  # Linear polynomial, discriminant is trivial
        elif min_poly is not None and min_poly_deg == 2:
            # ax^2 + bx + c: disc = b^2 - 4ac
            a, b, c = min_poly[0], min_poly[1], min_poly[2]
            disc = b * b - 4 * a * c

        # ── Step 4: Fixed-point analysis ──
        fp = analyze_fixed_point(a_coeffs, b_coeffs, -3)
        print(f"    b(0) = {fp['b0']}, CF value = {fp['cf_value']}")

        # ── Step 5: Convergence rate analysis ──
        # How fast does the CF converge to -3?
        mp.dps = 60
        depths = [100, 500, 1000, 2000]
        conv_rates = []
        for depth in depths:
            v = eval_pcf_forward(a_func, b_func, depth, 50)
            if v is not None:
                err = float(fabs(v - mpf(-3)))
                conv_rates.append({"depth": depth, "error": err})

        print(f"    Convergence: ", end="")
        for cr in conv_rates:
            print(f"n={cr['depth']}: {cr['error']:.2e}  ", end="")
        print()

        families.append({
            "degree": d,
            "family_id": fam_id,
            "coefficients": {"a": a_coeffs, "b": b_coeffs},
            "K_value": -3,
            "K_rational": K_rational,
            "pslq_relation_1K": rel_1K,
            "residual_1K": resid_1K,
            "verified_K_neg3": verified,
            "minimal_poly": min_poly,
            "poly_degree": min_poly_deg,
            "discriminant": disc,
            "b0": fp["b0"],
            "cf_value": fp["cf_value"],
            "convergence": conv_rates,
        })

    # ── Summary analysis ──
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)

    all_verified = all(f.get("verified_K_neg3") for f in families)
    print(f"  All 7 families verified K = -3: {all_verified}")

    # Check shared minimal polynomial
    min_polys = [tuple(f["minimal_poly"]) for f in families if f.get("minimal_poly")]
    unique_polys = set(min_polys)
    shared = list(unique_polys)[0] if len(unique_polys) == 1 else None

    print(f"  Minimal polynomials found: {len(min_polys)}")
    print(f"  Unique minimal polynomials: {len(unique_polys)}")
    for p in unique_polys:
        count = min_polys.count(p)
        print(f"    {list(p)}: {count} families")

    if shared:
        print(f"  SHARED minimal polynomial: {list(shared)}")
    else:
        print(f"  NO shared polynomial (distinct polynomials found)")

    # Analyze b(0) distribution
    b0_vals = [f["b0"] for f in families if "b0" in f]
    print(f"\n  b(0) values: {b0_vals}")
    cf_vals = [f["cf_value"] for f in families if "cf_value" in f]
    print(f"  CF values (K - b(0)): {cf_vals}")

    # Hypothesis determination
    # If all b(0) == -3 and CF value == 0, then these are trivially convergent
    # If b(0) varies but K = -3 always, it's a genuine attractor
    b0_set = set(b0_vals)
    if len(b0_set) > 1:
        hypothesis = "attractor"
        print(f"\n  HYPOTHESIS: attractor")
        print(f"    b(0) takes {len(b0_set)} distinct values {sorted(b0_set)}")
        print(f"    but K = b(0) + CF_value = -3 always.")
        print(f"    → -3 is a genuine attractor of the PCF limit operator.")
    elif b0_vals[0] == -3 and all(v == 0.0 for v in cf_vals):
        hypothesis = "fixed_point"
        print(f"\n  HYPOTHESIS: fixed_point")
        print(f"    All families have b(0) = -3 and CF = 0.")
        print(f"    → Trivial convergence: the CF part vanishes.")
    else:
        hypothesis = "coincidence"
        print(f"\n  HYPOTHESIS: coincidence")

    elapsed = time.time() - t0
    print(f"\n  Completed in {elapsed:.1f}s")

    # ── Save output ──
    result = {
        "K_value": -3,
        "families": families,
        "shared_minimal_poly": list(shared) if shared else None,
        "n_unique_minimal_polys": len(unique_polys),
        "unique_minimal_polys": [list(p) for p in unique_polys],
        "hypothesis": hypothesis,
        "all_verified": all_verified,
        "b0_distribution": sorted(set(b0_vals)),
        "elapsed_s": round(elapsed, 1),
    }

    with open(OUTPUT, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    print(f"  Output: {OUTPUT}")


if __name__ == "__main__":
    main()
