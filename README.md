# pcf-rational-contamination

Supplementary code and data for:

**"Trivial Rational Contamination in PSLQ-Based Polynomial Continued Fraction Searches: Diagnosis, Correction, and a Pre-Screening Protocol"**

Papanokechi, submitted to *Experimental Mathematics*, 2026.

## Overview

We identify two mechanisms that produce rational limits in random polynomial continued fractions (PCFs)—vanishing first partial numerator ($a(1)=0$) and finite termination ($a(k)=0$ for $k \ge 2$)—and show that together they account for all 43 rational limits in a 1000-family survey across degrees 2–6. We propose a pre-screening protocol that eliminates both classes before extended-basis PSLQ, and document the five-iteration discovery arc—including a Wallis initialization bug that temporarily masked the triviality of the first mechanism—as a case study in machine-assisted self-correction.

## Repository Structure

```
scripts/
  ln2_prefilter.py         — PCF family generator + rationality pre-screen
  iter11_bugfix_reeval.py   — Corrected Wallis re-evaluation of all 43 families
  rat_characterize.py       — K=−3 family coefficient characterization
  iter10_algebraic.py       — Algebraic classification (Phase A–E)
  _verify_state.py          — AEAL governance state checker
  smoke_test_patch.py       — Smoke test for the Wallis patch (golden ratio, a(1)=0)

results/
  ln2_universality.jsonl    — Raw Phase 1 PSLQ output (1000 families)
  claims.jsonl              — Full AEAL claims log (135 total, 133 active)
  governance.jsonl          — Retraction and verification records
  iter11_bugfix_reeval.json — Corrected classification of all 43 families

pcf_rational_contamination_2026.tex — LaTeX source
pcf_rational_contamination_2026.pdf — Compiled manuscript
```

## Reproducing the Main Results

### Pre-screening protocol (Section 4)

```bash
python scripts/ln2_prefilter.py --emit-rational --all-degrees
```

### Full re-evaluation with corrected Wallis (Table 1)

```bash
python scripts/iter11_bugfix_reeval.py
```

### Verify AEAL governance state

```bash
python scripts/_verify_state.py --all
```

### Requirements

Python 3.10+, mpmath ≥ 1.3.0, sympy ≥ 1.12. Install via:

```bash
pip install -r requirements.txt
```

## The Wallis Initialization Bug (Section 2, Remark 3)

Our original evaluator computed the tail CF with initialization `Pc = b(1), Qc = 1` and returned `b(0) + Pc/Qc`, which adds an offset of exactly `b(1)` to every PCF value. The correct standard Wallis recurrence initializes `A₀ = b(0), B₀ = 1` and returns `Aₙ/Bₙ`. The bug was invisible to convergence diagnostics and rationality checks, and was detected only when the apparent theorem `K = b(0) + b(1)` was audited against the golden ratio φ = 1 + 1/(1 + 1/(1 + ⋯)). See Remark 3 in the paper for the full account, including a concrete numerical example.

## AEAL Governance

All claims, retractions, and verifications are recorded in machine-readable JSONL files under `results/`. The file `claims.jsonl` contains 135 claims across 11 relay iterations; 133 remain active and 2 were retracted (the Iteration 10 "theorem" claims, superseded by the bug discovery in Iteration 11). Each retraction entry includes a timestamp, reason code, and pointer to the superseding claim. The governance framework is described in:

- [Paper 7 — AI Behavior Science](https://doi.org/10.5281/zenodo.19562751)
- [Paper 8 — AEAL Framework](https://doi.org/10.5281/zenodo.19565086)

## Citation

```bibtex
@article{papanokechi2026rational,
  author  = {Papanokechi},
  title   = {Trivial Rational Contamination in {PSLQ}-Based Polynomial
             Continued Fraction Searches: Diagnosis, Correction,
             and a Pre-Screening Protocol},
  journal = {Experimental Mathematics},
  year    = {2026},
  note    = {Submitted}
}
```

## License

MIT
