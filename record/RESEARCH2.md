# RESEARCH PHASE 2 — the six unclaimed data items

Opened 2026-08-16 on the validated stack (validate_stack.py 7/7, bugs #12/#13 fixed).
Rules in force: /bridge execution rules 1–9. Every item pre-registers its falsifier here
before its first run. Compute tiers: LOCAL (this laptop) / KAGGLE (GPU notebook).

## Priority order (by decision value per compute hour)

### P1 — KAGGLE RUN 5: "the clean floors" (one notebook, ~4–6 h)
Two mandatory repairs that are also the two most valuable measurements:
  (a) **Item 3 core** — guarded-solver (deflated) remeasurement of the attraction exponent
      at beta = 2.4 and 2.2, L = 6..12. The installed fix removes the ~25% small-lambda
      contamination; this turns the "conservative floor" -2.86(20) into a clean number and
      decides how big the factor-6 excess over SU(3) really is.
      F: (a1) clean exponent still below free by >3 sigma at both betas -> excess is real;
         quote epsilon(SU(2), beta) vs Sternbeck's 0.16 as the open discrepancy.
         (a2) clean exponent collapses toward free -> the excess was solver contamination;
         retract the "sharpest in-window" framing. Either way: decided.
  (b) **Item 5 robustness** — r(lambda, D) under best-of-5-copy gauge fixing, n >= 96.
      The referee objection from 1211.3057 (r inherited from copy landing).
      F: (b1) r_bc consistent with r_fc -> the config-level law survives; claimable with
         the mandatory citation. (b2) |r_bc| < 0.2 -> the law was copy-landing; retract.
  Also in the same notebook: guarded Part-A s-means at 3 matched-ell points (repairs the
  A2 quantitative demotion).

### P2 — LOCAL (runs now, CPU-light)
  (c) **Item 2** — dim of lambda_min(Delta) vs volume, L = 3..6, 6 betas, n = 6.
      Dense eigvalsh (no eigsh bug). Measures BOTH ev[0] and ev[3] — the historical +0.25
      used ev[3] and the object choice was never justified; the ambiguity may BE the anomaly.
      F: (c1) dim volume-stable and object-consistent -> real anomalous exponent of the
         orbit softness; open a /bridge on deriving it. (c2) drifts with L or flips between
         ev[0]/ev[3] -> finite-size/object artifact; close with one line.
  (d) **Item 4** — the 0.82 ratio, ensembles: tr_H(Ric)/tr(T) at L = 3,4, n = 8, 4 betas.
      F: (d1) converges to a beta/L-stable value -> derivation target (rank structure of
         rho contraction); open a /bridge. (d2) drifts -> it was small-lattice numerology;
         close.
  (e) **Item 6** — residual-dimension asymmetry: K*V and Ric_min dims re-fit with the
      exact power-law estimator (kills the +8.3% log-linear bias) on n = 6/point.
      F: (e1) signs still opposite at >3 sigma each -> genuine small dimensions; report.
         (e2) consistent with zero/each other -> residuals were fit bias; close.

### P3 — KAGGLE RUN 6: the divisor, finally with power (~3 h)
  (f) chi(3,3) and chi(4,4) at L = 12, 16, n >= 24: the empirical d log a / d beta with
      loop-size convergence — feeds items 3 and 6, and closes the last scale-setting caveat.
      F: slopes converge in loop size -> one divisor number with error bars, final.
         still not converged at chi(4,4)/L=16 -> quote divisor as bound, permanent caveat.

### P4 — LOCAL ANALYTIC (after P2/P3 numbers firm up)
  (g) **Item 1** — derive the 0.0801 plateau. Bridge target: prove volume-independence of
      the interacting K*V trace contraction and identify the constant. First cheap test:
      compare 0.0801 against the free-field formula evaluated with a SELF-CONSISTENT
      massive propagator (one-parameter test, no new simulation).
  (h) **Method item** — the full blindness sweep: per-config Ric vs plaquette, Polyakov
      loop, action density at L = 4, n = 24. Either total blindness (strengthens the
      derivation) or the one thing the metric sees.

### Parked (explicitly)
  - Item 6's beta_c = 0.28 statistics (Dobrushin) — after P1–P4.
  - Any arrow-3 attempt — out of scope by standing rule.

## Status board
| item | tier | state |
|---|---|---|
| (a) clean exponent | KAGGLE 5 | DONE: 2.4 clean 4.8 sigma; 2.2 unresolved; eps grows with beta |
| (b) r under best-copy | KAGGLE 5 | DONE: I1 — law survives, replicated, item 5 CLOSED |
| (c) dim lambda_min(Delta) | LOCAL | RUNNING |
| (d) 0.82 ensembles | LOCAL | queued after (c) |
| (e) residual asymmetry | LOCAL | queued after (d) |
| (f) divisor with power | KAGGLE 6 | after run 5 |
| (g) 0.0801 derivation | ANALYTIC | ANSWERED: Haar-disorder constant, closure checks passed |
| (h) blindness sweep | LOCAL | after (e) |
