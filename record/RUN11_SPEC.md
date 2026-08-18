# Run 11 — pre-registered specification: is r = F(f_IR) a law?

*Registered 2026-08-18, before any data. Prompted by hostile external review.
Nothing below may change after the first ensemble is generated.*

## What run 10 established, and what it did not

Run 10 established: **the infrared coherence of the lowest FP eigenmode predicts
whether the horizon–gluon correlation is active** (two out-of-sample points, both
hit). It did **not** establish `r = F(f_IR)` as a universal single-variable law —
the shape of F between the OFF and ON plateaus is inferred, not measured, and
`f_IR` is constructed from the same eigenmode whose eigenvalue enters `r`.

## Design (three legs, all required)

**Leg 1 — dense blinded scan of the transition.**
10–15 ensembles chosen (by the measured `f_IR(β, L)` grid) to span
`f_IR ≈ 0.30 … 0.95` with at least 5 points inside `0.45–0.70`. For each:
measure per-config `(λ_min, D(p_min))`, compute partial `r` (plaquette removed).
The predicted `r` for every ensemble is written down **before** any `r` is
computed (holdout: `r` values unblinded only after all predictions are logged).

**Leg 2 — independence of the coherence variable.** Three coherence measures per
ensemble:
- `f₁` = f_IR of ψ_min (the original),
- `f₂` = f_IR of ψ₂, ψ₃ (higher modes — same operator, different vector),
- `f₃` = a non-spectral IR coherence observable built from the gauge field alone
  (candidate: normalized IR fraction of the momentum-space link fluctuation power;
  fixed exactly before data).
The law survives only if `f₃` gates `r` the same way `f₁` does.

**Leg 3 — conditional independence.** Test `r ⊥ (β, L) | f_IR`: at matched
`f_IR` (different (β, L) reaching the same coherence), `r` must agree within
errors. Any residual dependence on (β, L) at fixed `f_IR` falsifies
single-variable control.

## Additionally demanded by review (same run)

- Replace the correlation coefficient with a physical regression:
  fit `D(p_min) = a(f_IR)·g(λ_min) + ε` per ensemble; report `a(f_IR)` with
  units. The "locked −0.52" must reappear as systematic behaviour of `a`,
  or it is a distributional artifact.
- Continuous-m minimization of χ²(m) for the massive-screening alternative
  (closes the grid-search objection against the disorder-constant claim).

## Falsifiers (fixed now)

- F1: any ensemble with `f₁ < 0.45` showing `|r| > 0.25` at >2σ → threshold
  picture dead.
- F2: `f₃` fails to gate `r` (transition absent or displaced by >0.15 in f) →
  the "law" is a spectral self-correlation; claim reverts permanently to
  "crossover-associated correlation".
- F3: matched-`f_IR` pairs differing in `r` by >3σ → single-variable control dead.

## Verdict language (fixed now)

All three legs pass → "coherence-governed law" is earned and restored.
Any leg fails → the paper's current wording ("coherence-gated correlation,
threshold behaviour consistent with single-variable control") is the ceiling.
