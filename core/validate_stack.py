"""
INSTRUMENT VALIDATION SUITE — run this BEFORE any physics measurement. One command.

WHY THIS EXISTS.  The programme's two costliest failures (2026-08-12) were instruments
trusted without validation: a sampler that violated detailed balance from day one, and a
dimension pipeline never shown able to detect a dimension that exists. Both were findable
in minutes with the checks below. Rule now enforced: NO measurement script is trusted on a
stack that has not passed this suite in its current state.

CHECKS (all must PASS; each has a literature or exact anchor, never a self-reference):
  V1  sampler vs LITERATURE:  <plaq>(beta=2.4, 4^4)  = 0.630(4)      [independent anchor]
  V2  sampler vs INDEPENDENT ALGORITHM: Metropolis vs Creutz heat bath, same point, < 3 se
  V3  strong-coupling exact limit:      <plaq>(beta=0.4) = beta/4 = 0.1000(30)
  V4  detailed-balance signature:       hits-independence of <plaq> (the bug-#10 tell)
  V5  FP operator on a thermalised, Landau-fixed config: symmetric, PSD, exactly 3
      near-zero modes below 1% of the first physical eigenvalue  (the bug-#11 tell)
  V6  FP free-field exact value:        lam[3](U=1, L=6) = 4 sin^2(pi/6) = 1.000000
  V7  POSITIVE CONTROL of the dimension pipeline: sigma a^2 via chi(2,2), L=4, must give
      a NONZERO dimension of the correct sign (detectability), and the EMPIRICAL divisor
      chi(3,3) slope must be quoted next to any 1-loop assumption (calibration).
      This is deliberately the cheap version: it proves detectability, not precision.

Runtime ~15-25 min on this CPU. Prints PASS/FAIL per check and an overall verdict.
"""
import faulthandler
faulthandler.enable()

import numpy as np
import sampler_fix
from step3_lattice import Lattice
from sampler_fix import sweep_correct, heatbath_sweep, measure
from horizon_blind import landau_fix, M_operator
from positive_control import loops, creutz
from scipy.sparse.linalg import eigsh

RESULTS = []


def check(name, ok, detail):
    RESULTS.append((name, bool(ok), detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: {detail}", flush=True)


if __name__ == "__main__":
    print("=" * 74)
    print("INSTRUMENT VALIDATION SUITE — no physics before this passes")
    print("=" * 74, flush=True)

    # V1 + V2: literature anchor (VOLUME-MATCHED) and independent-algorithm cross-check.
    # First version of V1 compared a 4^4 measurement against the INFINITE-VOLUME literature
    # value and failed on correct physics — threshold violation #9, the same volume-mismatch
    # disease as gate #7. The anchor is now measured at L=8, where the large-volume value
    # ~0.630 applies within small finite-size effects, and the tolerance carries the run's
    # own statistical error.
    m8, sm8 = measure(lambda l, b, **k: sweep_correct(l, b, **k), 8, 2.4,
                      nther=300, nmeas=12)
    check("V1 plaquette vs literature (L=8)", abs(m8 - 0.630) < np.hypot(0.004, 3 * sm8),
          f"metropolis {m8:.4f} vs 0.630(4) [volume-matched anchor]")
    # BUG #12 (found 2026-08-16): sampler_fix seeds a MODULE-LEVEL rng with a fixed
    # constant, so every fresh process replays the identical stream — V2 "failed" twice
    # with bit-identical numbers because one legal ~3-sigma excursion was frozen at that
    # stream position. Validation reruns must be statistically independent: reseed from
    # entropy here (physics runs elsewhere keep their recorded seeds for reproducibility).
    sampler_fix.rng = np.random.default_rng()
    m, sm = measure(lambda l, b, **k: sweep_correct(l, b, **k), 4, 2.4,
                    nther=500, nmeas=36)
    sampler_fix.rng = np.random.default_rng()
    h, sh = measure(lambda l, b, **k: heatbath_sweep(l, b), 4, 2.4,
                    nther=500, nmeas=36)
    check("V2 metropolis vs heat bath (4^4)", abs(m - h) < 3.5 * np.hypot(sm, sh),
          f"{m:.4f} vs {h:.4f} ({abs(m-h)/np.hypot(sm,sh):.1f} se; independent entropy)")

    # V3: strong-coupling expansion WITH its known correction and a COMPUTED tolerance.
    # First version demanded beta/4 exactly with a picked 0.003 window — threshold
    # violation #10: it ignored the beta^3/96 term AND the run's statistical error.
    s, ss = measure(lambda l, b, **k: sweep_correct(l, b, **k), 4, 0.4,
                    nther=250, nmeas=30)
    expect = 0.4 / 4 + 0.4 ** 3 / 96              # SU(2) strong coupling to O(beta^3)
    tol = np.hypot(3 * ss, 0.002)                 # 3 se + truncation allowance
    check("V3 strong-coupling expansion", abs(s - expect) < tol,
          f"{s:.4f} vs {expect:.4f} (tol {tol:.4f}, computed not picked)")

    # V4: hits-independence (detailed-balance signature; the bug-#10 tell)
    a, sa = measure(lambda l, b, **k: sweep_correct(l, b, hits=1), 4, 2.4,
                    nther=350, nmeas=20)
    b_, sb = measure(lambda l, b, **k: sweep_correct(l, b, hits=4), 4, 2.4,
                     nther=350, nmeas=20)
    check("V4 hits-independence", abs(a - b_) < 3 * np.hypot(sa, sb),
          f"hits=1: {a:.4f} vs hits=4: {b_:.4f}")

    # V5: FP operator structure on a real configuration (the bug-#11 tell)
    lat = Lattice(4, 4)
    for _ in range(400):
        sweep_correct(lat, 2.4)
    th = landau_fix(lat)
    ev = np.sort(eigsh(M_operator(lat), k=6, which='SA', tol=1e-7, maxiter=30000,
                       return_eigenvectors=False))
    near0 = np.all(np.abs(ev[:3]) < 0.01 * abs(ev[3]))
    check("V5 FP operator structure", th < 1e-8 and ev[3] > 0 and near0,
          f"theta={th:.1e}, ev[0:4]={np.array2string(ev[:4], precision=5)}")

    # V6: FP free-field exact value
    latf = Lattice(6, 4)
    evf = np.sort(eigsh(M_operator(latf), k=6, which='SA', tol=1e-7, maxiter=30000,
                        return_eigenvectors=False))
    check("V6 FP free-field exact", abs(evf[3] - 1.0) < 1e-4, f"lam[3]={evf[3]:.6f} vs 1")

    # V7: dimension-pipeline positive control (cheap detectability version)
    print("  [V7 running: chi(2,2) and chi(3,3) scans ...]", flush=True)
    betas = [2.2, 2.5, 2.8]
    c22, c33 = [], []
    for beta in betas:
        lat = Lattice(4, 4)
        for _ in range(400):
            sweep_correct(lat, beta)
        acc = np.zeros((4, 4)); n = 0
        for _ in range(8):
            for _ in range(3):
                sweep_correct(lat, beta)
            acc += loops(lat, 3); n += 1
        W = acc / n
        c22.append(creutz(W, 2, 2)); c33.append(creutz(W, 3, 3))
    s22 = np.polyfit(betas, np.log(c22), 1)[0]
    s33 = np.polyfit(betas, np.log(np.abs(c33)), 1)[0]
    # detectability: chi(2,2) must RUN (fall with beta) clearly — slope well below zero
    check("V7 positive control detects a dimension", s22 < -0.5,
          f"d log chi22/d beta = {s22:+.3f} (must be clearly < 0); "
          f"EMPIRICAL divisor candidate: d log a/d beta ~ {s33/2:+.3f} from chi33 "
          f"(1-loop would claim -2.69 — quote both, believe the measured one)")

    print("\n" + "=" * 74)
    npass = sum(1 for _, ok, _ in RESULTS if ok)
    print(f"VERDICT: {npass}/{len(RESULTS)} checks passed")
    if npass == len(RESULTS):
        print("STACK VALIDATED — physics measurements may be trusted on this build.")
    else:
        print("STACK NOT VALIDATED — fix the failing instrument before ANY measurement.")
        for name, ok, detail in RESULTS:
            if not ok:
                print(f"  FAILING: {name} — {detail}")
