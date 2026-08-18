"""
P2(c) — the +0.25 anomaly: dimension of lambda_min(Delta), volume scan, BOTH objects.

HISTORY. horizon_scale.py (2026-08-13) measured dim(lambda_min(Delta)) = +0.25 +- 0.02 at
L=4 using np.sort(eigvalsh(Delta))[3] — the [3] was inherited from the U=1 case (3 exact
zero modes) but for INTERACTING configs Delta has a trivial kernel and ev[0] is already
physical. The object choice was never justified; the ambiguity may BE the anomaly.

THIS SCAN. L = 3..6, beta = 2.0..3.0 (6 points), n = 6 configs per point, corrected
sampler, entropy seeding. Dense eigvalsh (no eigsh, no bug #13 exposure). For each config
record ev[0] AND ev[3]; fit both dimensions per L with the EXACT power-law estimator
(nonlinear fit of Q = A * exp(n * DLOGA_EMP * beta), which has none of the +8.3%
log-linear bias) using DLOGA_EMP = -1.17031.

PRE-REGISTERED
 (C1) dim volume-stable (spread across L within errors) AND ev[0]/ev[3] dims agree
      -> a real anomalous exponent of the gauge-orbit soft edge; open a /bridge on it.
 (C2) dim drifts with L, or ev[0] and ev[3] disagree
      -> finite-size/object artifact; the +0.25 closes with one line.
Prints numbers and the rule; no conclusion.
"""
import faulthandler
faulthandler.enable()

import numpy as np
import sampler_fix
from step3_lattice import Lattice, adjoint_matrix
from sampler_fix import sweep_correct

DLOGA_EMP = -1.17031


def delta_dense(lat):
    D, V = lat.D, lat.V
    R = np.stack([adjoint_matrix(lat.U[mu]) for mu in range(D)])
    NV = 3 * V
    G = np.zeros((D * NV, NV))
    for mu in range(D):
        f = lat.fwd[mu]
        for n in range(V):
            r0 = (mu * V + n) * 3
            G[r0:r0 + 3, 3 * f[n]:3 * f[n] + 3] += R[mu, n]
            G[r0:r0 + 3, 3 * n:3 * n + 3] -= np.eye(3)
    return G.T @ G


def fit_dim_exact(betas, vals, errs):
    """Nonlinear LS for Q = A exp(n * DLOGA * beta): exact for pure power laws."""
    b = np.array(betas, float); y = np.array(vals, float); e = np.array(errs, float)
    # init from log-linear
    sl = np.polyfit(b, np.log(y), 1)[0]
    n0 = sl / DLOGA_EMP
    A0 = np.exp(np.mean(np.log(y) - n0 * DLOGA_EMP * b))
    from scipy.optimize import curve_fit
    f = lambda bb, A, n: A * np.exp(n * DLOGA_EMP * bb)
    p, cov = curve_fit(f, b, y, p0=[A0, n0], sigma=e, absolute_sigma=True, maxfev=20000)
    return p[1], float(np.sqrt(cov[1, 1]))


if __name__ == "__main__":
    betas = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
    NCFG = 6
    print("P2(c): dim(lambda_min(Delta)) vs volume, ev[0] and ev[3], exact-fit estimator")
    out = {}
    for L in (3, 4, 5, 6):
        m0, e0, m3, e3 = [], [], [], []
        for beta in betas:
            v0, v3 = [], []
            for c in range(NCFG):
                sampler_fix.rng = np.random.default_rng()      # entropy per config
                lat = Lattice(L, 4)
                for _ in range(400):
                    sweep_correct(lat, beta)
                ev = np.sort(np.linalg.eigvalsh(delta_dense(lat)))
                v0.append(ev[0]); v3.append(ev[3])
            m0.append(np.mean(v0)); e0.append(np.std(v0, ddof=1) / np.sqrt(NCFG))
            m3.append(np.mean(v3)); e3.append(np.std(v3, ddof=1) / np.sqrt(NCFG))
            print(f"  L={L} beta={beta:.1f}  ev0={m0[-1]:.4f}+-{e0[-1]:.4f}"
                  f"  ev3={m3[-1]:.4f}+-{e3[-1]:.4f}", flush=True)
        d0, s0 = fit_dim_exact(betas, m0, e0)
        d3, s3 = fit_dim_exact(betas, m3, e3)
        out[L] = (d0, s0, d3, s3)
        print(f"  L={L}: dim(ev0) = {d0:+.4f}+-{s0:.4f}   dim(ev3) = {d3:+.4f}+-{s3:.4f}",
              flush=True)

    print("\n  SUMMARY (dim vs L):")
    for L, (d0, s0, d3, s3) in out.items():
        print(f"    L={L}:  ev0 {d0:+.4f}+-{s0:.4f}   ev3 {d3:+.4f}+-{s3:.4f}")
    print("\n  C1 volume-stable & objects agree -> real exponent, open a /bridge")
    print("  C2 drifts or objects disagree    -> artifact; close the +0.25 with one line")
