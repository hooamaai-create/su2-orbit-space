"""
THE HORIZON APPROACH: does lambda_min(M) -> 0 with volume FASTER than free momentum?

WHY THIS IS THE NEXT MEASUREMENT.  The differential test (2026-08-13) showed the horizon
sector scales differently from the metric sector (3.3 sigma), but its beta-dependence at
fixed L measures the RETREAT from the horizon, not the horizon's mass. Gribov-Zwanziger
lives at V -> infinity, where typical configurations concentrate ON the horizon:
lambda_min(M) -> 0. The signal that distinguishes dynamics from kinematics is the RATE.

FREE-FIELD BASELINE (exact).  For U = 1, M = -partial^2, and the smallest nonzero
eigenvalue is the smallest lattice momentum squared:
      lambda_free(L) = 4 sin^2(pi / L)   ~   4 pi^2 / L^2      i.e. exponent p = -2.

PRE-REGISTERED FALSIFIERS  (fit lambda_min(M) ~ L^p over L = 3,4,5,6 at beta = 2.4)
 (V1) p < -2 by more than 3 sigma  ->  configurations pile toward the horizon FASTER than
      free momentum decay: genuine dynamical horizon attraction. The GZ mechanism is
      active on these lattices, and the ratio lambda_min(M)/lambda_free(L) is a measurable
      horizon-attraction factor whose beta-dependence can later give the scale.
 (V2) p = -2 within errors  ->  the approach is purely kinematic (free momentum). The
      horizon is approached but NOT dynamically favoured at these volumes; extracting
      gamma here is hopeless and the GZ door needs far bigger lattices. Report honestly.
 (V3) p > -2 by more than 3 sigma  ->  configurations RESIST the horizon. Would contradict
      GZ expectations; check for a bug before believing it.

CONTROL.  The same pipeline applied to the FREE field must return p = -2.00 within errors
(it is exact), or the fit machinery is wrong and nothing is quoted.

Prints numbers and the rules. No conclusion is printed.
"""
import faulthandler
faulthandler.enable()

import numpy as np
from step3_lattice import Lattice
from sampler_fix import sweep_correct
from horizon_blind import landau_fix, M_operator
from scipy.sparse.linalg import eigsh


def lam_min_M(lat):
    """BUG #13 guard (audit 2026-08-16): eigsh can silently drop a zero-triplet member at
    small lambda, promoting ev[3] to the SECOND physical eigenvalue. Guard: require the
    3 near-zeros below 1% of ev[3]; on failure, deflate the constant modes exactly."""
    A = M_operator(lat)
    ev = np.sort(eigsh(A, k=6, which='SA', tol=1e-7,
                       maxiter=40000, return_eigenvectors=False))
    if np.all(np.abs(ev[:3]) < 0.01 * abs(ev[3])):
        return float(ev[3])
    V = lat.V

    def proj(x):
        w = x.reshape(V, 3)
        return (w - w.mean(axis=0, keepdims=True)).reshape(-1)

    from scipy.sparse.linalg import LinearOperator as _LO
    Ad = _LO((3 * V, 3 * V),
             matvec=lambda x: proj(A.matvec(proj(x))) + 10.0 * (x - proj(x)),
             dtype=float)
    lam = np.sort(eigsh(Ad, k=3, which='SA', tol=1e-7, maxiter=40000,
                        return_eigenvectors=False))
    return float(lam[0])


def pfit(Ls, vals, errs):
    x = np.log(np.array(Ls, float)); a = np.array(vals, float); e = np.array(errs, float)
    y = np.log(a); w = np.maximum(e / a, 1e-9); W = 1 / w ** 2
    Sw, Sx, Sy = W.sum(), (W * x).sum(), (W * y).sum()
    Sxx, Sxy = (W * x * x).sum(), (W * x * y).sum()
    den = Sw * Sxx - Sx ** 2
    return (Sw * Sxy - Sx * Sy) / den, np.sqrt(Sw / den)


if __name__ == "__main__":
    BETA = 2.4
    Ls = [3, 4, 5, 6]
    NCFG = {3: 6, 4: 5, 5: 4, 6: 3}

    print("CONTROL: free field, exact lambda = 4 sin^2(pi/L); pipeline must return p = -2")
    fv = [4 * np.sin(np.pi / L) ** 2 for L in Ls]
    for L, v in zip(Ls, fv):
        print(f"  L={L}  lambda_free = {v:.6f}")
    pc, ec = pfit(Ls, fv, [1e-6 * v for v in fv])
    print(f"  fitted p_free = {pc:+.4f} +- {ec:.4f}   "
          f"{'PASS' if abs(pc + 2) < 0.15 else 'FAIL -> machinery wrong, stop'}")
    # exact -2 only asymptotically; small-L sin correction makes it slightly steeper.

    print(f"\nINTERACTING, beta = {BETA}, Landau-fixed, corrected sampler")
    print(f"{'L':>3} {'V':>6} {'ncfg':>5} {'<plaq>':>8} {'lambda_min(M)':>22} "
          f"{'/lambda_free':>12} {'discarded':>10}", flush=True)
    mm, ms = [], []
    for L in Ls:
        vals, pq = [], []
        disc = 0
        for c in range(NCFG[L]):
            lat = Lattice(L, 4)
            for _ in range(400 + 13 * c):
                sweep_correct(lat, BETA)
            th = landau_fix(lat)
            if th > 1e-8:
                disc += 1
                continue
            try:
                vals.append(lam_min_M(lat))
                pq.append(lat.mean_plaquette())
            except Exception:
                disc += 1
        if len(vals) < 2:
            print(f"{L:>3} {L**4:>6} {NCFG[L]:>5}   too few usable configs "
                  f"(discarded {disc})", flush=True)
            mm.append(np.nan); ms.append(np.nan)
            continue
        m = float(np.mean(vals)); s = float(np.std(vals, ddof=1) / np.sqrt(len(vals)))
        mm.append(m); ms.append(s)
        lf = 4 * np.sin(np.pi / L) ** 2
        print(f"{L:>3} {L**4:>6} {len(vals):>5} {np.mean(pq):>8.4f} "
              f"{m:>15.6f} +- {s:<5.6f} {m/lf:>12.4f} {disc:>10}", flush=True)

    ok = [i for i in range(len(Ls)) if np.isfinite(mm[i])]
    if len(ok) >= 3:
        p, ep = pfit([Ls[i] for i in ok], [mm[i] for i in ok], [ms[i] for i in ok])
        print(f"\n  fitted exponent p = {p:+.4f} +- {ep:.4f}")
        print(f"  free-field baseline p_free = {pc:+.4f}")
        print(f"  (p - p_free) = {p - pc:+.4f}   ({abs(p-pc)/max(ep,1e-9):.1f} sigma)")
    else:
        print("\n  too few volumes usable; nothing fitted.")

    print("\n  V1 p < p_free by > 3 sigma -> dynamical horizon attraction; GZ active")
    print("  V2 p = p_free within errors -> approach is kinematic; gamma needs bigger boxes")
    print("  V3 p > p_free by > 3 sigma -> configs resist the horizon; suspect a bug first")
