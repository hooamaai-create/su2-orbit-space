"""
IS THE GEOMETRY BLIND TO THE GRIBOV HORIZON?  Direct test on configurations.

THE CHAIN BEING TESTED (derived 2026-08-12 from the user's observer/Unruh observation).
  L1  Any relation among quantities that scale together under a -> lambda a has the scale
      cancel. Pure YM has ONE scale, so every internal relation yields ratios. Confirmed
      four times today (KK, hbar c/L, transmutation, Unruh/CKS). To mint a number the theory
      must be FORCED onto a saturated condition.
  L2  Orbit space HAS one: Zwanziger's horizon condition, an EQUALITY stating the measure
      sits exactly on the boundary where M = -div D loses positivity. It generates a mass,
      the Gribov parameter gamma, entering D(k) = k^2/(k^4 + gamma^4).
  L3  But the Gram/Schur result (LEDGER 2026-08-10) gives Delta >= M^dag L^-1 M, so as
      M -> 0 at the horizon, Delta = D^dag D stays STRICTLY POSITIVE. The metric does not
      degenerate at the boundary; Ric, built from Delta and P_H alone, never learns it is
      there.
  L4  Therefore the scale lives at the horizon and the curvature is blind to it, so the
      curvature CANNOT carry the scale -- which DERIVES the measured |n| <= 0.024 instead
      of merely observing it.

WHAT THIS SCRIPT TESTS.  L3, directly, on real configurations rather than from an inequality.
Both operators live on the SAME gauge field; nothing a priori forbids them moving together.
Fix Landau gauge (maximise sum Re tr U, the gauge in which M = -div D is the FP operator),
then on each configuration measure:
      lambda_min(M)   -- distance to the Gribov horizon (small = near the boundary)
      lambda_min(Delta), K*V  -- the metric quantities
and correlate them ACROSS configurations at fixed beta and volume.

PRE-REGISTERED FALSIFIERS
 (H1) |Pearson r| < 0.3 between lambda_min(M) and each metric quantity, over >= 20 configs
      -> the metric is blind to horizon proximity. L3 CONFIRMED on configurations, and the
         chain L1-L4 stands as a derivation of the day's measured null.
 (H2) |r| > 0.6 -> the metric DOES track horizon proximity. L3 is wrong as a statement about
      configurations even though the inequality is proved, and the chain's conclusion does
      not follow. RETRACT the derivation.
 (H3) 0.3 <= |r| <= 0.6 -> inconclusive; report r and the sample size, claim nothing.

CONTROL. lambda_min(M) must actually VARY across configurations (spread > 20% of its mean),
otherwise there is no signal to correlate against and the test is vacuous. Also: Landau
gauge fixing must converge (theta < 1e-8) on every configuration used; non-converged ones
are DISCARDED and counted, not silently kept.
"""
import faulthandler
faulthandler.enable()

import numpy as np
from step3_lattice import Lattice, qmul, qconj, adjoint_matrix, EPS
from sampler_fix import sweep_correct
from mf_kv import kv_mf, ops
from scipy.sparse.linalg import LinearOperator, eigsh


def landau_fix(lat, tol=1e-9, itmax=6000, alpha=0.08):
    """Maximise F = sum_{n,mu} Re tr U_mu(n) by steepest descent.  Returns final theta."""
    V, D = lat.V, lat.D
    for _ in range(itmax):
        # g(n) = sum_mu [ A_mu(n) - A_mu(n-mu) ],  A_mu = traceless antihermitian part
        g = np.zeros((V, 3))
        for mu in range(D):
            A = lat.U[mu][:, 1:]                 # quaternion vector part ~ A_mu
            g += A - A[lat.bwd[mu]]
        theta = float((g ** 2).sum() / V)
        if theta < tol:
            return theta
        v = -alpha * g
        n = np.linalg.norm(v, axis=-1, keepdims=True)
        w = np.cos(n)
        s = np.where(n > 1e-14, np.sin(n) / np.maximum(n, 1e-14), 1.0)
        r = np.concatenate([w, s * v], axis=-1)
        for mu in range(D):
            lat.U[mu] = qmul(qmul(r, lat.U[mu]), qconj(r[lat.fwd[mu]]))
    return theta


def _grad_from_U(U, bwd, D):
    """grad of the Landau functional F = sum_{n,mu} Re tr U_mu(n), from raw link arrays.
    Takes U/bwd directly so the Hessian matvec never rebuilds a Lattice (which would
    regenerate every neighbour-index array on each of several hundred matvecs)."""
    g = np.zeros(U[0][:, 1:].shape)
    for mu in range(D):
        A = U[mu][:, 1:]
        g += A - A[bwd[mu]]
    return g


def gauge_grad(lat):
    """Stationary points of F are exactly the Landau-gauge configurations."""
    return _grad_from_U(lat.U, lat.bwd, lat.D)


def _transform_U(U, fwd, D, v):
    """Links gauge-transformed by g = exp(v); returns the array only, no Lattice object."""
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    w = np.cos(n)
    s = np.where(n > 1e-14, np.sin(n) / np.maximum(n, 1e-14), 1.0)
    r = np.concatenate([w, s * v], axis=-1)
    return np.stack([qmul(qmul(r, U[mu]), qconj(r[fwd[mu]])) for mu in range(D)])


def gauge_transform(lat, v):
    """Return a copy of lat gauge-transformed by g = exp(v), v shaped (V,3)."""
    new = Lattice(lat.L, lat.D)
    new.U = _transform_U(lat.U, lat.fwd, lat.D, v)
    return new


def M_operator(lat, h=1e-4):
    """Faddeev-Popov operator as the HESSIAN of the gauge functional, taken as a central
    difference of the gradient.

    A hand-derived M was WRONG (2026-08-12): ||M - M^T||/||M|| = 0.44, a 44% antisymmetric
    part, which produced spurious negative eigenvalues on configurations that are local
    maxima of F.  The FP operator is the second derivative of a real function and is
    therefore symmetric by construction; taking it as a difference of the gradient makes
    that automatic instead of relying on my algebra.  Same technique as bochner_hessV.py.

    Sign: at a local MAXIMUM of F the Hessian is negative semi-definite, so M = -Hessian
    is positive semi-definite -- the Gribov region.
    """
    V = lat.V

    def mv(x):
        # gauge_grad returns -dF/da, so d(gauge_grad)/da = -Hessian(F) = M.  Verified:
        # at a local maximum of F this comes out positive semi-definite with exactly 3
        # zero modes (the global colour rotations), i.e. inside the Gribov region.
        v = x.reshape(V, 3)
        gp = gauge_grad(gauge_transform(lat, +h * v))
        gm = gauge_grad(gauge_transform(lat, -h * v))
        return ((gp - gm) / (2 * h)).ravel()
    return LinearOperator((3 * V, 3 * V), matvec=mv, dtype=float)


def lam_min_M(lat, k=6):
    """Smallest eigenvalue of M above the trivial constant modes."""
    A = M_operator(lat)
    ev = np.sort(eigsh(A, k=k, which='SA', tol=1e-7, maxiter=20000,
                       return_eigenvectors=False))
    return float(ev[3])          # drop the 3 global colour zero modes


def ric_min(lat):
    G, Delta, Dp, E = lat.geometry()
    D, V, dimH = lat.D, lat.V, E.shape[1]
    Ee = E.reshape(D, V, 3, dimH)

    def ap(vec):
        a = (E @ vec).reshape(D, V, 3)
        P = np.einsum('mnaj,mnbi,abc->nci', a[..., None], Ee, EPS, optimize=True)
        Q = (Dp @ P.reshape(3 * V, dimH)).reshape(V, 3, dimH)
        out = np.einsum('mnbi,nci,cab->mna', Ee, Q, EPS, optimize=True)
        return 3.0 * (E.T @ out.reshape(-1))
    op = LinearOperator((dimH, dimH), matvec=ap, dtype=float)
    return float(np.min(eigsh(op, k=min(4, dimH - 2), which='SA', tol=1e-8,
                              maxiter=8000, return_eigenvectors=False)))


def pearson(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    x = x - x.mean(); y = y - y.mean()
    d = np.sqrt((x * x).sum() * (y * y).sum())
    return float((x * y).sum() / d) if d > 0 else np.nan


if __name__ == "__main__":
    L, BETA, NCFG = 4, 2.4, 24
    print(f"Landau-gauge-fixed SU(2), L={L}^4, beta={BETA}, target {NCFG} configs")
    print(f"{'cfg':>4} {'theta':>10} {'lam_min(M)':>12} {'lam_min(Delta)':>14} "
          f"{'Ric_min':>10} {'K*V':>10}", flush=True)

    lm, ld, rc, kv = [], [], [], []
    discarded = 0
    for c in range(NCFG):
        lat = Lattice(L, 4)
        for _ in range(400 + 7 * c):            # decorrelate configurations
            sweep_correct(lat, BETA)
        th = landau_fix(lat)
        if th > 1e-8:
            discarded += 1
            print(f"{c:>4} {th:>10.2e}   DISCARDED (gauge fixing did not converge)",
                  flush=True)
            continue
        try:
            mm = lam_min_M(lat)
            _, Delta, _, _ = lat.geometry()
            dd = float(np.sort(np.linalg.eigvalsh(Delta))[3])
            rr = ric_min(lat)
            kk = kv_mf(lat, 150, seed=7000 + c).mean() * lat.V
        except Exception as e:
            discarded += 1
            print(f"{c:>4} {th:>10.2e}   DISCARDED ({type(e).__name__})", flush=True)
            continue
        lm.append(mm); ld.append(dd); rc.append(rr); kv.append(kk)
        print(f"{c:>4} {th:>10.2e} {mm:>12.6f} {dd:>14.6f} {rr:>10.6f} {kk:>10.6f}",
              flush=True)

    n = len(lm)
    print(f"\n  usable configs: {n}   discarded: {discarded}")
    if n < 8:
        print("  too few configurations to correlate. Test inconclusive; nothing claimed.")
        raise SystemExit(0)

    lm = np.array(lm)
    spread = (lm.max() - lm.min()) / lm.mean()
    print(f"\n  CONTROL: lambda_min(M) spread = {100*spread:.1f}% of its mean "
          f"(needs > 20% for the test to be non-vacuous)  "
          f"{'PASS' if spread > 0.20 else 'FAIL -> test is vacuous'}")
    print(f"    lambda_min(M): min {lm.min():.6f}  mean {lm.mean():.6f}  max {lm.max():.6f}")

    print(f"\n  CORRELATION of horizon proximity with the metric quantities:")
    for nm, arr in (("lambda_min(Delta)", ld), ("Ric_min", rc), ("K*V", kv)):
        r = pearson(lm, arr)
        print(f"    r( lambda_min(M) , {nm:<18} ) = {r:+.4f}")

    rs = [abs(pearson(lm, a)) for a in (ld, rc, kv)]
    mx = max(rs)
    print(f"\n  max |r| = {mx:.4f}   (n = {n})")
    print("\n  H1 |r| < 0.3      -> metric blind to the horizon; chain L1-L4 stands")
    print("  H2 |r| > 0.6      -> metric TRACKS the horizon; RETRACT the derivation")
    print("  H3 0.3 <= |r| <= 0.6 -> inconclusive; report r and n, claim nothing")
