"""
THE POSITIVE CONTROL THE DIMENSION PIPELINE HAS NEVER HAD.

THE GAP (identified by the adversarial audit; the most serious in the programme).
Every "dimensionless" verdict today -- Ric_I in two formulations, K*V, D=3, the whole KK
series, eight in all -- was produced by one pipeline:

      measure Q in lattice units across beta, fit d log Q / d beta,
      divide by d log a / d beta = -1/(8 b0) = -2.6917,  call the result a mass dimension.

That pipeline has NEVER been shown capable of returning a nonzero answer. At fixed L the
physical box shrinks ~11x across the beta range and the L^4 boxes cross the finite-T
deconfinement crossover. **A genuine mass-carrier, clamped by finite volume, would also
read n ~ 0.** Until the pipeline finds a dimension that is known to be there, "we measured
zero" and "there is no scale" are different sentences.

THE CONTROL.  The string tension sigma has mass dimension 2, so sigma a^2 in lattice units
must return n = 2. It is extracted from planar Wilson loops by the Creutz ratio

      chi(R,T) = - log [ W(R,T) W(R-1,T-1) / ( W(R-1,T) W(R,T-1) ) ]   ->  sigma a^2

Run through the IDENTICAL fitting machinery, on the SAME lattice sizes and beta windows used
for the physics claims, with the CORRECTED sampler (sampler_fix.sweep_correct, validated
against an independent Creutz heat bath to 0.3 sigma).

PRE-REGISTERED FALSIFIERS
 (F1) chi(2,2) returns n consistent with 2  ->  THE PIPELINE WORKS at these volumes. The
      eight "n ~ 0" verdicts then mean what I have been saying they mean, once remeasured
      on the corrected ensemble.
 (F2) chi(2,2) returns n ~ 0, or anything far from 2  ->  THE PIPELINE CANNOT DETECT A
      DIMENSION THAT IS THERE at these volumes. Every "dimensionless" claim in this ledger
      is then uninformative, not wrong but unsupported, and the programme's central negative
      result loses its evidential basis until redone on larger lattices.
 (F3) n is nonzero but not 2 -> report; quantify how much of the pipeline's dynamic range
      survives, and treat the eight nulls as bounded rather than clean.

This is the single most important measurement remaining, because F2 would retract more than
any individual result today. Prints numbers and the rule; no conclusion is printed.
"""
import numpy as np
from step3_lattice import Lattice, qmul, qconj
from sampler_fix import sweep_correct

B0 = 11 * 2 / (48 * np.pi ** 2)
DLOGA = -1.0 / (8 * B0)


def wilson_loop(lat, mu, nu, R, T):
    """Planar RxT Wilson loop in the (mu,nu) plane, averaged over all sites."""
    U = lat.U
    W = np.tile(np.array([1.0, 0, 0, 0]), (lat.V, 1))
    pos = np.arange(lat.V)
    for _ in range(R):                       # +mu
        W = qmul(W, U[mu][pos]); pos = lat.fwd[mu][pos]
    for _ in range(T):                       # +nu
        W = qmul(W, U[nu][pos]); pos = lat.fwd[nu][pos]
    for _ in range(R):                       # -mu
        pos = lat.bwd[mu][pos]; W = qmul(W, qconj(U[mu][pos]))
    for _ in range(T):                       # -nu
        pos = lat.bwd[nu][pos]; W = qmul(W, qconj(U[nu][pos]))
    return float(W[..., 0].mean())


def loops(lat, rmax=2):
    """W[R,T] averaged over all planes."""
    out = np.zeros((rmax + 1, rmax + 1))
    npl = 0
    for mu in range(lat.D):
        for nu in range(mu + 1, lat.D):
            npl += 1
            for R in range(1, rmax + 1):
                for T in range(1, rmax + 1):
                    out[R, T] += wilson_loop(lat, mu, nu, R, T)
    return out / npl


def creutz(W, R=2, T=2):
    num = W[R, T] * W[R - 1, T - 1]
    den = W[R - 1, T] * W[R, T - 1]
    if num <= 0 or den <= 0:
        return np.nan
    return -np.log(num / den)


def fit_dim(betas, vals, errs):
    b = np.array(betas, float); a = np.array(vals, float); e = np.array(errs, float)
    m = np.isfinite(a) & (a > 0)
    if m.sum() < 3:
        return None, None, int(m.sum())
    b, a, e = b[m], a[m], e[m]
    y = np.log(a); w = np.maximum(e / a, 1e-6); W = 1 / w ** 2
    Sw, Sx, Sy = W.sum(), (W * b).sum(), (W * y).sum()
    Sxx, Sxy = (W * b * b).sum(), (W * b * y).sum()
    den = Sw * Sxx - Sx ** 2
    sl = (Sw * Sxy - Sx * Sy) / den
    return sl / DLOGA, np.sqrt(Sw / den) / abs(DLOGA), int(m.sum())


if __name__ == "__main__":
    betas = [2.0, 2.2, 2.4, 2.6, 2.8]
    for L in (4, 6):
        print("=" * 78)
        print(f"POSITIVE CONTROL  L = {L}^4,  corrected sampler.  sigma a^2 must give n = 2")
        print("=" * 78)
        print(f"{'beta':>6} {'<plaq>':>9} {'W(2,2)':>10} {'chi(2,2)':>22}", flush=True)
        cm, cs = [], []
        for beta in betas:
            ch, pq, w22 = [], [], []
            for c in range(3):
                lat = Lattice(L, 4)
                for _ in range(500):
                    sweep_correct(lat, beta)
                acc = np.zeros((3, 3)); n = 0
                for _ in range(12):
                    for _ in range(3):
                        sweep_correct(lat, beta)
                    acc += loops(lat, 2); n += 1
                W = acc / n
                ch.append(creutz(W)); pq.append(lat.mean_plaquette()); w22.append(W[2, 2])
            ch = np.array(ch)
            cm.append(float(np.nanmean(ch)))
            cs.append(float(np.nanstd(ch, ddof=1) / np.sqrt(np.isfinite(ch).sum())))
            print(f"{beta:>6.1f} {np.mean(pq):>9.4f} {np.mean(w22):>10.5f} "
                  f"{cm[-1]:>15.5f} +- {cs[-1]:<6.5f}", flush=True)

        n, e, used = fit_dim(betas, cm, cs)
        if n is None:
            print(f"\n  only {used} usable points -> no fit. F3.")
        else:
            print(f"\n  MEASURED MASS DIMENSION of sigma a^2 = {n:+.4f} +- {e:.4f}"
                  f"   (using {used}/{len(betas)} points)")
            print(f"    distance from 2 : {abs(n-2)/e:.1f} sigma")
            print(f"    distance from 0 : {abs(n)/e:.1f} sigma")
        print()

    print("=" * 78)
    print("PRE-REGISTERED RULE")
    print("=" * 78)
    print("  F1 n ~ 2  -> pipeline WORKS at these volumes; the eight nulls are meaningful")
    print("  F2 n ~ 0  -> pipeline CANNOT see a dimension that is there; the eight nulls")
    print("              are unsupported, and the central negative result loses its basis")
    print("  F3 other  -> report; treat the nulls as bounded, not clean")
    print("\n  for reference, the eight metric verdicts were all |n| < 0.03")
