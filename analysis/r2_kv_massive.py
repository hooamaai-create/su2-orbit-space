"""
P4(g) FIRST ATTACK — is the 0.0801 plateau a massive-screening number?

IDEA (pre-registered in RESEARCH2.md). The interacting K*V is volume-flat because the
interacting Delta^{-1} is short-ranged, while the massless free Delta^{-1} keeps
accumulating with volume. The cheapest quantitative version: evaluate the VALIDATED
free-field closed form (kv_constant.py, independently verified to 7 digits) with a MASSIVE
propagator, lambda(k) -> lambda(k) + m^2, and ask whether ONE m simultaneously:
   (i)  reproduces the measured plateau value 0.0801 at large L, and
   (ii) reproduces the measured flatness (interacting: 0.95% over L=3..6; free m=0: 7.4%).
Bonus: at m > 0 the k=0 mode contributes 1/m^2 — the U=1 discontinuity (K*V divergence)
is regulated by the mass, which is structurally the right cure for that singularity.

PRE-REGISTERED
 (M1) some m* gives BOTH |K*V(inf) - 0.0801| < 1% AND spread(L=8..16) < 1.5%
      -> the plateau IS an effective-screening number; report m* and compare it against
         the measured spectral scales of the interacting Delta (lambda_min(Delta) ~ 2.09
         at beta=2.4 in lattice units, i.e. m ~ 1.45). If m* ~ sqrt(lambda_min(Delta)),
         the "origin unknown" closes: 0.0801 = free-field K*V evaluated at the interacting
         Delta's own gap. A derivation with zero fitted inputs would then be the claim
         K*V_massive(m = sqrt(gap(Delta))) = plateau — testable at other betas.
 (M2) no single m does both -> the plateau is NOT massive screening; mechanism stays open.
Runs in seconds (pure FFT arithmetic). Prints numbers and the rule.
"""
import numpy as np

MEASURED_PLATEAU = 0.0801          # corrected-ensemble value, beta = 2.4, L = 3..6
MEASURED_GAP_DELTA = 2.09          # lambda_min(Delta), L=3, beta=2.4 (audit-verified)


def kv_massive(L, m2, D=4):
    """Free-field K*V with massive propagator: exact same derivation as kv_constant.py
    (validated to 7 digits) with lambda(k) -> lambda(k) + m^2 everywhere in Delta^{-1}.
    The transverse projector (gauge structure) is unchanged. Returns K*V (plane-average
    normalisation, the measured object: 18 V^2 S / (dimH (dimH - 1)))."""
    shape = (L,) * D
    V = L ** D
    k = 2.0 * np.pi * np.fft.fftfreq(L)
    grids = np.meshgrid(*([k] * D), indexing='ij')
    d = np.stack([np.exp(1j * g) - 1.0 for g in grids])
    lam = np.sum(np.abs(d) ** 2, axis=0)

    T = np.zeros((D, D) + shape, dtype=complex)
    for mu in range(D):
        T[mu, mu] = 1.0
    nz = lam > 1e-12
    for mu in range(D):
        for nu in range(D):
            T[mu, nu][nz] -= (d[mu] * np.conj(d[nu]))[nz] / lam[nz]

    dimH = float(np.real(np.trace(T.reshape(D, D, V), axis1=0, axis2=1).sum()))
    Tr = np.fft.ifftn(T, axes=tuple(range(2, 2 + D)))
    F = np.sum(np.real(Tr) ** 2 + np.imag(Tr) ** 2, axis=(0, 1))

    inv = 1.0 / (lam.ravel() + m2)                     # massive: k=0 included, finite
    Dp = np.real(np.fft.ifftn(inv.reshape(shape)))
    S = float(np.sum(Dp * F))
    return 18.0 * V ** 2 * S / (dimH * (dimH - 1.0))


if __name__ == "__main__":
    Ls = [4, 6, 8, 10, 12, 16]
    print("massive free-field K*V  (validated closed form + m^2)")
    print(f"{'m2':>7} {'m':>6} " + " ".join(f"L={L:<8}" for L in Ls) +
          f"{'spread L>=8':>12} {'KV(16)':>9}")
    best = None
    for m2 in (0.25, 0.5, 1.0, 1.5, 2.09, 2.5, 3.0, 4.0, 6.0):
        vals = [kv_massive(L, m2) for L in Ls]
        big = vals[2:]
        spread = 100 * (max(big) / min(big) - 1)
        row = " ".join(f"{v:.6f} " for v in vals)
        print(f"{m2:>7.2f} {np.sqrt(m2):>6.3f} {row}{spread:>11.2f}% {vals[-1]:>9.6f}",
              flush=True)
        if spread < 1.5:
            miss = abs(vals[-1] - MEASURED_PLATEAU) / MEASURED_PLATEAU
            if best is None or miss < best[1]:
                best = (m2, miss, vals[-1], spread)

    print(f"\n  measured plateau (interacting, beta=2.4): {MEASURED_PLATEAU}")
    print(f"  measured Delta gap: lambda_min(Delta) = {MEASURED_GAP_DELTA}"
          f"  ->  parameter-free candidate m2 = {MEASURED_GAP_DELTA}")
    pf = kv_massive(16, MEASURED_GAP_DELTA)
    print(f"  PARAMETER-FREE TEST: K*V_massive(m2 = 2.09, L=16) = {pf:.6f}"
          f"   vs 0.0801  ({100*(pf/MEASURED_PLATEAU-1):+.2f}%)")
    if best:
        m2s, miss, val, spread = best
        print(f"  best flat m2 = {m2s}: K*V(16) = {val:.6f}, miss {100*miss:.2f}%, "
              f"spread {spread:.2f}%")
    print("\n  M1 one m gives value <1% AND flat <1.5% -> plateau = screening number;")
    print("     compare m* to sqrt(lambda_min(Delta)) for the zero-parameter closure")
    print("  M2 no m does both -> not massive screening; mechanism stays open")
