"""
IS  K_median * V  A DERIVABLE CONSTANT?

CLAIM UNDER TEST.  Step 3 measured, over random horizontal 2-planes,
      K * V = 0.0793 / 0.0822 / 0.0794 / 0.0798      (L=2,3 x beta=2.0,2.4)
i.e. constant to ~3% across a 5.06x volume change.  Two candidate explanations:

  (H1)  K*V -> 1/(4 pi) = 0.0795775 ,  a continuum solid-angle constant.
  (H2)  K*V -> a LATTICE constant computable from the free-field transverse projector
        and the lattice scalar propagator.  Then it is DERIVABLE, carries no dynamics,
        and its agreement with 1/(4pi) (if any) is numerology.

DERIVATION FOR (H2).  Free field U = 1.  grad_mu w(n) = w(n+mu) - w(n) is colour-diagonal,
so in momentum space  d_mu(k) = e^{i k_mu} - 1 ,  lambda(k) = sum_mu |d_mu|^2 = Delta(k),
and the horizontal projector is the lattice transverse projector

      T_{mu nu}(k) = delta_{mu nu} - d_mu(k) conj(d_nu(k)) / lambda(k)      (k != 0)
      T(0) = identity                                                        (k  = 0)

Take x, y iid uniform on the unit sphere of H.  For large dim H their components have
covariance   E[x^a_mu(n) x^{a'}_nu(m)] = delta^{a a'} T_{mu nu}(n-m) / dimH.
With  rho^c(n) = eps^{abc} sum_mu x^a_mu(n) y^b_mu(n)  and  eps^{abc} eps^{abc'} = 2 delta^{cc'},

      E[rho^c(n) rho^{c'}(m)] = 2 delta^{cc'} F(n-m) / dimH^2 ,
      F(r) = sum_{mu nu} T_{mu nu}(r)^2 .

Hence   E[K] = 3 * E<rho, Delta^+ rho> = 18 V S / dimH^2 ,   S = sum_r Delta^+(r) F(r),

                    ***   K * V  =  18 V^2 S / dimH^2   ***

with  dimH = 3[ D + (D-1)(V-1) ]  exactly (grad has a 3-dim kernel: constants).
Everything on the right is a pure lattice number.  NO free parameter, NO fit.

PRE-REGISTERED FALSIFIERS
 (F1) free-field Monte-Carlo mean K*V disagrees with the analytic formula by > 2%
      => the derivation is wrong; report and stop, conclude nothing about H1/H2.
 (F2) analytic K*V agrees with free-field MC but the INTERACTING K*V (beta = 2.0-2.6,
      several configs) differs from it by more than its own error bar
      => K*V is NOT the free-field constant; it carries coupling dependence. H2 fails.
 (F3) analytic K*V is volume-DEPENDENT (drifts with L)
      => "K*V = const" was an accident of two nearby volumes; both H1 and H2 fail.
 (F4) analytic K*V is volume-independent and equals the interacting value
      => H2 CONFIRMED: K*V is a derivable free-field lattice constant.
      Then compare to 1/(4pi) as a SEPARATE, secondary question.

This script prints numbers and the pre-registered rule.  It does not print a conclusion.
"""
import numpy as np
from step3_lattice import Lattice, EPS

INV4PI = 1.0 / (4.0 * np.pi)


# --------------------------------------------------------------------------
# (A) ANALYTIC free-field K*V
# --------------------------------------------------------------------------
def kv_analytic(L, D=4):
    shape = (L,) * D
    V = L ** D
    k = 2.0 * np.pi * np.fft.fftfreq(L) * L / L          # k values 2 pi n / L
    grids = np.meshgrid(*([k] * D), indexing='ij')
    d = np.stack([np.exp(1j * g) - 1.0 for g in grids])   # (D, L,...,L)
    lam = np.sum(np.abs(d) ** 2, axis=0)                  # Delta(k)
    lam_flat = lam.ravel()

    # transverse projector T(k), with T(0) = I
    T = np.zeros((D, D) + shape, dtype=complex)
    for mu in range(D):
        T[mu, mu] = 1.0
    nz = lam > 1e-12
    for mu in range(D):
        for nu in range(D):
            T[mu, nu][nz] -= (d[mu] * np.conj(d[nu]))[nz] / lam[nz]

    # exact dim H  =  3 * sum_k tr T(k)
    dimH = 3.0 * float(np.real(np.trace(T.reshape(D, D, V), axis1=0, axis2=1).sum()))

    # T(r) = (1/V) sum_k e^{ikr} T(k)   ->  numpy ifft is exactly (1/V) sum_k e^{+ikr}
    Tr = np.fft.ifftn(T, axes=tuple(range(2, 2 + D)))
    F = np.sum(np.real(Tr) ** 2 + np.imag(Tr) ** 2, axis=(0, 1))   # |T_{mu nu}(r)|^2

    # Delta^+(r) = (1/V) sum_{k != 0} e^{ikr} / lambda(k)
    inv = np.zeros_like(lam_flat)
    inv[nz.ravel()] = 1.0 / lam_flat[nz.ravel()]
    Dp = np.real(np.fft.ifftn(inv.reshape(shape)))

    S = float(np.sum(Dp * F))
    return 18.0 * V ** 2 * S / dimH ** 2, dimH, S


# --------------------------------------------------------------------------
# (B) MONTE-CARLO K over random horizontal 2-planes, on any configuration
# --------------------------------------------------------------------------
def kv_mc(lat, nplanes=4000, rng=None):
    rng = rng or np.random.default_rng(0)
    G, Delta, Dp, E = lat.geometry()
    D, V, dimH = lat.D, lat.V, E.shape[1]
    c = rng.normal(size=(dimH, 2 * nplanes))
    X = (E @ c).reshape(D, V, 3, 2 * nplanes)
    Ks = np.empty(nplanes)
    for j in range(nplanes):
        x = X[..., 2 * j]; y = X[..., 2 * j + 1]
        x = x / np.linalg.norm(x)
        y = y - (x * y).sum() * x
        y = y / np.linalg.norm(y)
        r = np.cross(x, y).sum(axis=0).ravel()
        Ks[j] = 3.0 * r @ (Dp @ r)
    return Ks, dimH


if __name__ == "__main__":
    D = 4
    print("=" * 78)
    print("(A) ANALYTIC free-field prediction   K*V = 18 V^2 S / dimH^2")
    print("=" * 78)
    print(f"{'L':>3} {'V':>6} {'dimH':>8} {'S':>14} {'K*V (analytic)':>16} {'ratio to 1/4pi':>15}")
    ana = {}
    for L in (2, 3, 4, 5, 6, 8):
        kv, dimH, S = kv_analytic(L, D)
        ana[L] = kv
        print(f"{L:>3} {L**D:>6} {dimH:>8.0f} {S:>14.6e} {kv:>16.7f} {kv/INV4PI:>15.4f}")
    print(f"\n   1/(4 pi) = {INV4PI:.7f}")

    print("\n" + "=" * 78)
    print("(B) FREE-FIELD MONTE CARLO  (validates the derivation)   [F1: |diff| > 2% kills it]")
    print("=" * 78)
    rng = np.random.default_rng(7)
    print(f"{'L':>3} {'dimH_mc':>8} {'dimH_ana':>9} {'mean K*V':>12} {'median K*V':>12} "
          f"{'analytic':>11} {'mean/ana':>9}")
    for L in (2, 3):
        lat = Lattice(L, D)                      # cold start = free field U = 1
        Ks, dimH = kv_mc(lat, nplanes=3000, rng=rng)
        kvm = Ks.mean() * lat.V
        kvmed = np.median(Ks) * lat.V
        print(f"{L:>3} {dimH:>8} {kv_analytic(L, D)[1]:>9.0f} {kvm:>12.6f} {kvmed:>12.6f} "
              f"{ana[L]:>11.6f} {kvm/ana[L]:>9.4f}")

    print("\n" + "=" * 78)
    print("(C) INTERACTING  (several configs per point)   [F2: differs from analytic]")
    print("=" * 78)
    print(f"{'L':>3} {'beta':>5} {'n':>3} {'<plaq>':>8} {'mean K*V':>20} {'median K*V':>20}")
    inter = {}
    for L, beta in [(2, 2.0), (2, 2.4), (3, 2.0), (3, 2.4), (3, 2.6)]:
        means, meds, plaqs = [], [], []
        for cfg in range(4):
            lat = Lattice(L, D)
            for _ in range(300):
                lat.metropolis_sweep(beta)
            Ks, _ = kv_mc(lat, nplanes=2000, rng=np.random.default_rng(100 + cfg))
            means.append(Ks.mean() * lat.V)
            meds.append(np.median(Ks) * lat.V)
            plaqs.append(lat.mean_plaquette())
        m, s = np.mean(means), np.std(means, ddof=1) / np.sqrt(len(means))
        md, sd = np.mean(meds), np.std(meds, ddof=1) / np.sqrt(len(meds))
        inter[(L, beta)] = (m, s)
        print(f"{L:>3} {beta:>5.1f} {len(means):>3} {np.mean(plaqs):>8.4f} "
              f"{m:>13.6f} +- {s:<6.6f} {md:>13.6f} +- {sd:<6.6f}")

    print("\n" + "=" * 78)
    print("PRE-REGISTERED RULE (state which falsifier fires; do not editorialise)")
    print("=" * 78)
    print("  F1  free-field MC vs analytic differs > 2%      -> derivation wrong, stop")
    print("  F2  interacting differs from analytic > 1 sigma -> K*V is not the free constant")
    print("  F3  analytic K*V drifts with L                  -> 'K*V = const' was accidental")
    print("  F4  analytic flat AND matches interacting       -> K*V is a derivable constant")
    print(f"\n  analytic spread over L=2..8: "
          f"{min(ana.values()):.6f} .. {max(ana.values()):.6f}  "
          f"({100*(max(ana.values())/min(ana.values())-1):.2f}% )")
    for key, (m, s) in inter.items():
        L = key[0]
        print(f"  L={key[0]} beta={key[1]}: interacting {m:.6f}+-{s:.6f}   "
              f"analytic {ana[L]:.6f}   diff {(m-ana[L])/max(s,1e-12):+.1f} sigma")
