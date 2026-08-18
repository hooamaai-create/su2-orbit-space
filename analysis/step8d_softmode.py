"""
BOX 1d — is the T/C cancellation tied to the SOFT SPECTRUM of Delta, or accidental?

Both T and Ric are LINEAR in the explicit propagator G = Delta^{-1}:

    Ric(a,a) = 3 Sum eps eps a a  G^{cc'}(n,n')  P_H[...]          <- G appears once, linearly
    T(a,a)   = 3 Sum_n a Theta(n) a,   Theta from G(n,n)            <- likewise

so writing G = Sum_k psi_k psi_k^T / lambda_k gives an EXACT additive decomposition

    T = Sum_k T_k ,   Ric = Sum_k Ric_k ,   C = Sum_k C_k = Sum_k (T_k - Ric_k).

P_H is held EXACT throughout (it is a geometric projector; truncating it would stop it being
a projector and would confound the test).  Note P_V = Sum_k D psi_k psi_k^T D^dag / lambda_k
with ||D psi_k||^2 = lambda_k, so every term is a unit-norm rank-1 piece: P_V is finite no
matter how soft Delta gets.  The 1/lambda_k blow-up lives ONLY in the explicit G.

THE QUESTION (the reviewer's T_soft =? C_soft, made precise):
    for the softest modes, is  ||Ric_soft|| << ||T_soft||  ?
    i.e. does the soft-mode contribution cancel between T and C, leaving Ric untouched?

    yes  -> the cancellation is spectrally located and is a real mechanism
    no   -> the cancellation is accidental over the tested smearing range

FALSIFIER, fixed in advance: the soft-mode cancellation is claimed only if
||Ric_soft||/||T_soft|| < 0.25 for the softest mode AND that ratio falls as lambda_1 falls
(i.e. the cancellation gets BETTER as the mode gets softer). Anything else = accidental.

CAVEAT: n = 1 configuration per smearing level.
"""
import numpy as np
from step3_lattice import Lattice, EPS
from step6_scheme import geom, ape_smear


def build(lat):
    Delta, Dp, E = geom(lat)
    lam, psi = np.linalg.eigh(Delta)
    return Delta, Dp, E, lam, psi


def ops_from_G(lat, G, E):
    """Ric and T built from an arbitrary propagator matrix G (P_H kept exact via E)."""
    D, V = lat.D, lat.V
    dimH = E.shape[1]
    Ee = E.reshape(D, V, 3, dimH)
    Ric = np.empty((dimH, dimH))
    Iv = np.eye(dimH)
    for j in range(dimH):
        a = (E @ Iv[:, j]).reshape(D, V, 3)
        P = np.einsum('mnaj,mnbi,abc->nci', a[..., None], Ee, EPS, optimize=True)
        Q = (G @ P.reshape(3 * V, dimH)).reshape(V, 3, dimH)
        out = np.einsum('mnbi,nci,cab->mna', Ee, Q, EPS, optimize=True)
        Ric[:, j] = 3.0 * (E.T @ out.reshape(-1))
    Ric = 0.5 * (Ric + Ric.T)

    NL = D * 3 * V
    Tf = np.zeros((NL, NL))
    for n in range(V):
        Gn = G[3 * n:3 * n + 3, 3 * n:3 * n + 3]
        Th = 3.0 * (np.trace(Gn) * np.eye(3) - Gn.T)
        for mu in range(D):
            r = (mu * V + n) * 3
            Tf[r:r + 3, r:r + 3] = Th
    T = E.T @ Tf @ E
    T = 0.5 * (T + T.T)
    return Ric, T


def report(lat, label):
    Delta, Dp, E, lam, psi = build(lat)
    nz = lam > 1e-8
    lam_nz = lam[nz]
    print(f"\n--- {label}  <plaq>={lat.mean_plaquette():.4f}  "
          f"lam_1={lam_nz[0]:.5f}  lam_2={lam_nz[1]:.5f} ---")
    print(f"  {'modes kept':>14} {'||T_s||':>10} {'||Ric_s||':>10} {'||C_s||':>10} "
          f"{'Ric/T':>8} {'tr T_s/dim':>11} {'tr Ric_s/dim':>13}")
    idx = np.where(nz)[0]
    out = {}
    for tag, sel in (("softest 1", idx[:1]), ("softest 3", idx[:3]),
                     ("softest 10", idx[:10]), ("ALL", idx)):
        Gs = (psi[:, sel] / lam[sel]) @ psi[:, sel].T
        Ric_s, T_s = ops_from_G(lat, Gs, E)
        C_s = T_s - Ric_s
        nT, nR, nC = (np.linalg.norm(x) for x in (T_s, Ric_s, C_s))
        d = T_s.shape[0]
        print(f"  {tag:>14} {nT:>10.4f} {nR:>10.4f} {nC:>10.4f} {nR/nT:>8.4f} "
              f"{np.trace(T_s)/d:>11.5f} {np.trace(Ric_s)/d:>13.5f}")
        out[tag] = (nR / nT, lam_nz[0])
    return out


if __name__ == "__main__":
    res = {}
    base = Lattice(3, 4)
    for _ in range(500):
        base.metropolis_sweep(2.5)
    U0 = base.U.copy()
    for ns in (0, 2, 8):
        lat = Lattice(3, 4)
        lat.U = U0.copy()
        if ns:
            ape_smear(lat, 0.5, ns)
        res[ns] = report(lat, f"smear={ns}")

    print("\n================ VERDICT ================")
    print(f"  {'smear':>7} {'lambda_1':>10} {'||Ric||/||T|| for softest mode':>32}")
    for ns in (0, 2, 8):
        r, l1 = res[ns]["softest 1"]
        print(f"  {ns:>7} {l1:>10.5f} {r:>32.4f}")
    r0 = res[0]["softest 1"][0]; r8 = res[8]["softest 1"][0]
    print(f"\n  falls as lambda_1 falls?  {r0:.4f} -> {r8:.4f}  "
          f"({'YES' if r8 < r0 else 'NO'})")
    print(f"  below 0.25 at the softest point? {'YES' if r8 < 0.25 else 'NO'}")
    print("  (pre-registered rule: BOTH must hold to claim a spectral cancellation;")
    print("   otherwise the cancellation is accidental over the tested range.)")
    print("  CAVEAT: n = 1 configuration per smearing level.")
