"""
STEP 6 — CONTINUUM LIMIT / scheme independence.

Step 5 showed the un-renormalised geometric gap has NO finite continuum limit. Mondal's
construction survives only if the divergence CANCELS between the two sides of

        Delta_E  >=  (hbar^2/2) * Delta_Ricci

his stated reason being that "the singular parts of Delta_E and Delta vary uniformly as one
changes the cut-off parameter chi". That is a claim of SCHEME INDEPENDENCE of the ratio.

TEST. Change the UV regulator without changing the IR, by APE smearing the links, and ask
whether the dimensionless ratio

        R  =  lambda_min(Ric_I)  /  lambda_min(Delta)

is invariant. Smearing is an ideal probe here: it is a pure short-distance filter with a
tunable radius, it leaves long-distance physics intact, and for SU(2) the projection back
into the group is exact (any real combination of quaternions is a multiple of a unit one).

If both quantities are IR-dominated, R is stable. If they are UV-dominated by DIFFERENT
amounts, R drifts -- and "the singular parts are exactly the same" fails for this pair.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh
from step3_lattice import Lattice, adjoint_matrix, qmul, qconj, EPS

rng = np.random.default_rng(41)


def ape_smear(lat, alpha=0.5, n=1):
    """APE smearing; SU(2) projection = quaternion normalisation (exact for SU(2)).

    NOTE ON ORIENTATION (bug found 2026-08-10 20:58): this code's staple_sum S is defined
    so that the action is -beta * w(U.S), i.e. plaquette = Re Tr(U S) = Re Tr(U V^dag) with
    V the conventional staple. So V = S^dag, and APE must add qconj(S), NOT S. Adding S
    anti-aligns the links: <plaq> went 0.612 -> 0.025 -> -0.220 instead of rising toward 1.
    """
    for _ in range(n):
        new = np.empty_like(lat.U)
        for mu in range(lat.D):
            V = qconj(lat.staple_sum(mu))
            q = (1 - alpha) * lat.U[mu] + (alpha / (2 * (lat.D - 1))) * V
            nrm = np.linalg.norm(q, axis=-1, keepdims=True)
            nrm = np.where(nrm < 1e-12, 1.0, nrm)
            new[mu] = q / nrm
        lat.U = new
    return lat


def _smear_sanity(L=4, beta=2.5, sweeps=200, alpha=0.5, nmax=8):
    """Guard: smearing must drive <plaq> monotonically toward 1."""
    lat = Lattice(L, 4)
    for _ in range(sweeps):
        lat.metropolis_sweep(beta)
    seq = [lat.mean_plaquette()]
    for _ in range(nmax):
        ape_smear(lat, alpha, 1)
        seq.append(lat.mean_plaquette())
    ok = all(seq[i + 1] >= seq[i] - 1e-9 for i in range(len(seq) - 1)) and seq[-1] > seq[0]
    print("  smearing sanity <plaq>: " + " -> ".join(f"{v:.4f}" for v in seq))
    print(f"  monotone increasing toward 1: {'PASS' if ok else 'FAIL'}")
    return ok


def geom(lat):
    D, V = lat.D, lat.V
    R = np.stack([adjoint_matrix(lat.U[mu]) for mu in range(D)])
    NV, NL = 3 * V, D * 3 * V
    G = np.zeros((NL, NV))
    for mu in range(D):
        f = lat.fwd[mu]
        for n in range(V):
            r0 = (mu * V + n) * 3
            G[r0:r0+3, 3*f[n]:3*f[n]+3] += R[mu, n]
            G[r0:r0+3, 3*n:3*n+3] -= np.eye(3)
    Delta = G.T @ G
    Dp = np.linalg.pinv(Delta, rcond=1e-10)
    PH = np.eye(NL) - G @ Dp @ G.T
    w, Vv = np.linalg.eigh(PH)
    E = Vv[:, w > 0.5]
    return Delta, Dp, E


def ric_min(lat, Dp, E):
    D, V = lat.D, lat.V
    dimH = E.shape[1]
    Ee = E.reshape(D, V, 3, dimH)

    def mv(vec):
        a = (E @ vec).reshape(D, V, 3)
        P = np.einsum('mnaj,mnbi,abc->nci', a[..., None], Ee, EPS, optimize=True)
        Q = (Dp @ P.reshape(3 * V, dimH)).reshape(V, 3, dimH)
        out = np.einsum('mnbi,nci,cab->mna', Ee, Q, EPS, optimize=True)
        return 3.0 * (E.T @ out.reshape(-1))

    # k=1 does not converge for 'SA' on this operator; k=6 (as in step3) does.
    op = LinearOperator((dimH, dimH), matvec=mv, dtype=float)
    k = min(6, dimH - 2)
    lo = eigsh(op, k=k, which='SA', tol=1e-6, maxiter=20000, return_eigenvectors=False)
    hi = eigsh(op, k=2, which='LA', tol=1e-6, maxiter=20000, return_eigenvectors=False)
    return float(np.min(lo)), float(np.max(hi))


def run(L, beta, sweeps=500, alpha=0.5, levels=(0, 1, 2, 4, 8)):
    base = Lattice(L, 4)
    for _ in range(sweeps):
        base.metropolis_sweep(beta)
    U0 = base.U.copy()
    print(f"\n=== SU(2) L={L}^4  beta={beta}   APE alpha={alpha} ===")
    print(f"{'n_smear':>8} {'<plaq>':>8} {'lam_min(Delta)':>15} {'Ric_min':>10} "
          f"{'Ric_max':>10} {'R=Ric_min/lamD':>15}")
    out = []
    for ns in levels:
        lat = Lattice(L, 4)
        lat.U = U0.copy()
        if ns:
            ape_smear(lat, alpha, ns)
        Delta, Dp, E = geom(lat)
        lamD = float(np.linalg.eigvalsh(Delta)[0])
        rmin, rmax = ric_min(lat, Dp, E)
        R = rmin / lamD
        out.append((ns, lat.mean_plaquette(), lamD, rmin, rmax, R))
        print(f"{ns:>8} {lat.mean_plaquette():>8.4f} {lamD:>15.5f} {rmin:>10.5f} "
              f"{rmax:>10.5f} {R:>15.5f}")
    ns = np.array([o[0] for o in out], float)
    lamD = np.array([o[2] for o in out])
    rmin = np.array([o[3] for o in out])
    Rr = np.array([o[5] for o in out])
    print(f"  change from n=0 to n={int(ns[-1])}:")
    print(f"    lambda_min(Delta): {lamD[0]:.5f} -> {lamD[-1]:.5f}   "
          f"factor {lamD[-1]/lamD[0]:.4f}")
    print(f"    Ric_min          : {rmin[0]:.5f} -> {rmin[-1]:.5f}   "
          f"factor {rmin[-1]/rmin[0]:.4f}")
    print(f"    RATIO R          : {Rr[0]:.5f} -> {Rr[-1]:.5f}   "
          f"factor {Rr[-1]/Rr[0]:.4f}   drift {100*abs(Rr[-1]/Rr[0]-1):.1f}%")
    return out


if __name__ == "__main__":
    run(4, 2.5)
    run(6, 2.5)
