"""
STEP 3 — TEST ON LATTICE.

Measure orbit-space curvature on REAL SU(2) gauge configurations, with the emphasis the
Step-2 audit demands: the LOWER tail, not the upper. Theorem 1.1 of arXiv:2301.06996 needs

        Ric^{B.E}(alpha,alpha) >= Delta * G(alpha,alpha)   uniformly,  Delta > 0

so the object of interest is  lambda_min( Ric ) , not <K> and not K_max.

WHAT IS COMPUTED
  - SU(2) Wilson lattice gauge theory on L^4, Metropolis, thermalised.
  - Adjoint covariant difference   (grad w)_{mu}(n) = R_mu(n) w(n+mu) - w(n),  R = Ad(U).
  - FP / vertical Laplacian        Delta = grad^T grad   (positive semi-definite).
  - Horizontal space               H = ker(grad^T)   (= Landau/Coulomb-type transversality).
  - Sectional curvature            K(x,y) = 3 <rho, Delta^+ rho> / (|x|^2|y|^2 - <x,y>^2),
                                   rho(x,y)(n) = sum_mu x_mu(n) X y_mu(n)   (cross product).
  - Term-I Ricci form              Ric(a,a) = 3 sum_i <rho(a,e_i), Delta^+ rho(a,e_i)>
                                   over an orthonormal basis {e_i} of H, matrix-free.

SCOPE / HONESTY
  This is the CONTINUUM-FORM estimator evaluated with lattice-regularised operators. The
  exact curvature of the compact lattice orbit space G^{|L|}/G^{|V|} carries an ADDITIONAL
  non-negative term from the bi-invariant curvature of the total space. Hence what is
  computed here is a LOWER BOUND on the exact lattice sectional curvature -- which is the
  correct side for testing a lower-bound hypothesis.
  All quantities are in lattice units (a = 1); no continuum extrapolation is claimed here.
"""
import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

rng = np.random.default_rng(11)

# ---------------- SU(2) as unit quaternions ----------------
def qmul(p, q):
    w1, x1, y1, z1 = p[..., 0], p[..., 1], p[..., 2], p[..., 3]
    w2, x2, y2, z2 = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    return np.stack([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                     w1*x2 + x1*w2 + y1*z2 - z1*y2,
                     w1*y2 - x1*z2 + y1*w2 + z1*x2,
                     w1*z2 + x1*y2 - y1*x2 + z1*w2], axis=-1)

def qconj(p):
    return np.stack([p[..., 0], -p[..., 1], -p[..., 2], -p[..., 3]], axis=-1)

def qrand_near_identity(shape, eps):
    v = rng.normal(size=shape + (3,)) * eps
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    w = np.cos(n)
    s = np.where(n > 1e-12, np.sin(n) / np.maximum(n, 1e-12), 1.0)
    return np.concatenate([w, s * v], axis=-1)

def adjoint_matrix(q):
    """SO(3) matrix of the adjoint action of a unit quaternion (batched)."""
    w, x, y, z = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    R = np.empty(q.shape[:-1] + (3, 3))
    R[..., 0, 0] = 1 - 2*(y*y + z*z); R[..., 0, 1] = 2*(x*y - z*w); R[..., 0, 2] = 2*(x*z + y*w)
    R[..., 1, 0] = 2*(x*y + z*w); R[..., 1, 1] = 1 - 2*(x*x + z*z); R[..., 1, 2] = 2*(y*z - x*w)
    R[..., 2, 0] = 2*(x*z - y*w); R[..., 2, 1] = 2*(y*z + x*w); R[..., 2, 2] = 1 - 2*(x*x + y*y)
    return R


class Lattice:
    def __init__(self, L, D=4):
        self.L, self.D, self.V = L, D, L ** D
        self.shape = (L,) * D
        idx = np.arange(self.V).reshape(self.shape)
        self.fwd = [np.roll(idx, -1, axis=mu).ravel() for mu in range(D)]   # n -> n+mu
        self.bwd = [np.roll(idx, +1, axis=mu).ravel() for mu in range(D)]   # n -> n-mu
        self.U = np.zeros((D, self.V, 4)); self.U[..., 0] = 1.0            # cold start

    def plaquette_action(self, beta):
        """S = -beta * sum_p (1/2)Tr U_p = -beta * sum_p w(U_p)  [SU(2): (1/2)Tr = w]."""
        tot = 0.0
        for mu in range(self.D):
            for nu in range(mu + 1, self.D):
                P = qmul(qmul(self.U[mu], self.U[nu][self.fwd[mu]]),
                         qmul(qconj(self.U[mu][self.fwd[nu]]), qconj(self.U[nu])))
                tot += P[..., 0].sum()
        return -beta * tot

    def staple_sum(self, mu):
        """Sum of staples for link mu at every site, as a quaternion-valued (V,4) array."""
        S = np.zeros((self.V, 4))
        for nu in range(self.D):
            if nu == mu:
                continue
            # forward staple: U_nu(n+mu) U_mu(n+nu)^dag U_nu(n)^dag
            S += qmul(qmul(self.U[nu][self.fwd[mu]], qconj(self.U[mu][self.fwd[nu]])),
                      qconj(self.U[nu]))
            # backward staple: U_nu(n+mu-nu)^dag U_mu(n-nu)^dag U_nu(n-nu)
            S += qmul(qmul(qconj(self.U[nu][self.bwd[nu]][self.fwd[mu]]),
                           qconj(self.U[mu][self.bwd[nu]])), self.U[nu][self.bwd[nu]])
        return S

    def metropolis_sweep(self, beta, eps=0.35, hits=3):
        acc = 0; tot = 0
        for mu in range(self.D):
            S = self.staple_sum(mu)
            for _ in range(hits):
                r = qrand_near_identity((self.V,), eps)
                Unew = qmul(r, self.U[mu])
                # dS = -beta * ( w(Unew S) - w(Uold S) )
                dS = -beta * (qmul(Unew, S)[..., 0] - qmul(self.U[mu], S)[..., 0])
                take = (dS <= 0) | (rng.random(self.V) < np.exp(-np.clip(dS, 0, 50)))
                self.U[mu] = np.where(take[:, None], Unew, self.U[mu])
                acc += take.sum(); tot += self.V
        return acc / tot

    def mean_plaquette(self):
        tot = 0.0; c = 0
        for mu in range(self.D):
            for nu in range(mu + 1, self.D):
                P = qmul(qmul(self.U[mu], self.U[nu][self.fwd[mu]]),
                         qmul(qconj(self.U[mu][self.fwd[nu]]), qconj(self.U[nu])))
                tot += P[..., 0].sum(); c += self.V
        return tot / c

    # ---------------- orbit-space geometry ----------------
    def geometry(self):
        D, V = self.D, self.V
        R = np.stack([adjoint_matrix(self.U[mu]) for mu in range(D)])       # (D,V,3,3)
        # grad : R^{3V} -> R^{D*3V}     (grad w)_{mu,n} = R_mu(n) w(n+mu) - w(n)
        NV, NL = 3 * V, D * 3 * V
        G = np.zeros((NL, NV))
        for mu in range(D):
            for n in range(V):
                r0 = (mu * V + n) * 3
                G[r0:r0+3, 3*self.fwd[mu][n]:3*self.fwd[mu][n]+3] += R[mu, n]
                G[r0:r0+3, 3*n:3*n+3] -= np.eye(3)
        Delta = G.T @ G
        # horizontal projector = I - G (G^T G)^+ G^T
        Dp = np.linalg.pinv(Delta, rcond=1e-10)
        PH = np.eye(NL) - G @ Dp @ G.T
        w, Vv = np.linalg.eigh(PH)
        E = Vv[:, w > 0.5]                                                 # (NL, dimH)
        return G, Delta, Dp, E

    @staticmethod
    def rho(x, y, D, V):
        """rho(n) = sum_mu x_mu(n) X y_mu(n);  x,y shaped (D,V,3) -> (V,3)."""
        return np.cross(x, y).sum(axis=0)


def analyse(L, beta, D=4, nsweep=200, nplanes=4000, label=""):
    lat = Lattice(L, D)
    for i in range(nsweep):
        lat.metropolis_sweep(beta)
    plaq = lat.mean_plaquette()
    G, Delta, Dp, E = lat.geometry()
    V, dimH = lat.V, E.shape[1]
    ev = np.linalg.eigvalsh(Delta)
    print(f"\n=== {label} L={L}^{D}, beta={beta} ===")
    print(f"  <(1/2)Tr U_plaq> = {plaq:.4f}   dim A = {D*3*V}, dim orbit = {3*V}, dim H = {dimH}")
    print(f"  Delta spectrum: lambda_min = {ev[0]:.3e} (exact zero modes expected: 0 for "
          f"generic A up to global), next = {ev[1]:.4f}, lambda_max = {ev[-1]:.3f}")

    # ---- sectional curvature over random horizontal 2-planes ----
    c = rng.normal(size=(dimH, 2 * nplanes))
    Xall = (E @ c).reshape(D, V, 3, 2 * nplanes)
    Ks = np.empty(nplanes)
    for k in range(nplanes):
        x = Xall[..., 2*k]; y = Xall[..., 2*k+1]
        nx = np.linalg.norm(x); ny = np.linalg.norm(y)
        x = x / nx; y = y / ny
        y = y - (x*y).sum() * x
        y /= np.linalg.norm(y)
        r = Lattice.rho(x, y, D, V).ravel()
        Ks[k] = 3.0 * r @ (Dp @ r)
    qs = np.percentile(Ks, [0, 1, 5, 25, 50, 75, 95, 100])
    print(f"  sectional K over {nplanes} random horizontal 2-planes (lattice units):")
    print(f"    min={qs[0]:.3e}  1%={qs[1]:.3e}  5%={qs[2]:.3e}  median={qs[4]:.3e}  "
          f"95%={qs[6]:.3e}  max={qs[7]:.3e}")
    print(f"    mean={Ks.mean():.3e}   fraction with K < 1% of median = "
          f"{(Ks < 0.01*qs[4]).mean():.3f}")

    # ---- term-I Ricci form, matrix-free lambda_min ----
    Ee = E.reshape(D, V, 3, dimH)

    def ric_apply(vec):
        a = (E @ vec).reshape(D, V, 3)
        # P[n,c,i] = sum_mu ( a_mu(n) X e_{i,mu}(n) )_c
        P = np.einsum('mnaj,mnbi,abc->nci', a[..., None], Ee, EPS, optimize=True)
        P = P.reshape(3*V, dimH)
        Q = Dp @ P                                   # (3V, dimH)
        Q = Q.reshape(V, 3, dimH)
        # out_{mu,n,a} = sum_i eps_{cab} e_{i,mu,n,b} Q_{n,c,i}
        out = np.einsum('mnbi,nci,cab->mna', Ee, Q, EPS, optimize=True)
        return 3.0 * (E.T @ out.reshape(-1))

    op = LinearOperator((dimH, dimH), matvec=ric_apply, dtype=float)
    k = min(6, dimH - 2)
    lo = eigsh(op, k=k, which='SA', tol=1e-8, maxiter=5000, return_eigenvectors=False)
    hi = eigsh(op, k=2, which='LA', tol=1e-8, maxiter=5000, return_eigenvectors=False)
    lo = np.sort(lo); hi = np.sort(hi)
    print(f"  term-I Ricci form on H (this is Mondal's Delta, kinetic part, lattice units):")
    print(f"    lambda_min  = {lo[0]:.6e}      lowest {k}: {np.array2string(lo, precision=4)}")
    print(f"    lambda_max  = {hi[-1]:.6e}")
    print(f"    ratio lambda_min/lambda_max = {lo[0]/hi[-1]:.3e}")
    return dict(L=L, beta=beta, plaq=plaq, dimH=dimH, K=Ks, ric_min=lo[0], ric_max=hi[-1],
                delta_min=ev[0], delta_next=ev[1])


EPS = np.zeros((3, 3, 3))
for i, j, k_ in [(0,1,2),(1,2,0),(2,0,1)]:
    EPS[i, j, k_] = 1.0; EPS[i, k_, j] = -1.0

if __name__ == "__main__":
    res = []
    for L, beta in [(2, 2.0), (2, 2.4), (3, 2.0), (3, 2.4)]:
        res.append(analyse(L, beta, nsweep=300, nplanes=3000,
                           label="SU(2) Wilson"))
    print("\n================ SUMMARY ================")
    print(f"{'L':>3} {'beta':>6} {'<plaq>':>8} {'dimH':>6} {'K_med':>11} {'K_min':>11} "
          f"{'Ric_min':>12} {'Ric_max':>12}")
    for r in res:
        print(f"{r['L']:>3} {r['beta']:>6.2f} {r['plaq']:>8.4f} {r['dimH']:>6} "
              f"{np.median(r['K']):>11.3e} {r['K'].min():>11.3e} "
              f"{r['ric_min']:>12.4e} {r['ric_max']:>12.4e}")
