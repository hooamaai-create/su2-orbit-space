"""
MATRIX-FREE lambda_min(Ric) — the rewrite needed for large volumes.

WHY THE OLD CODE COULD NOT SCALE
  geom() built a dense (12V)x(12V) horizontal projector: 75 MB at L=4, 1.9 GB at L=6,
  19 GB at L=8. Everything downstream inherited that wall.

THE RESTRUCTURING
  Ric(a,a) = 3 Sum_i <rho(a,e_i), Delta^-1 rho(a,e_i)>  over an ON basis {e_i} of H.
  Define the linear map  M_a : Omega^1 -> Omega^0,  (M_a e)^c(n) = Sum_mu (a_mu(n) x e_mu(n))^c.
  Since Sum_i e_i e_i^T = P_H,

        Ric(a,a) = 3 tr( P_H M_a^dag Delta^-1 M_a )                    <-- a TRACE

  which is exactly the shape Hutchinson estimation wants:
        tr(P_H B) = tr(P_H B P_H) = E_z[ (P_H z)^dag B (P_H z) ],   z = random +-1
  P_H is APPLIED, never formed:  P_H z = z - D Delta^{-1} D^dag z   (one CG solve).

  With w_z := P_H z FROZEN, Ric becomes a deterministic quadratic form in a:
        Ric ~ (3/N) Sum_z L_{w_z}^dag Delta^-1 L_{w_z},   (L_w a)(n) = Sum_mu a_mu(n) x w_mu(n)
  so Lanczos is valid. Adjoint: (L_w^dag u)_mu(n) = w_mu(n) x u(n), from (axb).c = a.(bxc).

  Cost per Ric matvec: N CG solves on the 3V-dimensional Delta. No 12V x 12V object ever exists.

RESTRICTION TO H
  lambda_min must be taken over H only; the vertical directions would contribute spurious
  zeros below the answer. Handled by working with  Ric + c*(1 - P_H)  with c above the
  expected lambda_min, so vertical directions sit at c and never pollute the bottom.

BIAS WARNING
  finite N makes the estimated operator noisy, and noise biases lambda_min DOWNWARD.
  The convergence of lambda_min in N is therefore tested explicitly against the exact dense
  answer at L=3 before any large-volume run is trusted.
"""
import numpy as np
from step3_lattice import Lattice, adjoint_matrix


class MF:
    def __init__(self, lat, cg_tol=1e-10, cg_max=2000):
        self.lat = lat
        self.D, self.V = lat.D, lat.V
        self.R = np.stack([adjoint_matrix(lat.U[mu]) for mu in range(lat.D)])  # (D,V,3,3)
        self.fwd = lat.fwd
        self.bwd = lat.bwd
        self.cg_tol, self.cg_max = cg_tol, cg_max
        self.nsolve = 0

    # ---- D : Omega^0 -> Omega^1 ,  (Dw)_mu(n) = R_mu(n) w(n+mu) - w(n)
    def Dop(self, w):
        w = w.reshape(self.V, 3)
        out = np.empty((self.D, self.V, 3))
        for mu in range(self.D):
            out[mu] = np.einsum('vij,vj->vi', self.R[mu], w[self.fwd[mu]]) - w
        return out

    # ---- D^dag : Omega^1 -> Omega^0 ,  (D^dag b)(n) = Sum_mu [ R_mu(n-mu)^T b_mu(n-mu) - b_mu(n) ]
    def Ddag(self, b):
        b = b.reshape(self.D, self.V, 3)
        out = np.zeros((self.V, 3))
        for mu in range(self.D):
            bm = b[mu][self.bwd[mu]]
            Rm = self.R[mu][self.bwd[mu]]
            out += np.einsum('vji,vj->vi', Rm, bm) - b[mu]
        return out

    def Delta(self, w):
        return self.Ddag(self.Dop(w)).ravel()

    def cg(self, rhs, x0=None):
        """CG for Delta x = rhs on Omega^0 (3V)."""
        b = rhs.ravel()
        x = np.zeros_like(b) if x0 is None else x0.copy()
        r = b - self.Delta(x)
        p = r.copy()
        rs = r @ r
        nb = np.linalg.norm(b)
        if nb == 0:
            return x
        for _ in range(self.cg_max):
            Ap = self.Delta(p)
            al = rs / (p @ Ap)
            x += al * p
            r -= al * Ap
            rs_new = r @ r
            if np.sqrt(rs_new) / nb < self.cg_tol:
                break
            p = r + (rs_new / rs) * p
            rs = rs_new
        self.nsolve += 1
        return x

    def P_H(self, z):
        """z - D Delta^{-1} D^dag z, applied never formed."""
        z = z.reshape(self.D, self.V, 3)
        return (z - self.Dop(self.cg(self.Ddag(z)))).ravel()

    # ---- the bilinear maps
    @staticmethod
    def L(w, a):
        """(L_w a)(n) = Sum_mu a_mu(n) x w_mu(n)   ->  Omega^0"""
        return np.cross(a, w).sum(axis=0)

    @staticmethod
    def Ladj(w, u):
        """(L_w^dag u)_mu(n) = w_mu(n) x u(n)      ->  Omega^1"""
        return np.cross(w, u[None, :, :])

    def freeze_probes(self, N, seed=0):
        rng = np.random.default_rng(seed)
        self.W = []
        for _ in range(N):
            z = rng.choice([-1.0, 1.0], size=self.D * self.V * 3)
            self.W.append(self.P_H(z).reshape(self.D, self.V, 3))
        self.N = N

    def ric(self, a):
        """(3/N) Sum_z L_wz^dag Delta^-1 L_wz  applied to a."""
        a = a.reshape(self.D, self.V, 3)
        out = np.zeros_like(a)
        for w in self.W:
            u = self.cg(self.L(w, a)).reshape(self.V, 3)
            out += self.Ladj(w, u)
        return (3.0 / self.N) * out.ravel()

    def op(self, a, c=10.0):
        """Ric restricted to H, with vertical directions lifted to c."""
        ah = self.P_H(a)
        return self.ric(ah) + c * (a - ah)


def lanczos_min(apply, n, iters=120, seed=1):
    """Symmetric Lanczos; returns the smallest Ritz value."""
    rng = np.random.default_rng(seed)
    q = rng.normal(size=n)
    q /= np.linalg.norm(q)
    Q = [q]
    alphas, betas = [], []
    beta = 0.0
    qprev = np.zeros(n)
    for j in range(iters):
        v = apply(Q[-1])
        alpha = Q[-1] @ v
        v = v - alpha * Q[-1] - beta * qprev
        # full reorthogonalisation (cheap relative to the CG solves)
        for u in Q:
            v -= (u @ v) * u
        beta = np.linalg.norm(v)
        alphas.append(alpha)
        if beta < 1e-12 or j == iters - 1:
            break
        betas.append(beta)
        qprev = Q[-1]
        Q.append(v / beta)
    T = np.diag(alphas) + np.diag(betas, 1) + np.diag(betas, -1)
    return np.linalg.eigvalsh(T)[0]
