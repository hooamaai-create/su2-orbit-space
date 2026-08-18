"""
KAGGLE GPU RUN 4 — decide the sign: does horizon proximity strengthen or weaken IR gluons?

Run 3 found r(lambda_min, D(p_min)) = -0.32 at 2.2 sigma, n = 48 — SIGN OPPOSITE the naive
Gribov-Zwanziger intuition (near-horizon configs showed STRONGER low-momentum gluons).
Unresolved. This run cuts the noise three ways instead of only adding statistics:
  1. n = 144 configs (3-sigma sensitivity |r| >= 0.24; if the true r is 0.32, expect ~4 sigma)
  2. D(p_min) averaged over ALL FOUR lattice axes (run 3 used one axis)
  3. Spearman rank correlation alongside Pearson (lambda_min is long-tailed), and
     PARTIAL correlation controlling the per-config plaquette — if r survives with the UV
     confound removed, the link is genuinely infrared.

PRE-REGISTERED (on the axis-averaged D at the lowest momentum, Fisher on the PARTIAL r):
 (G1sign) partial r < 0, |sig| >= 3           -> anomalous sign CONFIRMED: near-horizon
          configs have stronger IR gluons. A genuine, unpublished per-config result.
 (G2sign) partial r > 0, |sig| >= 3           -> naive GZ restored; run 3's sign was noise.
 (G3sign) |partial r| < 0.15                  -> no config-local link (the E2 outcome).
 (G4sign) else                                -> quote the bound |r| <= value and STOP:
          this question is then declared closed-by-bound, not chased further.
"""
import json
import time
import numpy as np
import torch
from scipy.sparse.linalg import LinearOperator, eigsh

FAST = False                      # True = smoke test only
DEV = "cuda" if torch.cuda.is_available() else "cpu"
DT = torch.float64
torch.manual_seed(20260813)
print(f"device = {DEV}  ({torch.cuda.get_device_name(0) if DEV=='cuda' else 'no GPU'})")

# a*sqrt(sigma) from chi(3,3) at L=8, corrected sampler (2026-08-12); 2.0/3.2 extrapolated
A_SQS = {2.0: 0.662, 2.4: 0.412, 2.8: 0.260, 3.2: 0.193}


# ---------------- quaternion SU(2), batched (B, V, 4) ----------------
def qmul(p, q):
    w1, x1, y1, z1 = p.unbind(-1)
    w2, x2, y2, z2 = q.unbind(-1)
    return torch.stack([w1*w2 - x1*x2 - y1*y2 - z1*z2,
                        w1*x2 + x1*w2 + y1*z2 - z1*y2,
                        w1*y2 - x1*z2 + y1*w2 + z1*x2,
                        w1*z2 + x1*y2 - y1*x2 + z1*w2], dim=-1)


def qconj(p):
    out = p.clone()
    out[..., 1:] = -out[..., 1:]
    return out


def qexp(v):
    """exp of pure quaternion, v shaped (..., 3)."""
    n = v.norm(dim=-1, keepdim=True)
    w = torch.cos(n)
    s = torch.where(n > 1e-14, torch.sin(n) / n.clamp(min=1e-14), torch.ones_like(n))
    return torch.cat([w, s * v], dim=-1)


class Lat:
    """B independent SU(2) lattices of shape (L,)*4, links (B, D, V, 4)."""

    def __init__(self, L, B):
        self.L, self.B, self.D = L, B, 4
        self.shape = (L,) * 4
        self.V = L ** 4
        idx = np.arange(self.V).reshape(self.shape)
        self.fwd = [torch.tensor(np.roll(idx, -1, mu).ravel(), device=DEV) for mu in range(4)]
        self.bwd = [torch.tensor(np.roll(idx, +1, mu).ravel(), device=DEV) for mu in range(4)]
        self.U = torch.zeros(B, 4, self.V, 4, dtype=DT, device=DEV)
        self.U[..., 0] = 1.0
        # colour classes per direction (parity if orthogonal extents even, else mod-k)
        coords = np.indices(self.shape)
        self.classes = []
        for mu in range(4):
            others = [nu for nu in range(4) if nu != mu]
            if all(self.shape[nu] % 2 == 0 for nu in others):
                key = sum(coords[nu] for nu in others) % 2
            else:
                key = np.zeros(self.shape, dtype=np.int64)
                for nu in others:
                    k = 2 if self.shape[nu] % 2 == 0 else 3
                    key = key * k + (coords[nu] % k)
            self.classes.append([torch.tensor(idx[key == c].ravel(), device=DEV)
                                 for c in np.unique(key)])

    def staple(self, mu):
        S = torch.zeros(self.B, self.V, 4, dtype=DT, device=DEV)
        U = self.U
        for nu in range(4):
            if nu == mu:
                continue
            S += qmul(qmul(U[:, nu][:, self.fwd[mu]], qconj(U[:, mu][:, self.fwd[nu]])),
                      qconj(U[:, nu]))
            S += qmul(qmul(qconj(U[:, nu][:, self.bwd[nu]][:, self.fwd[mu]]),
                           qconj(U[:, mu][:, self.bwd[nu]])), U[:, nu][:, self.bwd[nu]])
        return S

    def sweep(self, beta, eps=0.35, hits=2):
        """Metropolis, staples recomputed per colour class (detailed balance holds)."""
        for mu in range(4):
            for _ in range(hits):
                for cls in self.classes[mu]:
                    S = self.staple(mu)[:, cls]
                    Uo = self.U[:, mu][:, cls]
                    v = torch.randn(self.B, len(cls), 3, dtype=DT, device=DEV) * eps
                    Un = qmul(qexp(v), Uo)
                    dS = -beta * (qmul(Un, S)[..., 0] - qmul(Uo, S)[..., 0])
                    acc = (dS <= 0) | (torch.rand(self.B, len(cls), dtype=DT, device=DEV)
                                       < torch.exp(-dS.clamp(min=0, max=50)))
                    upd = torch.where(acc.unsqueeze(-1), Un, Uo)
                    self.U[:, mu].index_copy_(1, cls, upd)

    def plaq(self):
        tot = torch.zeros(self.B, dtype=DT, device=DEV)
        for mu in range(4):
            for nu in range(mu + 1, 4):
                P = qmul(qmul(self.U[:, mu], self.U[:, nu][:, self.fwd[mu]]),
                         qmul(qconj(self.U[:, mu][:, self.fwd[nu]]), qconj(self.U[:, nu])))
                tot += P[..., 0].mean(dim=1)
        return tot / 6.0

    # -------- Landau gauge + Faddeev-Popov operator --------
    def grad(self, U=None):
        """-dF/da of the Landau functional; (B, V, 3)."""
        U = self.U if U is None else U
        g = torch.zeros(self.B, self.V, 3, dtype=DT, device=DEV)
        for mu in range(4):
            A = U[:, mu][..., 1:]
            g += A - A[:, self.bwd[mu]]
        return g

    def transform(self, v):
        """links gauge-transformed by exp(v), v (B, V, 3) — returns new link tensor."""
        r = qexp(v)
        return torch.stack([qmul(qmul(r, self.U[:, mu]),
                                 qconj(r[:, self.fwd[mu]])) for mu in range(4)], dim=1)

    def landau_fix(self, tol=1e-9, itmax=30000, alpha=0.08):
        for _ in range(itmax):
            g = self.grad()
            theta = (g ** 2).sum(dim=(1, 2)) / self.V
            if theta.max() < tol:
                return theta
            self.U = self.transform(-alpha * g)
        return theta

    def lam_min_M(self, b, h=1e-4, k=6):
        """Smallest FP eigenvalue above the 3 zero modes, config b. scipy eigsh, GPU matvec.
        M = finite-difference Hessian of the gauge functional (symmetric by construction)."""
        Ub = self.U[b:b+1]
        sub = Lat.__new__(Lat)
        sub.__dict__.update(self.__dict__)
        sub.B, sub.U = 1, Ub

        def mv(x):
            v = torch.tensor(x, dtype=DT, device=DEV).reshape(1, self.V, 3)
            gp = sub.grad(sub.transform(+h * v))
            gm = sub.grad(sub.transform(-h * v))
            return ((gp - gm) / (2 * h)).reshape(-1).cpu().numpy()

        A = LinearOperator((3 * self.V, 3 * self.V), matvec=mv, dtype=np.float64)
        ev = np.sort(eigsh(A, k=k, which='SA', tol=1e-7, maxiter=50000,
                           return_eigenvectors=False))
        return ev



def fisher_sig(r, n):
    z = 0.5 * np.log((1 + r) / (1 - r))
    return z * np.sqrt(max(n - 3, 1))


def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))


def spearman(a, b):
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    return pearson(ra, rb)


def partial_r(rxy, rxz, ryz):
    return (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))


def gluon_D_allaxes(lat, k):
    """Per-config D at momentum 2 pi k / L, AVERAGED over the four axis orientations."""
    B, V, L = lat.B, lat.V, lat.L
    A = lat.U[..., 1:].reshape(B, 4, L, L, L, L, 3)
    At = torch.fft.fftn(A, dim=(2, 3, 4, 5))
    acc = torch.zeros(B, dtype=DT, device=DEV)
    for ax in range(4):
        idx = [0, 0, 0, 0]
        idx[ax] = k
        v = At[:, :, idx[0], idx[1], idx[2], idx[3], :]
        acc += (v.real ** 2 + v.imag ** 2).sum(dim=(1, 2))
    return (acc / 4.0).cpu().numpy() / V


if __name__ == "__main__":
    t0 = time.time()
    results = {"gates": {}, "run4": {}}

    print("\n=== GATES (G1 + G3; full set passed three times on this hardware) ===",
          flush=True)
    lat = Lat(8, 4 if not FAST else 2)
    for _ in range(300 if not FAST else 30):
        lat.sweep(2.4)
    p = lat.plaq().mean().item()
    g1 = abs(p - 0.6285) < 0.006 if not FAST else True
    print(f"G1 <plaq> = {p:.4f}  {'PASS' if g1 else 'FAIL'}")
    latf = Lat(6, 1)
    evf = latf.lam_min_M(0)
    g3 = abs(evf[3] - 4 * np.sin(np.pi / 6) ** 2) < 1e-4
    print(f"G3 free = {evf[3]:.6f}  {'PASS' if g3 else 'FAIL'}")
    results["gates"] = {"plaq": p, "free": float(evf[3]), "all_pass": bool(g1 and g3)}
    if not (g1 and g3):
        raise SystemExit("GATE FAILED - no physics number from this run is valid.")

    L, N = (10, 144) if not FAST else (6, 6)
    print(f"\n=== ensemble: beta = 2.4, L = {L}, n = {N} ===", flush=True)
    lat = Lat(L, N)
    for i in range(500 if not FAST else 40):
        lat.sweep(2.4)
        if (i + 1) % 100 == 0:
            print(f"  sweep {i+1}", flush=True)
    plq = lat.plaq().cpu().numpy()                       # per-config UV control variable
    print(f"  <plaq> = {plq.mean():.4f}", flush=True)

    lat.landau_fix()
    D1 = gluon_D_allaxes(lat, 1)                          # lowest nonzero momentum
    D2 = gluon_D_allaxes(lat, 2)
    lams = np.full(N, np.nan)
    for b in range(N):
        try:
            lams[b] = float(lat.lam_min_M(b)[3])
        except Exception as e:
            print(f"  cfg {b} discarded ({type(e).__name__})", flush=True)
        if (b + 1) % 24 == 0:
            print(f"  lambda_min: {b+1}/{N}  ({time.time()-t0:.0f}s)", flush=True)

    ok = np.isfinite(lams)
    n = int(ok.sum())
    lam = lams[ok]; d1 = D1[ok]; d2 = D2[ok]; pl = plq[ok]
    print(f"\n=== RUN 4 VERDICT INPUTS (n = {n}) ===")
    for tag, d in (("D(p1)", d1), ("D(p2)", d2)):
        r_p = pearson(lam, d)
        r_s = spearman(lam, d)
        r_lp = pearson(lam, pl)
        r_dp = pearson(d, pl)
        r_part = partial_r(r_p, r_lp, r_dp)
        sig_raw = fisher_sig(r_p, n)
        sig_part = fisher_sig(r_part, n)
        results["run4"][tag] = {"pearson": r_p, "spearman": r_s,
                                "r_lam_plaq": r_lp, "r_D_plaq": r_dp,
                                "partial": r_part,
                                "sigma_raw": float(sig_raw),
                                "sigma_partial": float(sig_part)}
        print(f"  {tag}:  pearson {r_p:+.4f} ({sig_raw:+.2f} sig)   spearman {r_s:+.4f}")
        print(f"          confounds: r(lam,plaq) {r_lp:+.3f}  r(D,plaq) {r_dp:+.3f}")
        print(f"          PARTIAL r (plaq removed) = {r_part:+.4f}  ({sig_part:+.2f} sig)")

    print("\n  rules on partial r at D(p1):")
    print("  G1sign  r<0, |sig|>=3  -> anomalous sign CONFIRMED")
    print("  G2sign  r>0, |sig|>=3  -> naive GZ restored; run-3 sign was noise")
    print("  G3sign  |r|<0.15       -> no config-local link")
    print("  G4sign  else           -> quote the bound and STOP (closed-by-bound)")

    with open("results4.json", "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nwrote results4.json   total {time.time() - t0:.0f}s")
