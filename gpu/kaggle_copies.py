"""
KAGGLE GPU RUN 2 — copy systematic + second coupling for the horizon attraction.

Previous run: lambda_min(M) ~ L^-2.86(20) vs free L^-1.90 at beta = 2.4 — 4.8 sigma
dynamical Gribov-horizon attraction. Two caveats were logged with that result; this run
closes both:
  PART C  single-copy gauge fixing (copy effects grow with volume): measure first-copy vs
          best-of-5-copies lambda_min at L = 6 and 10; compare the two exponents.
  PART D  one coupling only: repeat the L = 6..12 scan at beta = 2.2.
Gates identical to run 1 (they passed on this hardware). FAST = True first: ~10 min.
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



def fit_p(Ls, ms, es):
    x = np.log(Ls); y = np.log(ms); w = np.array(es) / np.array(ms)
    W = 1 / w ** 2
    den = W.sum() * (W * x * x).sum() - (W * x).sum() ** 2
    p = (W.sum() * (W * x * y).sum() - (W * x).sum() * (W * y).sum()) / den
    return p, np.sqrt(W.sum() / den)


def rand_gauge(lat):
    """Random gauge transform (large rotations) of the current links."""
    v = torch.randn(lat.B, lat.V, 3, dtype=DT, device=DEV) * 2.0
    lat.U = lat.transform(v)


def lam_of_batch(lat, tag):
    out = np.full(lat.B, np.nan)
    for b in range(lat.B):
        try:
            out[b] = float(lat.lam_min_M(b)[3])
        except Exception as e:
            print(f"    {tag} cfg {b} discarded ({type(e).__name__})", flush=True)
    return out


if __name__ == "__main__":
    t0 = time.time()
    results = {"gates": {}, "partC": {}, "partD": {}}

    # ------------- GATES (same three, validated on the previous run) -------------
    print("\n=== GATES ===", flush=True)
    lat = Lat(8, 4 if not FAST else 2)
    for _ in range(300 if not FAST else 30):
        lat.sweep(2.4)
    p = lat.plaq().mean().item()
    g1 = abs(p - 0.6285) < 0.006 if not FAST else True
    print(f"G1 <plaq>(2.4, L=8) = {p:.4f}  target 0.6285+-0.006  {'PASS' if g1 else 'FAIL'}")
    lat.landau_fix()
    x = torch.randn(1, lat.V, 3, dtype=DT, device=DEV)
    y = torch.randn(1, lat.V, 3, dtype=DT, device=DEV)
    h = 1e-4
    def Mv(v):
        return (lat.grad(lat.transform(+h * v)) - lat.grad(lat.transform(-h * v))) / (2 * h)
    sym = abs((x * Mv(y)).sum().item() - (y * Mv(x)).sum().item()) / abs((x * Mv(y)).sum().item())
    ev = lat.lam_min_M(0)
    near0 = bool(np.all(np.abs(ev[:3]) < 0.01 * abs(ev[3])))
    g2 = sym < 1e-4 and ev[3] > 0 and near0
    print(f"G2 sym={sym:.2e}  ev[0:4]={np.array2string(ev[:4], precision=6)}  "
          f"{'PASS' if g2 else 'FAIL'}")
    latf = Lat(6, 1)
    evf = latf.lam_min_M(0)
    g3 = abs(evf[3] - 4 * np.sin(np.pi / 6) ** 2) < 1e-4
    print(f"G3 free lam[3] = {evf[3]:.6f}  {'PASS' if g3 else 'FAIL'}")
    results["gates"] = {"plaq": p, "sym": sym, "all_pass": bool(g1 and g2 and g3)}
    if not (g1 and g2 and g3):
        raise SystemExit("GATE FAILED - no physics number from this run is valid.")

    # ------------- PART C: Gribov-copy systematic on the exponent -------------
    # The 4.8-sigma attraction used SINGLE-COPY steepest-descent gauge fixing. Copy effects
    # grow with volume and could inflate the exponent (Sternbeck et al.: lambda_1 differs
    # between copies). NC extra random-start copies per config at L = 6 and 10; compare
    # first-copy vs best-F-copy lambda_min and the 2-point exponents.
    #  (C1) exponents fc vs bc agree < 2 sigma -> attraction is NOT a copy artifact.
    #  (C2) differ > 3 sigma                   -> copy artifact drives it; DOWNGRADE.
    #  (C3) else -> report; needs more copies/volumes.
    NC = 4 if not FAST else 1
    NCFG_C = 12 if not FAST else 2
    print("\n=== PART C: copy systematic, beta = 2.4, L = 6 & 10 ===", flush=True)
    exps = {}
    for L in ([6, 10] if not FAST else [4, 6]):
        lat = Lat(L, NCFG_C)
        for _ in range(500 if not FAST else 40):
            lat.sweep(2.4)
        U0 = lat.U.clone()
        lat.landau_fix()
        lam_fc = lam_of_batch(lat, f"L{L} fc")
        F_best = lat.U[..., 0].sum(dim=(1, 2)).cpu().numpy()
        lam_bc = lam_fc.copy()
        for c in range(NC):
            lat.U = U0.clone()
            rand_gauge(lat)
            lat.landau_fix()
            F_c = lat.U[..., 0].sum(dim=(1, 2)).cpu().numpy()
            lam_c = lam_of_batch(lat, f"L{L} copy{c}")
            better = F_c > F_best
            lam_bc = np.where(better, lam_c, lam_bc)
            F_best = np.maximum(F_c, F_best)
        nf = np.isfinite(lam_fc).sum(); nb_ = np.isfinite(lam_bc).sum()
        m_fc = np.nanmean(lam_fc); e_fc = np.nanstd(lam_fc) / np.sqrt(nf)
        m_bc = np.nanmean(lam_bc); e_bc = np.nanstd(lam_bc) / np.sqrt(nb_)
        chg = np.nanmean((lam_bc - lam_fc) / lam_fc)
        exps[L] = (m_fc, e_fc, m_bc, e_bc)
        results["partC"][L] = {"fc": [m_fc, e_fc], "bc": [m_bc, e_bc],
                               "rel_change": float(chg)}
        print(f"  L={L}: fc lam = {m_fc:.4f}+-{e_fc:.4f}   bc lam = {m_bc:.4f}+-{e_bc:.4f}"
              f"   mean rel change = {100*chg:+.1f}%", flush=True)
    Ls = sorted(exps)
    if len(Ls) == 2:
        (a1, ea1, b1, eb1), (a2, ea2, b2, eb2) = exps[Ls[0]], exps[Ls[1]]
        dx = np.log(Ls[1] / Ls[0])
        p_fc = np.log(a2 / a1) / dx; ep_fc = np.hypot(ea1 / a1, ea2 / a2) / dx
        p_bc = np.log(b2 / b1) / dx; ep_bc = np.hypot(eb1 / b1, eb2 / b2) / dx
        sig = abs(p_fc - p_bc) / np.hypot(ep_fc, ep_bc)
        results["partC"]["exponents"] = {"fc": [p_fc, ep_fc], "bc": [p_bc, ep_bc],
                                         "sigma": float(sig)}
        print(f"  2-point exponent: fc {p_fc:+.3f}+-{ep_fc:.3f}   "
              f"bc {p_bc:+.3f}+-{ep_bc:.3f}   differ by {sig:.1f} sigma")
        print("  C1 < 2 sigma -> not a copy artifact   C2 > 3 sigma -> DOWNGRADE")

    # ------------- PART D: the exponent at a second coupling -------------
    #  (D1) p(2.2) below free by > 3 sigma too -> attraction robust in beta.
    #  (D2) p(2.2) consistent with free       -> the 2.4 result is beta-fragile; flag.
    BETA_D = 2.2
    print(f"\n=== PART D: exponent at beta = {BETA_D} ===", flush=True)
    PLAN = [(6, 20), (8, 16), (10, 12), (12, 8)] if not FAST else [(4, 2), (6, 2)]
    Ls, ms, es = [], [], []
    for L, n in PLAN:
        lat = Lat(L, n)
        for _ in range(500 if not FAST else 40):
            lat.sweep(BETA_D)
        pq = lat.plaq().mean().item()
        lat.landau_fix()
        lam = lam_of_batch(lat, f"D L{L}")
        lam = lam[np.isfinite(lam)]
        Ls.append(L); ms.append(float(lam.mean()))
        es.append(float(lam.std(ddof=1) / np.sqrt(len(lam))))
        results["partD"][L] = {"lam": ms[-1], "err": es[-1], "n": len(lam), "plaq": pq}
        print(f"  L={L} n={len(lam)} <plaq>={pq:.4f} lam={ms[-1]:.4f}+-{es[-1]:.4f}",
              flush=True)
    pD, epD = fit_p(Ls, ms, es)
    pfree = np.polyfit(np.log(Ls), np.log([4 * np.sin(np.pi / L) ** 2 for L in Ls]), 1)[0]
    results["partD"]["fit"] = {"p": pD, "err": epD, "p_free": pfree}
    print(f"\n  p({BETA_D}) = {pD:+.4f} +- {epD:.4f}   free = {pfree:+.4f}"
          f"   excess {(pD - pfree):+.3f} ({abs(pD - pfree)/epD:.1f} sigma)")
    print(f"  previous run: p(2.4) = -2.8595 +- 0.2017 (excess -0.96, 4.8 sigma)")
    print("  D1 below free > 3 sigma -> robust in beta   D2 within errors -> beta-fragile")

    with open("results2.json", "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nwrote results2.json   total {time.time() - t0:.0f}s")
