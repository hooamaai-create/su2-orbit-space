"""
KAGGLE GPU RUN 6 — POWER RUN: the divisor, the third beta, and the coherence second leg.

PART K  THE EMPIRICAL DIVISOR, FINALLY WITH POWER.  chi(3,3) and chi(4,4) from Wilson
        loops at L = 12 and 16, beta = 2.2 / 2.5 / 2.8, n = 20/12.  Every dimension in
        the programme divides by d log a/d beta; the current value (-1.17) comes from
        chi(3,3) at L=8 with unproven loop-size convergence.
         (K1) chi(3,3) and chi(4,4) slopes agree < 2 sigma at L = 16
              -> loop-size CONVERGED; one final divisor number with error bars.
         (K2) still apart > 3 sigma -> divisor remains a bound; permanent caveat stated.
PART L  ATTRACTION EXPONENT AT beta = 2.6 (clean solver), L = 6..12.
        With run 5's clean 2.2 (-2.31(21)) and 2.4 (-3.08(25)):
         (L1) p(2.6) steeper than p(2.4) -> eps(beta) RISES toward continuum at 3 points;
              the coherence mechanism's first leg is a trend, not a two-point accident.
         (L2) p(2.6) ~ p(2.2) -> the 2.4 point was the outlier; re-examine.
PART M  COHERENCE SECOND LEG: partial r(lambda, D(p1)) at beta = 2.2 and 2.8, L = 10,
        n = 48 each.  Known: r(2.4) = -0.50(n=72)/-0.44(n=144).
         (M1) |r(2.2)| < |r(2.4)| < |r(2.8)| (ordering holds within errors)
              -> coherence confirmed by a second, independent observable chain.
         (M2) ordering violated -> the crossover story needs revision; report plainly.

GATES: G1 plaquette anchor; G3 deflated free-field exact. FAST=True smoke ~10 min;
full run ~5-7 h. Writes results6.json.
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
        """Smallest FP eigenvalue above the 3 zero modes, config b.

        BUG #13 FIX (audit 2026-08-16): plain eigsh(k=6,'SA',random v0) silently drops a
        member of the exact zero-triplet on ~27% of solves when lambda < 0.15 (all L>=10
        ensembles live there), making ev[3] the SECOND physical eigenvalue. Audit-validated
        fix: DEFLATE the 3 constant colour modes exactly — solve P M P + 10 (1-P) whose
        low spectrum is purely physical — and return in the legacy layout
        [0,0,0, lam1, lam2, lam3] so all callers using ev[3] keep working."""
        Ub = self.U[b:b+1]
        sub = Lat.__new__(Lat)
        sub.__dict__.update(self.__dict__)
        sub.B, sub.U = 1, Ub
        V = self.V

        def proj(x):                       # remove the 3 constant colour modes exactly
            w = x.reshape(V, 3)
            return (w - w.mean(axis=0, keepdims=True)).reshape(-1)

        def mv(x):
            xp = proj(x)
            v = torch.tensor(xp, dtype=DT, device=DEV).reshape(1, V, 3)
            gp = sub.grad(sub.transform(+h * v))
            gm = sub.grad(sub.transform(-h * v))
            z = ((gp - gm) / (2 * h)).reshape(-1).cpu().numpy()
            return proj(z) + 10.0 * (x - xp)

        A = LinearOperator((3 * V, 3 * V), matvec=mv, dtype=np.float64)
        lam = np.sort(eigsh(A, k=3, which='SA', tol=1e-7, maxiter=50000,
                            return_eigenvectors=False))
        return np.concatenate([[0.0, 0.0, 0.0], lam])



def wilson_loop(lat, mu, nu, R, T):
    """Planar R x T Wilson loop in the (mu,nu) plane, batched, averaged over sites."""
    B, V = lat.B, lat.V
    W = torch.zeros(B, V, 4, dtype=DT, device=DEV)
    W[..., 0] = 1.0
    pos = torch.arange(V, device=DEV)
    for _ in range(R):
        W = qmul(W, lat.U[:, mu][:, pos])
        pos = lat.fwd[mu][pos]
    for _ in range(T):
        W = qmul(W, lat.U[:, nu][:, pos])
        pos = lat.fwd[nu][pos]
    for _ in range(R):
        pos = lat.bwd[mu][pos]
        W = qmul(W, qconj(lat.U[:, mu][:, pos]))
    for _ in range(T):
        pos = lat.bwd[nu][pos]
        W = qmul(W, qconj(lat.U[:, nu][:, pos]))
    return W[..., 0].mean(dim=1)                      # (B,) per-config loop average


def loops_upto(lat, rmax):
    """W[R,T] averaged over all planes and configs; returns (rmax+1, rmax+1) numpy."""
    out = np.zeros((rmax + 1, rmax + 1))
    npl = 0
    for mu in range(4):
        for nu in range(mu + 1, 4):
            npl += 1
            for R in range(1, rmax + 1):
                for T in range(1, rmax + 1):
                    out[R, T] += wilson_loop(lat, mu, nu, R, T).mean().item()
    return out / npl


def creutz(W, R, T):
    num = W[R, T] * W[R - 1, T - 1]
    den = W[R - 1, T] * W[R, T - 1]
    if num <= 0 or den <= 0:
        return np.nan
    return -np.log(num / den)


def wfit(betas, vals):
    b = np.array(betas, float); y = np.log(np.abs(np.array(vals, float)))
    return np.polyfit(b, y, 1)[0]


def fit_p(Ls, ms, es):
    x = np.log(Ls); y = np.log(ms); w = np.array(es) / np.array(ms)
    W = 1 / w ** 2
    den = W.sum() * (W * x * x).sum() - (W * x).sum() ** 2
    p = (W.sum() * (W * x * y).sum() - (W * x).sum() * (W * y).sum()) / den
    return p, np.sqrt(W.sum() / den)


def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))


def partial_r(rxy, rxz, ryz):
    return (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))


def gluon_D_p1(lat):
    B, V, L = lat.B, lat.V, lat.L
    A = lat.U[..., 1:].reshape(B, 4, L, L, L, L, 3)
    At = torch.fft.fftn(A, dim=(2, 3, 4, 5))
    acc = torch.zeros(B, dtype=DT, device=DEV)
    for ax in range(4):
        idx = [0, 0, 0, 0]; idx[ax] = 1
        v = At[:, :, idx[0], idx[1], idx[2], idx[3], :]
        acc += (v.real ** 2 + v.imag ** 2).sum(dim=(1, 2))
    return (acc / 4.0).cpu().numpy() / V


def lam_batch(lat, tag):
    out = np.full(lat.B, np.nan)
    for b in range(lat.B):
        try:
            out[b] = float(lat.lam_min_M(b)[3])
        except Exception as e:
            print(f"    {tag} cfg {b} discarded ({type(e).__name__})", flush=True)
    return out


if __name__ == "__main__":
    t0 = time.time()
    results = {"gates": {}, "partK": {}, "partL": {}, "partM": {}}

    print("\n=== GATES ===", flush=True)
    lat = Lat(8, 4 if not FAST else 2)
    for _ in range(300 if not FAST else 30):
        lat.sweep(2.4)
    p = lat.plaq().mean().item()
    g1 = abs(p - 0.6285) < 0.006 if not FAST else True
    print(f"G1 <plaq> = {p:.4f}  {'PASS' if g1 else 'FAIL'}")
    latf = Lat(6, 1)
    evf = latf.lam_min_M(0)
    g3 = abs(evf[3] - 1.0) < 1e-4
    print(f"G3 deflated free = {evf[3]:.6f}  {'PASS' if g3 else 'FAIL'}")
    results["gates"] = {"plaq": p, "free": float(evf[3]), "all_pass": bool(g1 and g3)}
    if not (g1 and g3):
        raise SystemExit("GATE FAILED - no physics number from this run is valid.")

    # ---------------- PART K: the divisor with power ----------------
    KBETAS = [2.2, 2.5, 2.8]
    KPLAN = [(12, 20), (16, 12)] if not FAST else [(6, 2)]
    for L, n in KPLAN:
        c33, c44 = [], []
        for beta in KBETAS:
            lat = Lat(L, n)
            for _ in range(400 if not FAST else 30):
                lat.sweep(beta)
            acc = np.zeros((5, 5)); m = 0
            for _ in range(8 if not FAST else 2):
                for _ in range(3):
                    lat.sweep(beta)
                acc += loops_upto(lat, 4); m += 1
            W = acc / m
            c33.append(creutz(W, 3, 3)); c44.append(creutz(W, 4, 4))
            results["partK"][f"L{L},b{beta}"] = {"chi33": c33[-1], "chi44": c44[-1],
                                                "plaq": lat.plaq().mean().item()}
            print(f"  K L={L} beta={beta}: chi33={c33[-1]:.5f} chi44={c44[-1]:.5f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
        s33 = wfit(KBETAS, c33); s44 = wfit(KBETAS, c44)
        results["partK"][f"L{L}_slopes"] = {"s33": s33, "s44": s44,
                                            "div33": s33 / 2, "div44": s44 / 2}
        print(f"  K L={L}: d log chi33/dbeta = {s33:+.4f} -> div {s33/2:+.4f}; "
              f"chi44 -> div {s44/2:+.4f}", flush=True)
    print("  K1 chi33/chi44 divisors agree at L=16 -> CONVERGED; final number")
    print("  K2 apart -> divisor stays a bound")

    # ---------------- PART L: exponent at beta = 2.6 ----------------
    BETA_L = 2.6
    PLAN = [(6, 20), (8, 16), (10, 12), (12, 8)] if not FAST else [(4, 2), (6, 2)]
    print(f"\n=== PART L: clean exponent, beta = {BETA_L} ===", flush=True)
    Ls, ms, es = [], [], []
    for L, n in PLAN:
        lat = Lat(L, n)
        for _ in range(500 if not FAST else 40):
            lat.sweep(BETA_L)
        th = lat.landau_fix()
        assert th.max().item() < 1e-8
        lam = lam_batch(lat, f"L{L}")
        lam = lam[np.isfinite(lam)]
        Ls.append(L); ms.append(float(lam.mean()))
        es.append(float(lam.std(ddof=1) / np.sqrt(len(lam))))
        results["partL"][L] = {"lam": ms[-1], "err": es[-1], "n": len(lam)}
        print(f"  L={L} lam={ms[-1]:.5f}+-{es[-1]:.5f}", flush=True)
    pL, eL = fit_p(Ls, ms, es)
    pfree = np.polyfit(np.log(Ls), np.log([4 * np.sin(np.pi / L) ** 2 for L in Ls]), 1)[0]
    results["partL"]["fit"] = {"p": pL, "err": eL, "p_free": pfree}
    print(f"  p(2.6) = {pL:+.4f} +- {eL:.4f}  free {pfree:+.4f}"
          f"  excess {(pL-pfree):+.3f} ({abs(pL-pfree)/eL:.1f} sigma)")
    print("  known clean: p(2.2) = -2.309(211), p(2.4) = -3.081(246)")
    print("  L1 p(2.6) steeper than p(2.4) -> eps rises at 3 points   L2 ~p(2.2) -> outlier")

    # ---------------- PART M: coherence second leg ----------------
    print("\n=== PART M: r(lambda, D) at beta = 2.2 and 2.8 ===", flush=True)
    for beta in ([2.2, 2.8] if not FAST else [2.4]):
        L, N = (10, 48) if not FAST else (6, 4)
        lat = Lat(L, N)
        for _ in range(500 if not FAST else 40):
            lat.sweep(beta)
        plq = lat.plaq().cpu().numpy()
        th = lat.landau_fix()
        assert th.max().item() < 1e-8
        D1 = gluon_D_p1(lat)
        lam = lam_batch(lat, f"M{beta}")
        ok = np.isfinite(lam)
        l, d, pl = lam[ok], D1[ok], plq[ok]
        r = pearson(l, d)
        rp = partial_r(r, pearson(l, pl), pearson(d, pl))
        results["partM"][beta] = {"r": r, "partial": rp, "n": int(ok.sum())}
        print(f"  beta={beta}: partial r = {rp:+.4f} (n={ok.sum()})", flush=True)
    print("  known: r(2.4) = -0.50 (n=72) / -0.44 (n=144)")
    print("  M1 |r(2.2)| < |r(2.4)| < |r(2.8)| -> coherence confirmed by second leg")
    print("  M2 ordering violated -> crossover story needs revision")

    with open("results6.json", "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nwrote results6.json   total {time.time() - t0:.0f}s")
