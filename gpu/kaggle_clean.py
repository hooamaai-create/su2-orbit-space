"""
KAGGLE GPU RUN 5 — THE CLEAN FLOORS. Deflated solver throughout (bug #13 fixed).

PART H  Attraction exponent, clean: L = 6,8,10,12 at beta = 2.4 AND 2.2.
        Old (contaminated, conservative floors): -2.86(20) and -2.68(23) vs free -1.90.
         (H1) clean exponents below free by >3 sigma at both betas -> the excess is real;
              quote epsilon(SU(2)) vs the published SU(3) 0.16 as an open discrepancy.
         (H2) clean exponents collapse toward free -> the excess was solver contamination;
              RETRACT the "sharpest in-window" framing. Either way: decided.
PART I  The config-level law under best-copy gauge fixing (referee objection, 1211.3057):
        n = 72 at L = 10, beta = 2.4; per config measure (lambda, D(p1)) on the first copy
        AND on the best-of-4-random-restart copy; partial r (plaquette removed) for each.
         (I1) r_bc consistent with r_fc (< 2 sigma) -> law survives; claimable with the
              mandatory citation.  (I2) |r_bc| < 0.2 -> the law was copy landing; RETRACT.
         (I3) else -> unresolved, report.
PART J  Guarded Part-A repair: s at the matched-ell triple (2.0,L3) (2.4,L5) (2.8,L8),
        n = 24, clean solver.  (J1) collapse still fails >3 sigma -> A2 verdict upgraded
        back to quantitative.  (J2) collapse now HOLDS -> the 11.8-sigma kill was solver
        contamination; RETRACT A2 and reopen the universality question.

GATES: G1 plaquette (L=8, 0.6285+-0.006+3se), G3 free-field lam = 1.000000 via the
DEFLATED path (validates the fix end-to-end). Abort on failure.
FAST=True smoke first (~10 min). Full run ~4-6 h. Writes results5.json.
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


def fisher_sig(r, n):
    return 0.5 * np.log((1 + r) / (1 - r)) * np.sqrt(max(n - 4, 1))


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


def rand_gauge(lat):
    v = torch.randn(lat.B, lat.V, 3, dtype=DT, device=DEV) * 2.0
    lat.U = lat.transform(v)


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
    results = {"gates": {}, "partH": {}, "partI": {}, "partJ": {}}

    print("\n=== GATES (deflated-solver validation) ===", flush=True)
    lat = Lat(8, 4 if not FAST else 2)
    for _ in range(300 if not FAST else 30):
        lat.sweep(2.4)
    p = lat.plaq().mean().item()
    g1 = abs(p - 0.6285) < 0.006 if not FAST else True
    print(f"G1 <plaq> = {p:.4f}  {'PASS' if g1 else 'FAIL'}")
    latf = Lat(6, 1)
    evf = latf.lam_min_M(0)
    g3 = abs(evf[3] - 1.0) < 1e-4
    print(f"G3 free lam (deflated path) = {evf[3]:.6f}  {'PASS' if g3 else 'FAIL'}")
    results["gates"] = {"plaq": p, "free": float(evf[3]), "all_pass": bool(g1 and g3)}
    if not (g1 and g3):
        raise SystemExit("GATE FAILED - no physics number from this run is valid.")

    # ---------------- PART H: clean exponents at two betas ----------------
    PLAN = [(6, 20), (8, 16), (10, 12), (12, 8)] if not FAST else [(4, 2), (6, 2)]
    for BETA in (2.4, 2.2):
        print(f"\n=== PART H: clean exponent, beta = {BETA} ===", flush=True)
        Ls, ms, es = [], [], []
        for L, n in PLAN:
            lat = Lat(L, n)
            for _ in range(500 if not FAST else 40):
                lat.sweep(BETA)
            pq = lat.plaq().mean().item()
            th = lat.landau_fix()
            assert th.max().item() < 1e-8, f"landau nonconvergence at L={L}"
            lam = lam_batch(lat, f"H{BETA} L{L}")
            lam = lam[np.isfinite(lam)]
            Ls.append(L); ms.append(float(lam.mean()))
            es.append(float(lam.std(ddof=1) / np.sqrt(len(lam))))
            results["partH"][f"{BETA},{L}"] = {"lam": ms[-1], "err": es[-1],
                                               "n": len(lam), "plaq": pq}
            print(f"  L={L} n={len(lam)} <plaq>={pq:.4f} lam={ms[-1]:.5f}+-{es[-1]:.5f}",
                  flush=True)
        pH, eH = fit_p(Ls, ms, es)
        pfree = np.polyfit(np.log(Ls), np.log([4 * np.sin(np.pi / L) ** 2 for L in Ls]), 1)[0]
        results["partH"][f"fit{BETA}"] = {"p": pH, "err": eH, "p_free": pfree}
        print(f"  p({BETA}) = {pH:+.4f} +- {eH:.4f}   free {pfree:+.4f}"
              f"   excess {(pH - pfree):+.3f} ({abs(pH - pfree)/eH:.1f} sigma)")
    print("  old contaminated: -2.8595(2017) @2.4, -2.6765(2253) @2.2")
    print("  H1 both below free >3 sigma -> excess real   H2 collapse -> RETRACT framing")

    # ---------------- PART I: the law under best-copy fixing ----------------
    L, N, NC = (10, 72, 4) if not FAST else (6, 6, 1)
    print(f"\n=== PART I: r under fc vs bc, L={L}, n={N}, {NC} restarts ===", flush=True)
    lat = Lat(L, N)
    for _ in range(500 if not FAST else 40):
        lat.sweep(2.4)
    plq = lat.plaq().cpu().numpy()
    U0 = lat.U.clone()
    th = lat.landau_fix()
    assert th.max().item() < 1e-8
    lam_fc = lam_batch(lat, "I fc")
    D_fc = gluon_D_p1(lat)
    F_best = lat.U[..., 0].sum(dim=(1, 2)).cpu().numpy()
    U_best = lat.U.clone()
    for c in range(NC):
        lat.U = U0.clone()
        rand_gauge(lat)
        th = lat.landau_fix()
        F_c = lat.U[..., 0].sum(dim=(1, 2)).cpu().numpy()
        better = torch.tensor(F_c > F_best, device=DEV)
        U_best[better] = lat.U[better]
        F_best = np.maximum(F_c, F_best)
        print(f"  copy {c+1}/{NC} fixed ({time.time()-t0:.0f}s)", flush=True)
    lat.U = U_best
    lam_bc = lam_batch(lat, "I bc")
    D_bc = gluon_D_p1(lat)

    def partial_block(lam, D, tag):
        ok = np.isfinite(lam)
        l, d, pl = lam[ok], D[ok], plq[ok]
        r = pearson(l, d)
        rp = partial_r(r, pearson(l, pl), pearson(d, pl))
        sig = fisher_sig(rp, int(ok.sum()))
        results["partI"][tag] = {"r": r, "partial": rp, "sigma": float(sig),
                                 "n": int(ok.sum())}
        print(f"  {tag}: r={r:+.4f}  partial={rp:+.4f}  ({sig:+.2f} sigma, n={ok.sum()})")
        return rp, int(ok.sum())

    rp_fc, n_fc = partial_block(lam_fc, D_fc, "fc")
    rp_bc, n_bc = partial_block(lam_bc, D_bc, "bc")
    dz = abs(0.5 * np.log((1 + rp_fc) / (1 - rp_fc)) -
             0.5 * np.log((1 + rp_bc) / (1 - rp_bc))) / np.sqrt(1 / (n_fc - 3) + 1 / (n_bc - 3))
    results["partI"]["fc_vs_bc_sigma"] = float(dz)
    print(f"  fc vs bc difference: {dz:.2f} sigma")
    print("  I1 <2 sigma & both nonzero -> law survives   I2 |r_bc|<0.2 -> RETRACT   I3 else")

    # ---------------- PART J: guarded matched-ell triple ----------------
    print(f"\n=== PART J: guarded s at matched ell ~2.0 ===", flush=True)
    NJ = 24 if not FAST else 3
    for beta, L in ([(2.0, 3), (2.4, 5), (2.8, 8)] if not FAST else [(2.4, 4)]):
        lat = Lat(L, NJ)
        for _ in range(500 if not FAST else 40):
            lat.sweep(beta)
        th = lat.landau_fix()
        assert th.max().item() < 1e-8
        lam = lam_batch(lat, f"J {beta},{L}")
        lam = lam[np.isfinite(lam)]
        s = lam / (4 * np.sin(np.pi / L) ** 2)
        results["partJ"][f"{beta},{L}"] = {"s": float(s.mean()),
                                           "err": float(s.std(ddof=1) / np.sqrt(len(s))),
                                           "n": len(s)}
        print(f"  beta={beta} L={L}: s = {s.mean():.4f}+-{s.std(ddof=1)/np.sqrt(len(s)):.4f}",
              flush=True)
    print("  old (contaminated): 0.1239/0.1613/0.4028 -> 11.8 sigma kill")
    print("  J1 still fails >3 sigma -> A2 back to quantitative   J2 collapses -> RETRACT A2")

    with open("results5.json", "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nwrote results5.json   total {time.time() - t0:.0f}s")
