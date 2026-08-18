"""
KAGGLE GPU RUN 11 — IS r = F(f_IR) A LAW? (pre-registered: record/RUN11_SPEC.md, 2026-08-18)

Prompted by hostile external review. Run 10 established that coherence predicts WHETHER
the horizon-gluon correlation is active (2 out-of-sample points). It did NOT establish
a universal law. Three legs, all required:

  Leg 1  dense blinded scan: 12 ensembles spanning f_IR ~ 0.30-0.92, six inside the
         transition band 0.45-0.70. Predicted r for EVERY ensemble is computed from the
         frozen activation curve and written to the results file BEFORE any r is measured.
  Leg 2  independence of the coherence variable: per ensemble report
           f1 = IR fraction of psi_min          (the original, spectral)
           f2 = IR fraction of psi_2, psi_3     (same operator, different vectors)
           f3 = IR fraction of the GLUON FIELD power (non-spectral: built from the links
                alone; answers the circularity attack — f_IR shares an eigenmode with r)
  Leg 3  conditional independence r | f_IR: matched-coherence pairs at different (beta,L)
         exist inside the plan; analyzed after both parts return.

  Also: physical regression slope a = cov(lam,D)/var(lam) per ensemble (units of D/lam),
  so the "locked -0.52" must reappear as systematic behaviour of a, not of a bare
  correlation coefficient.

FALSIFIERS (fixed in RUN11_SPEC.md before any data):
  F1  any ensemble with f1 < 0.45 showing |r| > 0.25 at >2 sigma  -> threshold dead
  F2  f3 fails to gate r (transition absent/displaced by >0.15)   -> spectral artifact
  F3  matched-f_IR pairs differing in r by >3 sigma               -> not single-variable

Two parts (each ~5-7 h on Kaggle): PART A = transition band, PART B = extremes + new
couplings beta = 2.2 / 2.5 / 2.6. Writes results11A.json / results11B.json.
"""
import json
import time
import numpy as np
import torch
from scipy.sparse.linalg import LinearOperator, eigsh

PART = "A"                        # "A" or "B" — set per notebook
FAST = False                      # True = smoke test only
DEV = "cuda" if torch.cuda.is_available() else "cpu"
DT = torch.float64
torch.manual_seed(20260818 + (0 if PART == "A" else 1))
print(f"device = {DEV}  ({torch.cuda.get_device_name(0) if DEV=='cuda' else 'no GPU'})")

# measured f_IR(beta, L) grid (runs 6-10), L in {8,10,12,14,16,18}; 2.3/2.5 interpolated
FIR_GRID = {2.2: [.3806, .3829, .3839, .3008, .3564, .2802],
            2.4: [.687, .622, .560, .464, .444, .409],
            2.6: [.916, .829, .868, .883, .869, .860]}
GRID_L = [8, 10, 12, 14, 16, 18]


def fir_predicted(beta, L):
    """bilinear interpolation on the measured grid (explicit keys, never lo+0.2)."""
    def at(b):
        arr = FIR_GRID[b]
        if L <= GRID_L[0]:
            return arr[0]
        for i in range(1, len(GRID_L)):
            if L <= GRID_L[i]:
                t = (L - GRID_L[i-1]) / (GRID_L[i] - GRID_L[i-1])
                return arr[i-1] + t * (arr[i] - arr[i-1])
        return arr[-1]
    if beta <= 2.2:
        return at(2.2)
    if beta >= 2.6:
        return at(2.6)
    lo = 2.2 if beta < 2.4 else 2.4
    hi = 2.4 if lo == 2.2 else 2.6
    t = (beta - lo) / (hi - lo)
    return at(lo) + t * (at(hi) - at(lo))


def r_predicted(f):
    """the FROZEN activation curve (viz/law.html, fixed 2026-08-17 — before this run)."""
    return -0.52 / (1.0 + np.exp(-(f - 0.52) / 0.045))


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

    def grad(self, U=None):
        U = self.U if U is None else U
        g = torch.zeros(self.B, self.V, 3, dtype=DT, device=DEV)
        for mu in range(4):
            A = U[:, mu][..., 1:]
            g += A - A[:, self.bwd[mu]]
        return g

    def transform(self, v):
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
        """Smallest FP eigenvalue above the 3 zero modes (BUG #13 FIX: exact deflation
        of the constant colour triplet; solve P M P + 10 (1-P))."""
        Ub = self.U[b:b+1]
        sub = Lat.__new__(Lat)
        sub.__dict__.update(self.__dict__)
        sub.B, sub.U = 1, Ub
        V = self.V

        def proj(x):
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


def pearson(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))


def partial_r(rxy, rxz, ryz):
    return (rxy - rxz * ryz) / np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))


def fisher_sig(r, n):
    return 0.5 * np.log((1 + r) / (1 - r)) * np.sqrt(max(n - 4, 1))


def k2_shells(L):
    k = np.fft.fftfreq(L) * L
    kk = np.zeros((L,) * 4)
    for ax in range(4):
        sh = [1] * 4; sh[ax] = L
        kk = kk + (k.reshape(sh) ** 2)
    return kk


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


def gluon_f3(lat):
    """Leg 2 — NON-SPECTRAL coherence per config: IR fraction of the gluon-field power,
    f3 = P(k^2=1) / P(k^2>0), built from the links alone (no FP operator anywhere).
    Definition frozen in RUN11_SPEC.md before data."""
    B, L = lat.B, lat.L
    A = lat.U[..., 1:].reshape(B, 4, L, L, L, L, 3)
    At = torch.fft.fftn(A, dim=(2, 3, 4, 5))
    w = (At.real ** 2 + At.imag ** 2).sum(dim=(1, 6)).cpu().numpy()   # (B, L,L,L,L)
    kk = k2_shells(L)
    out = np.zeros(B)
    for b in range(B):
        wb = w[b]
        out[b] = wb[kk == 1].sum() / wb[kk > 0].sum()
    return out


def psi_fir_k3(lat, b, h=1e-4):
    """Per config: (lam_min, f1, f2) — f1 = IR fraction of psi_min, f2 = mean IR
    fraction of psi_2 and psi_3 (same deflated operator, k=3 with vectors)."""
    Ub = lat.U[b:b+1]
    sub = Lat.__new__(Lat)
    sub.__dict__.update(lat.__dict__)
    sub.B, sub.U = 1, Ub
    V, L = lat.V, lat.L

    def proj(x):
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
    vals, vecs = eigsh(A, k=3, which='SA', tol=1e-6, maxiter=60000)
    order = np.argsort(vals)
    kk = k2_shells(L)

    def fir_of(vec):
        p = proj(vec); p /= np.linalg.norm(p)
        ft = np.fft.fftn(p.reshape((L,) * 4 + (3,)), axes=(0, 1, 2, 3))
        w = (np.abs(ft) ** 2).sum(axis=-1); w /= w.sum()
        return float(w[kk == 1].sum())

    f1 = fir_of(vecs[:, order[0]])
    f2 = 0.5 * (fir_of(vecs[:, order[1]]) + fir_of(vecs[:, order[2]]))
    return float(vals[order[0]]), f1, f2


# ---------------- the pre-registered plan ----------------
# (tag, beta, L, n) — chosen by the measured grid to span f_IR 0.30-0.92 with six
# points in the 0.45-0.70 transition band. Matched-coherence pairs across beta for
# Leg 3: (2.3,10)~(2.4,14); (2.4,8)~(2.5,12); (2.6,12)~(2.6,10 from run 6).
PLAN_A = [("T1", 2.3,  8, 40), ("T2", 2.3, 10, 32), ("T3", 2.4,  8, 40),
          ("T4", 2.4, 10, 40), ("T5", 2.4, 12, 24), ("T6", 2.4, 14, 16)]
PLAN_B = [("O1", 2.2, 10, 32), ("O2", 2.2, 14, 16), ("M1", 2.5, 10, 32),
          ("M2", 2.5, 12, 24), ("S1", 2.6,  8, 40), ("S2", 2.6, 12, 24)]

if __name__ == "__main__":
    t0 = time.time()
    PLAN = PLAN_A if PART == "A" else PLAN_B
    if FAST:
        PLAN = [("smoke", 2.4, 6, 4)]
    results = {"part": PART, "gates": {}, "predictions": {}, "ensembles": {}}

    # ---- BLIND: predictions for every ensemble, logged before any measurement ----
    print("\n=== PRE-REGISTERED PREDICTIONS (frozen curve, logged before data) ===")
    for tag, BETA, L, N in PLAN:
        fp = fir_predicted(BETA, L)
        rp = r_predicted(fp)
        results["predictions"][tag] = {"beta": BETA, "L": L, "f_pred": round(fp, 4),
                                       "r_pred": round(rp, 4)}
        print(f"  {tag}  (beta={BETA}, L={L}):  f_pred = {fp:.3f}  ->  r_pred = {rp:+.3f}")

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

    for tag, BETA, L, N in PLAN:
        pr = results["predictions"][tag]
        print(f"\n=== {tag}: beta={BETA}, L={L}, n={N}  "
              f"(pred f={pr['f_pred']}, r={pr['r_pred']:+.3f}) ===", flush=True)
        lat = Lat(L, N)
        for _ in range(500 if not FAST else 40):
            lat.sweep(BETA)
        plq = lat.plaq().cpu().numpy()
        th = lat.landau_fix()
        assert th.max().item() < 1e-8
        D1 = gluon_D_p1(lat)
        f3s = gluon_f3(lat)
        lams = np.full(N, np.nan); f1s = np.full(N, np.nan); f2s = np.full(N, np.nan)
        for b in range(N):
            try:
                lams[b], f1s[b], f2s[b] = psi_fir_k3(lat, b)
            except Exception as e:
                print(f"  cfg {b} discarded ({type(e).__name__})", flush=True)
            if (b + 1) % 8 == 0:
                print(f"  {b+1}/{N}  ({time.time()-t0:.0f}s)", flush=True)
        ok = np.isfinite(lams)
        l, d, pl = lams[ok], D1[ok], plq[ok]
        f1, f2, f3 = f1s[ok], f2s[ok], f3s[ok]
        n = int(ok.sum())
        r = pearson(l, d)
        rp = partial_r(r, pearson(l, pl), pearson(d, pl))
        sig = fisher_sig(rp, n)
        # physical regression: D = a*lam + c per ensemble (a carries units of D/lam)
        a_slope = float(np.cov(l, d)[0, 1] / np.var(l, ddof=1))
        resid = d - a_slope * l - (d.mean() - a_slope * l.mean())
        a_err = float(np.sqrt(np.sum(resid ** 2) / (n - 2) / np.sum((l - l.mean()) ** 2)))
        results["ensembles"][tag] = {
            "beta": BETA, "L": L, "n": n,
            "f1_mean": float(f1.mean()), "f1_err": float(f1.std(ddof=1) / np.sqrt(n)),
            "f2_mean": float(f2.mean()), "f2_err": float(f2.std(ddof=1) / np.sqrt(n)),
            "f3_mean": float(f3.mean()), "f3_err": float(f3.std(ddof=1) / np.sqrt(n)),
            "lam_mean": float(l.mean()), "D_mean": float(d.mean()),
            "r": r, "partial": rp, "sigma": float(sig),
            "slope_a": a_slope, "slope_a_err": a_err}
        print(f"  {tag}: f1 = {f1.mean():.3f}+-{f1.std(ddof=1)/np.sqrt(n):.3f}"
              f"  f2 = {f2.mean():.3f}  f3 = {f3.mean():.4f}"
              f"  partial r = {rp:+.4f} ({sig:+.2f} sigma, n={n})"
              f"  a = {a_slope:.3g}+-{a_err:.2g}", flush=True)

    print("\n  Verdict logic lives in RUN11_SPEC.md (F1-F3). Analysis after both parts.")
    out = f"results11{PART}.json" if not FAST else "results11_smoke.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nwrote {out}   total {time.time() - t0:.0f}s")
