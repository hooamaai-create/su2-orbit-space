"""
KAGGLE GPU RUN 10 — IS COHERENCE THE ORDER PARAMETER OF THE LAW?

The r-collapse (2026-08-17): the lambda-D law's strength reads off the carrier's
coherence f_IR alone — activation near 0.5, saturation |r| ~ 0.51. Five points, no
counterexample. PRE-REGISTERED PREDICTIONS, made before this run:

  P-OFF  (beta=2.4, L=16): f_IR = 0.444 (below activation) -> the law must be OFF,
         |r| < 0.2 — at the SAME coupling where it runs -0.44/-0.50 at L=10.
  P-ON   (beta=2.6, L=10): f_IR = 0.829 (saturated) -> the law must be ON, r ~ -0.5.

  (V1) both predictions hit -> f_IR is the order parameter; the collapse is a law.
  (V2) P-OFF fails (law still on at L=16) -> the collapse is broken; f_IR is not the
       controlling variable; conjecture joins the dead at #14.
  (V3) mixed/marginal -> report; the collapse stays a hypothesis.
Partial r (plaquette confound removed) throughout; n = 48 / 40. ~5-7 h. results10.json.
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


def psi_and_fir(lat, b, h=1e-4):
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
    val, vec = eigsh(A, k=1, which='SA', tol=1e-6, maxiter=60000)
    p = proj(vec[:, 0]); p /= np.linalg.norm(p)
    ft = np.fft.fftn(p.reshape((L,) * 4 + (3,)), axes=(0, 1, 2, 3))
    w = (np.abs(ft) ** 2).sum(axis=-1); w /= w.sum()
    k = np.fft.fftfreq(L) * L
    kk = np.zeros((L,) * 4)
    for ax in range(4):
        sh = [1] * 4; sh[ax] = L
        kk = kk + (k.reshape(sh) ** 2)
    return float(val[0]), float(w[kk == 1].sum())


if __name__ == "__main__":
    t0 = time.time()
    results = {"gates": {}, "partQ": {}}

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

    PLAN = ([("P-OFF", 2.4, 16, 48, "predict |r| < 0.2 (f_IR ~ 0.44)"),
             ("P-ON", 2.6, 10, 40, "predict r ~ -0.5 (f_IR ~ 0.83)")]
            if not FAST else [("smoke", 2.4, 6, 4, "smoke")])
    for tag, BETA, L, N, pred in PLAN:
        print(f"\n=== {tag}: beta={BETA}, L={L}, n={N} — {pred} ===", flush=True)
        lat = Lat(L, N)
        for _ in range(500 if not FAST else 40):
            lat.sweep(BETA)
        plq = lat.plaq().cpu().numpy()
        th = lat.landau_fix()
        assert th.max().item() < 1e-8
        D1 = gluon_D_p1(lat)
        lams = np.full(N, np.nan); firs = np.full(N, np.nan)
        for b in range(N):
            try:
                lams[b], firs[b] = psi_and_fir(lat, b)
            except Exception as e:
                print(f"  cfg {b} discarded ({type(e).__name__})", flush=True)
            if (b + 1) % 12 == 0:
                print(f"  {b+1}/{N}  ({time.time()-t0:.0f}s)", flush=True)
        ok = np.isfinite(lams)
        l, d, pl, f = lams[ok], D1[ok], plq[ok], firs[ok]
        r = pearson(l, d)
        rp = partial_r(r, pearson(l, pl), pearson(d, pl))
        sig = fisher_sig(rp, int(ok.sum()))
        results["partQ"][tag] = {
            "beta": BETA, "L": L, "n": int(ok.sum()),
            "f_ir_mean": float(f.mean()),
            "f_ir_err": float(f.std(ddof=1) / np.sqrt(len(f))),
            "r": r, "partial": rp, "sigma": float(sig)}
        print(f"  {tag}: measured f_IR = {f.mean():.3f}+-{f.std(ddof=1)/np.sqrt(len(f)):.3f}"
              f"   partial r = {rp:+.4f}  ({sig:+.2f} sigma, n={ok.sum()})", flush=True)

    print("\n  V1 both predictions hit  -> f_IR IS the order parameter; collapse = law")
    print("  V2 P-OFF fails           -> collapse broken; conjecture #14 dies")
    print("  V3 mixed/marginal        -> stays a hypothesis; report plainly")

    with open("results10.json", "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nwrote results10.json   total {time.time() - t0:.0f}s")
