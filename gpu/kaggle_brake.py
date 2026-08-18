"""
KAGGLE GPU RUN 8 — WHERE DOES THE FALL BRAKE?  (conjecture #12's test)

Two measured facts cannot both hold forever: our window (L=6-12) shows the vacuum falling
to the horizon FASTER than kinematics (p ~ -2.3 at beta=2.2); Cucchieri-Mendes at L=64-128
(same coupling) found it SLOWER (p ~ -1.5). Somewhere between, the local exponent must
cross the kinematic value: a TURNING SIZE ell*. Conjecture #12: the brake is the coherent
IR wave fragmenting when the box outgrows the medium's correlation length — in which case
ell* is a physical length, possibly tied to the mass scale.

PART O  lambda_min AND f_IR of the lowest deflated FP mode, beta = 2.2, L = 8..18.
        Local exponents p_eff from overlapping L-triplets.
         (T1) p_eff falls through the free value in-window -> THE TURN IS LOCATED;
              quote ell* and compare against 1/m (0.28/sqrt(sigma)), 1/sqrt(sigma),
              and the CM window (~35+).
         (T2) p_eff stays super-kinematic through L=18 -> no turn in reach; ell* > ~8
              in sqrt(sigma) units, and the naive ell* ~ 1/m version of #12 is REFUTED.
         (T3) MECHANISM: if f_IR falls with L where p_eff falls -> the brake is wave
              fragmentation (mechanism confirmed). If p turns but f_IR stays high ->
              the mechanism claim dies even if the length lives.
Gates: G1 plaquette anchor, G3 deflated free-field. ~6-9 h full. Writes results8.json.
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



def psi_min(lat, b, h=1e-4):
    Ub = lat.U[b:b+1]
    sub = Lat.__new__(Lat)
    sub.__dict__.update(lat.__dict__)
    sub.B, sub.U = 1, Ub
    V = lat.V

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
    p = proj(vec[:, 0])
    p /= np.linalg.norm(p)
    return float(val[0]), p.reshape(V, 3)


def f_ir_of(L, psi):
    ft = np.fft.fftn(psi.reshape((L,) * 4 + (3,)), axes=(0, 1, 2, 3))
    w = (np.abs(ft) ** 2).sum(axis=-1)
    w /= w.sum()
    k = np.fft.fftfreq(L) * L
    kk = np.zeros((L,) * 4)
    for ax in range(4):
        sh = [1] * 4
        sh[ax] = L
        kk = kk + (k.reshape(sh) ** 2)
    return float(w[kk == 1].sum()), float(w[kk == 0].sum())


if __name__ == "__main__":
    t0 = time.time()
    BETA = 2.2
    ASQS = 0.561                      # a*sqrt(sigma) at beta=2.2 (chi33, L=16, run 6)
    results = {"gates": {}, "partO": {}}

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

    PLAN = ([(8, 12), (10, 10), (12, 8), (14, 7), (16, 6), (18, 5)]
            if not FAST else [(4, 2), (6, 2)])
    Ls, lam_m, lam_e, fir_m, fir_e = [], [], [], [], []
    print(f"\n=== PART O: lambda_min and f_IR vs L at beta = {BETA} ===", flush=True)
    for L, n in PLAN:
        lat = Lat(L, n)
        for _ in range(500 if not FAST else 40):
            lat.sweep(BETA)
        pq = lat.plaq().mean().item()
        th = lat.landau_fix()
        assert th.max().item() < 1e-8, f"landau nonconvergence at L={L}"
        lv, fv = [], []
        for b in range(n):
            try:
                val, psi = psi_min(lat, b)
            except Exception as e:
                print(f"    L={L} cfg {b} discarded ({type(e).__name__})", flush=True)
                continue
            fir, f0 = f_ir_of(L, psi)
            assert f0 < 1e-3, "deflation leakage"
            lv.append(val)
            fv.append(fir)
        Ls.append(L)
        lam_m.append(float(np.mean(lv)))
        lam_e.append(float(np.std(lv, ddof=1) / np.sqrt(len(lv))))
        fir_m.append(float(np.mean(fv)))
        fir_e.append(float(np.std(fv, ddof=1) / np.sqrt(len(fv))))
        results["partO"][L] = {"lam": lam_m[-1], "lam_err": lam_e[-1],
                               "f_ir": fir_m[-1], "f_ir_err": fir_e[-1],
                               "ell": round(L * ASQS, 2), "n": len(lv), "plaq": pq}
        print(f"  L={L} (ell={L*ASQS:.1f}) n={len(lv)} lam={lam_m[-1]:.5f}"
              f"+-{lam_e[-1]:.5f}  f_IR={fir_m[-1]:.3f}+-{fir_e[-1]:.3f}"
              f"  ({time.time()-t0:.0f}s)", flush=True)

    print("\n  local exponents from consecutive L-triplets:", flush=True)
    for i in range(len(Ls) - 2):
        x = np.log(Ls[i:i+3]); y = np.log(lam_m[i:i+3])
        w = np.array(lam_e[i:i+3]) / np.array(lam_m[i:i+3])
        W = 1 / w ** 2
        den = W.sum() * (W * x * x).sum() - (W * x).sum() ** 2
        pe = (W.sum() * (W * x * y).sum() - (W * x).sum() * (W * y).sum()) / den
        epe = np.sqrt(W.sum() / den)
        pfree = np.polyfit(x, np.log([4 * np.sin(np.pi / l) ** 2 for l in Ls[i:i+3]]), 1)[0]
        mid_ell = Ls[i + 1] * ASQS
        results["partO"][f"win{Ls[i]}-{Ls[i+2]}"] = {
            "p_eff": pe, "err": epe, "p_free": pfree, "mid_ell": round(mid_ell, 2)}
        print(f"    L={Ls[i]}..{Ls[i+2]} (ell~{mid_ell:.1f}): p_eff = {pe:+.3f}+-{epe:.3f}"
              f"   free {pfree:+.3f}   excess {(pe-pfree):+.3f}", flush=True)

    print("\n  scale table (sqrt(sigma) units): 1/m(0++) = 0.28 | 1/sqrt(sigma) = 1 | "
          "our window 4.5-10 | CM regime 35+")
    print("  T1 p_eff crosses free in-window -> turn LOCATED, quote ell*")
    print("  T2 stays super-kinematic       -> ell* > ~8; naive ell*~1/m REFUTED")
    print("  T3 f_IR falls where p_eff falls -> fragmentation mechanism CONFIRMED;")
    print("     p turns with f_IR flat       -> mechanism dies even if length lives")

    with open("results8.json", "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nwrote results8.json   total {time.time() - t0:.0f}s")
