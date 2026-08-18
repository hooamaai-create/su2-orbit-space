"""
MATRIX-FREE  K * V  TO LARGE VOLUME.

WHY.  kv_controls.py (C2) left the decisive question open:

      L        2        3        4        5      6
      inter  0.0837   0.0804   0.0803    ?      ?
      free   0.0709   0.0773   0.0803  0.0820  0.0830     (analytic, still rising)

The interacting curve is falling-then-flat; the free-field curve is rising through it.
They cross at L=4.  If the interacting curve stays flat past L=4 while the free-field
one keeps rising, K*V is a genuine coupling-dependent plateau.  If it turns and tracks
the free-field rise, "K*V = const" is a small-volume accident and dies.

The dense path (step3_lattice.geometry) forms a 3DV x 3DV projector and cannot reach
L >= 5.  Here everything is applied matrix-free:

      P_H z = z - D (Delta^+ D^dag z)          one CG solve
      K     = 3 <rho, Delta^+ rho>             one CG solve,  rho = sum_mu x_mu X y_mu

VALIDATION GATE (must pass before any L>=5 number is quoted).
      matrix-free K*V at L=3, beta=2.4 must agree with the dense value 0.080408 +- 0.000402
      to within 1%.  If not, the implementation is wrong and no result is reported.

PRE-REGISTERED RULE
   (R1) interacting K*V flat (within errors) over L = 3,4,5,6 while analytic rises >7%
        -> the plateau is real and coupling-generated.
   (R2) interacting K*V rises with L, tracking the analytic
        -> "K*V = const" was a two-volume accident.  Dead.
   (R3) interacting K*V falls with L
        -> K decays faster than 1/V.  Also dead, differently.
"""
import numpy as np
from step3_lattice import Lattice, adjoint_matrix
from kv_constant import kv_analytic


def ops(lat):
    D, V = lat.D, lat.V
    R = np.stack([adjoint_matrix(lat.U[mu]) for mu in range(D)])

    def Dop(w):
        w = w.reshape(V, 3)
        return np.stack([np.einsum('vij,vj->vi', R[mu], w[lat.fwd[mu]]) - w
                         for mu in range(D)])

    def Ddag(b):
        b = b.reshape(D, V, 3)
        out = np.zeros((V, 3))
        for mu in range(D):
            out += np.einsum('vji,vj->vi', R[mu][lat.bwd[mu]], b[mu][lat.bwd[mu]]) - b[mu]
        return out

    def Delta(w):
        return Ddag(Dop(w)).ravel()

    def cg(rhs, tol=1e-11, itmax=20000):
        """CG on Delta.  RAISES rather than returning a silently wrong answer.

        Delta is singular at U = 1 (global colour rotations are exact zero modes).  For
        rhs = D^dag z the right-hand side lies in range(D^dag) perp ker(Delta) and CG from
        x0 = 0 converges to the pseudo-inverse solution.  But rho = sum_mu x_mu X y_mu has
        no such guarantee: at exact U = 1 it overlaps ker(Delta) and CG diverges, previously
        returning K ~ 5e32 after 13435 iterations with residual 5e19 and no complaint.
        Found by adversarial audit 2026-08-12.  Every published script thermalises first and
        so never hit this path, but an unguarded solver is a bug regardless of who calls it.
        """
        b = np.asarray(rhs).ravel()
        nb = np.linalg.norm(b)
        if nb == 0.0:
            return np.zeros_like(b)
        x = np.zeros_like(b); r = b.copy(); p = r.copy(); rs = r @ r
        converged = False
        for _ in range(itmax):
            Ap = Delta(p)
            pAp = p @ Ap
            if pAp <= 0:
                raise FloatingPointError(
                    f"CG breakdown: p.Delta.p = {pAp:.3e} <= 0. Delta is not positive "
                    f"definite on the Krylov space; rhs probably overlaps ker(Delta).")
            al = rs / pAp
            x += al * p; r -= al * Ap
            rsn = r @ r
            if np.sqrt(rsn) / nb < tol:
                converged = True
                break
            p = r + (rsn / rs) * p; rs = rsn
        if not converged:
            raise FloatingPointError(
                f"CG failed to converge in {itmax} iterations: relative residual "
                f"{np.sqrt(r @ r)/nb:.3e} > tol {tol:.1e}. Do NOT use this result.")
        # true residual, not the recurrence estimate (which drifts on ill-conditioned systems)
        true_res = np.linalg.norm(Delta(x) - b) / nb
        if true_res > 1e-6:
            raise FloatingPointError(
                f"CG recurrence residual converged but TRUE residual is {true_res:.3e}. "
                f"The recurrence has drifted; the solution is wrong.")
        return x

    def P_H(z):
        z = z.reshape(D, V, 3)
        return z - Dop(cg(Ddag(z)))

    return Dop, Ddag, cg, P_H


def kv_mf(lat, nplanes, seed=0):
    _, _, cg, P_H = ops(lat)
    rng = np.random.default_rng(seed)
    D, V = lat.D, lat.V
    Ks = np.empty(nplanes)
    for j in range(nplanes):
        x = P_H(rng.normal(size=(D, V, 3)))
        y = P_H(rng.normal(size=(D, V, 3)))
        x = x / np.linalg.norm(x)
        y = y - (x * y).sum() * x
        y = y / np.linalg.norm(y)
        rho = np.cross(x, y).sum(axis=0)
        Ks[j] = 3.0 * rho.ravel() @ cg(rho)
    return Ks


if __name__ == "__main__":
    D, BETA = 4, 2.4
    DENSE_L3 = 0.080408

    print("=" * 72)
    print("VALIDATION GATE: matrix-free vs dense at L=3, beta=2.4")
    print("=" * 72)
    lat = Lattice(3, D)
    for _ in range(300):
        lat.metropolis_sweep(BETA)
    Ks = kv_mf(lat, 300, seed=901)
    mf3 = Ks.mean() * lat.V
    print(f"  dense       0.080408 +- 0.000402")
    print(f"  matrix-free {mf3:.6f}    ratio {mf3/DENSE_L3:.4f}")
    ok = abs(mf3 / DENSE_L3 - 1) < 0.01
    print(f"  GATE {'PASS' if ok else 'FAIL — no L>=5 number will be quoted'}")
    if not ok:
        raise SystemExit(1)

    print("\n" + "=" * 72)
    print(f"K*V vs VOLUME at beta = {BETA}   [R1 flat / R2 rises / R3 falls]")
    print("=" * 72)
    print(f"{'L':>3} {'V':>6} {'<plaq>':>8} {'n':>3} {'mean K*V':>21} {'analytic':>10} {'ratio':>7}")
    for L, ncfg, npl in ((3, 3, 300), (4, 3, 250), (5, 2, 150), (6, 2, 100)):
        ana = kv_analytic(L, D)[0]
        vals, plaqs = [], []
        for c in range(ncfg):
            lt = Lattice(L, D)
            for _ in range(300):
                lt.metropolis_sweep(BETA)
            vals.append(kv_mf(lt, npl, seed=2000 + 13 * L + c).mean() * lt.V)
            plaqs.append(lt.mean_plaquette())
        m = float(np.mean(vals))
        s = float(np.std(vals, ddof=1) / np.sqrt(len(vals))) if len(vals) > 1 else 0.0
        print(f"{L:>3} {L**D:>6} {np.mean(plaqs):>8.4f} {ncfg:>3} "
              f"{m:>14.6f} +- {s:<6.6f} {ana:>10.6f} {m/ana:>7.4f}")

    print("\n  1/(4 pi) = %.7f" % (1 / (4 * np.pi)))
    print("  R1 flat while analytic rises -> real plateau")
    print("  R2 tracks analytic rise      -> two-volume accident, dead")
    print("  R3 falls                     -> decays faster than 1/V, dead")
