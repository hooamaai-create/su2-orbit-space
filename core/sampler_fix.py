"""
BUG #10 — THE MONTE CARLO SAMPLER VIOLATES DETAILED BALANCE.  Fix and validation.

THE BUG (step3_lattice.Lattice.metropolis_sweep).
    S = self.staple_sum(mu)          # computed ONCE
    for _ in range(hits):
        ... update ALL V links of direction mu simultaneously against S ...
staple_sum(mu) contains U[mu][fwd[nu]] -- links of the SAME direction mu at neighbouring
sites.  Those links are each other's environment.  Updating them simultaneously against a
frozen staple field, and then re-using that stale field for further hits, breaks detailed
balance.  The chain therefore samples no Boltzmann distribution at all: <plaq> at beta=2.4
came out ~0.588 against the literature value ~0.630.

THE FIX.  Update only a set of links that are mutually non-interacting, and recompute the
staples before each such set.  Two mu-links at sites n and m interact iff n - m = +-e_nu for
some nu != mu.  So a valid colouring must ensure that no two sites in a class differ by one
unit step in any direction other than mu.

  - all extents orthogonal to mu EVEN:  parity of sum_{nu != mu} n_nu.  2 classes.
  - otherwise: the tuple ( n_nu mod k_nu )_{nu != mu} with k_nu = 2 if L_nu is even else 3.
    Two sites in one class agree mod k_nu in every nu != mu; a +-1 step is nonzero mod k_nu
    since k_nu >= 2, so no two members of a class interact.  Valid for ANY extents.

The odd case costs up to 3^(D-1) staple recomputations per hit, which is why the naive
parity colouring is used whenever it is legal.  Parity colouring is ILLEGAL for odd L
(periodic: sites L-1 and 0 are neighbours and share parity), which is exactly the L=3 and
L=5 lattices used throughout this programme.

VALIDATION (pre-registered).
 (V1) SU(2) Wilson <plaq> on 4^4 at beta = 2.4 must come out ~0.630 (literature), not 0.588.
 (V2) the corrected sampler must agree with an independent Creutz SU(2) HEAT BATH, which is
      exact and shares no code with the Metropolis path, to within statistical error.
 (V3) strong-coupling check: <plaq> -> beta/4 as beta -> 0 for SU(2) (leading order).
 (V4) the corrected sampler must be independent of `hits`, since each hit now sees fresh
      staples.  The buggy one is not.
If any of these fails the fix is not adopted and nothing is remeasured.
"""
import numpy as np
from step3_lattice import Lattice, qmul, qconj, qrand_near_identity

rng = np.random.default_rng(20260812)


def colour_classes(shape, mu):
    """Index arrays of mutually non-interacting sites for updating direction mu."""
    D = len(shape)
    idx = np.arange(int(np.prod(shape))).reshape(shape)
    coords = np.indices(shape)
    others = [nu for nu in range(D) if nu != mu]
    if all(shape[nu] % 2 == 0 for nu in others):
        key = sum(coords[nu] for nu in others) % 2
    else:
        key = np.zeros(shape, dtype=np.int64)
        for nu in others:
            k = 2 if shape[nu] % 2 == 0 else 3
            key = key * k + (coords[nu] % k)
    out = []
    for c in np.unique(key):
        out.append(idx[key == c].ravel())
    return out


def sweep_correct(lat, beta, eps=0.35, hits=2):
    """Metropolis with fresh staples per colour class.  Satisfies detailed balance."""
    if not hasattr(lat, '_classes'):
        lat._classes = [colour_classes(lat.shape, mu) for mu in range(lat.D)]
    acc = tot = 0
    for mu in range(lat.D):
        for _ in range(hits):
            for cls in lat._classes[mu]:
                S = lat.staple_sum(mu)[cls]
                U = lat.U[mu][cls]
                r = qrand_near_identity((len(cls),), eps)
                Un = qmul(r, U)
                dS = -beta * (qmul(Un, S)[..., 0] - qmul(U, S)[..., 0])
                take = (dS <= 0) | (rng.random(len(cls)) < np.exp(-np.clip(dS, 0, 50)))
                lat.U[mu][cls] = np.where(take[:, None], Un, U)
                acc += int(take.sum()); tot += len(cls)
    return acc / tot


def heatbath_sweep(lat, beta):
    """Creutz SU(2) heat bath. Exact; shares no acceptance logic with Metropolis."""
    if not hasattr(lat, '_classes'):
        lat._classes = [colour_classes(lat.shape, mu) for mu in range(lat.D)]
    for mu in range(lat.D):
        for cls in lat._classes[mu]:
            S = lat.staple_sum(mu)[cls]
            k = np.sqrt(np.maximum((S ** 2).sum(axis=-1), 1e-300))       # |S|
            Vd = qconj(S) / k[:, None]                                    # S^dag/|S| in SU(2)
            a = beta * k
            n = len(cls)
            a0 = np.empty(n); todo = np.arange(n)
            while len(todo) > 0:
                x = rng.random(len(todo))
                lam = 1.0 + np.log(np.maximum(x, 1e-300)) / np.maximum(a[todo], 1e-300)
                ok = (rng.random(len(todo)) ** 2) <= (1.0 - lam ** 2 / 1.0).clip(min=0) + 1
                ok &= (np.abs(lam) <= 1.0)
                acc2 = rng.random(len(todo)) <= np.sqrt(np.clip(1 - lam ** 2, 0, 1))
                good = ok & acc2
                a0[todo[good]] = lam[good]
                todo = todo[~good]
            rad = np.sqrt(np.clip(1 - a0 ** 2, 0, 1))
            ct = 2 * rng.random(n) - 1
            st = np.sqrt(np.clip(1 - ct ** 2, 0, 1))
            ph = 2 * np.pi * rng.random(n)
            vec = np.stack([a0, rad * st * np.cos(ph), rad * st * np.sin(ph), rad * ct], -1)
            lat.U[mu][cls] = qmul(vec, Vd)
    return 1.0


def measure(sweeper, L, beta, D=4, nther=600, nmeas=40, gap=4, **kw):
    lat = Lattice(L, D)
    for _ in range(nther):
        sweeper(lat, beta, **kw)
    vals = []
    for _ in range(nmeas):
        for _ in range(gap):
            sweeper(lat, beta, **kw)
        vals.append(lat.mean_plaquette())
    v = np.array(vals)
    return v.mean(), v.std(ddof=1) / np.sqrt(len(v))


def buggy(lat, beta, **kw):
    return lat.metropolis_sweep(beta, **kw)


if __name__ == "__main__":
    print("=" * 78)
    print("V1 + V2:  4^4 SU(2) Wilson.  literature <plaq> at beta=2.4 is ~0.630")
    print("=" * 78)
    print(f"{'beta':>6} {'BUGGY (in use)':>22} {'FIXED metropolis':>22} {'heat bath':>22}",
          flush=True)
    for beta in (2.0, 2.4, 2.8):
        b = measure(buggy, 4, beta)
        f = measure(sweep_correct, 4, beta)
        h = measure(heatbath_sweep, 4, beta)
        print(f"{beta:>6.1f} {b[0]:>15.5f} +- {b[1]:<5.5f} {f[0]:>15.5f} +- {f[1]:<5.5f} "
              f"{h[0]:>15.5f} +- {h[1]:<5.5f}", flush=True)

    print("\n" + "=" * 78)
    print("V3: strong coupling.  SU(2) leading order  <plaq> -> beta/4")
    print("=" * 78)
    for beta in (0.2, 0.4):
        f = measure(sweep_correct, 4, beta, nther=300, nmeas=25)
        print(f"  beta={beta:>4.1f}  fixed = {f[0]:.5f} +- {f[1]:.5f}   beta/4 = {beta/4:.5f}")

    print("\n" + "=" * 78)
    print("V4: independence of `hits` (fresh staples each hit)")
    print("=" * 78)
    for h_ in (1, 2, 4):
        f = measure(sweep_correct, 4, 2.4, nther=500, nmeas=30, hits=h_)
        b = measure(buggy, 4, 2.4, nther=500, nmeas=30, hits=h_)
        print(f"  hits={h_}   fixed = {f[0]:.5f} +- {f[1]:.5f}      buggy = {b[0]:.5f} "
              f"+- {b[1]:.5f}")

    print("\n" + "=" * 78)
    print("ODD L: parity colouring is ILLEGAL there; check the mod-3 fallback")
    print("=" * 78)
    for L in (3, 5):
        cls = colour_classes((L,) * 4, 0)
        print(f"  L={L}: {len(cls)} colour classes, sizes {sorted(set(len(c) for c in cls))}")
        f = measure(sweep_correct, L, 2.4, nther=400, nmeas=25)
        h = measure(heatbath_sweep, L, 2.4, nther=400, nmeas=25)
        print(f"        fixed = {f[0]:.5f} +- {f[1]:.5f}   heat bath = {h[0]:.5f} +- {h[1]:.5f}")
