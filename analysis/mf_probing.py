"""
DETERMINISTIC PROBING (colouring) for the diagonal blocks — replaces random Hutchinson.

Random probes gave 18% error on min(r) at N=256, against a physical spread of 5%. Useless.

Colouring: assign colour(n) = n mod c componentwise, giving c^4 spatial colours. A probe
vector supported on one colour class picks up, at site n, the true diagonal PLUS terms
G(n,m) for same-coloured m, all of which are at distance >= c. So the estimate is EXACT up
to O(G(c)) -- no statistical noise at all, only a controlled truncation that shrinks as c
grows. Cost: 3 * c^4 solves (3 colour components), independent of N.
"""
import numpy as np
from mf_ricci import MF


def probing_diag(mf, c=2):
    D, V, L = mf.D, mf.V, mf.lat.L
    coords = np.stack(np.unravel_index(np.arange(V), (L,) * D), axis=-1)   # (V,D)
    lab = np.zeros(V, dtype=int)
    for d in range(D):
        lab = lab * c + (coords[:, d] % c)
    ncol = c ** D
    Gd = np.zeros((V, 3, 3))
    Hd = np.zeros((D, V, 3, 3))
    for k in range(ncol):
        mask = (lab == k).astype(float)
        for b in range(3):
            z0 = np.zeros((V, 3)); z0[:, b] = mask
            x = mf.cg(z0.ravel()).reshape(V, 3)
            Gd[:, :, b] += x * mask[:, None]

            # BUG FIX 2026-08-10: the probe previously set component b on ALL D link
            # directions at once, so the result was sum_mu' P_H[(mu,n,a),(mu',n,b)] rather
            # than the diagonal block. Probes must be separated per (mu, b).
            for mu in range(D):
                z1 = np.zeros((D, V, 3)); z1[mu, :, b] = mask
                y = mf.P_H(z1.ravel()).reshape(D, V, 3)
                Hd[mu, :, :, b] += y[mu] * mask[:, None]
    Gd = 0.5 * (Gd + np.transpose(Gd, (0, 2, 1)))
    Hd = 0.5 * (Hd + np.transpose(Hd, (0, 1, 3, 2)))
    return Gd, Hd, 3 * ncol * (1 + D)
