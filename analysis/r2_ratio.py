"""
ITEM 4 — is the 0.82 ratio a second disorder constant?

step10 (2026-08-11) established: tr(Ric)/tr(T) = G-weighted mean of diag P_H; unweighted
mean is EXACTLY 3/4 (rank identity); the excess ~0.07 is entirely a G-H correlation.
Left unexplained. Item 1's closure (K*V = Haar constant) suggests the same origin here.

PRE-REGISTERED
 (R1) Haar value == thermalized values within errors, beta-flat
      -> item 4 CLOSED: the 0.82 is the SECOND disorder constant; the G-H correlation is
         a property of Haar randomness (Weingarten-computable in principle).
 (R2) beta-dependent -> dynamic content; map the trend, keep open.
L=3 (n=6) and L=4 (n=4); beta = Haar, 2.0, 2.4, 2.8. Entropy seeding.
"""
import faulthandler; faulthandler.enable()
import numpy as np, sampler_fix
from step3_lattice import Lattice
from sampler_fix import sweep_correct
from step6_scheme import geom
from step8d_softmode import ops_from_G
rng = np.random.default_rng()

def ratio(lat):
    Delta, Dp, E = geom(lat)
    Ric, T = ops_from_G(lat, Dp, E)
    return float(np.trace(Ric) / np.trace(T))

for L, n in ((3, 6), (4, 4)):
    for tag in ("haar", 2.0, 2.4, 2.8):
        vals = []
        for c in range(n):
            sampler_fix.rng = np.random.default_rng()
            lat = Lattice(L, 4)
            if tag == "haar":
                q = rng.normal(size=(4, lat.V, 4))
                lat.U = q / np.linalg.norm(q, axis=-1, keepdims=True)
            else:
                for _ in range(400):
                    sweep_correct(lat, tag)
            vals.append(ratio(lat))
        print(f"  L={L} {str(tag):>5}: ratio = {np.mean(vals):.5f} "
              f"+- {np.std(vals, ddof=1)/np.sqrt(n):.5f}", flush=True)
print("R1 Haar==thermal, beta-flat -> second disorder constant, item 4 CLOSED")
print("R2 beta-dependent -> dynamic; map it")
