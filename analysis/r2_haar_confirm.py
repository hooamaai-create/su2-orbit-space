"""Item-1 closure checks. (K1) Haar flat over L=3..6 within ~1%. (K2) K*V(beta) rides
monotonically from the Haar value toward the plateau with no detour. Either failing
reopens the disorder explanation."""
import faulthandler; faulthandler.enable()
import numpy as np, sampler_fix
from step3_lattice import Lattice
from sampler_fix import sweep_correct
from mf_kv import kv_mf
rng = np.random.default_rng()

print("K1: Haar value vs volume")
for L, npl, nc in ((3, 300, 4), (4, 250, 4), (5, 150, 3), (6, 100, 2)):
    vals = []
    for c in range(nc):
        lat = Lattice(L, 4)
        q = rng.normal(size=(4, lat.V, 4))
        lat.U = q / np.linalg.norm(q, axis=-1, keepdims=True)
        vals.append(kv_mf(lat, npl, seed=int(rng.integers(1e6))).mean() * lat.V)
    print(f"  L={L}: {np.mean(vals):.6f} +- {np.std(vals,ddof=1)/np.sqrt(nc):.6f}", flush=True)

print("K2: interpolation Haar -> plateau, L=4")
for beta in (0.1, 0.5, 1.0, 1.5, 2.0):
    vals = []
    for c in range(3):
        sampler_fix.rng = np.random.default_rng()
        lat = Lattice(L=4, D=4) if False else Lattice(4, 4)
        q = rng.normal(size=(4, lat.V, 4))
        lat.U = q / np.linalg.norm(q, axis=-1, keepdims=True)   # hot start
        for _ in range(300):
            sweep_correct(lat, beta)
        vals.append(kv_mf(lat, 200, seed=int(rng.integers(1e6))).mean() * lat.V)
    print(f"  beta={beta:>4.1f}: {np.mean(vals):.6f} +- {np.std(vals,ddof=1)/np.sqrt(3):.6f}",
          flush=True)
print("ref: Haar 0.0795-0.0797, plateau(2.4) 0.0801")
