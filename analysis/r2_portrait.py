"""
THE PHASE PORTRAIT — every IR measurement of the programme in the (proximity, coherence)
plane, and the first fit of the pull's law as a function on that plane.

CLAIM UNDER TEST (formulated 2026-08-17 from elimination): the state of the vacuum near
the wall is captured by TWO coordinates — proximity s = lam/lam_free and carrier
coherence f_IR — and the pull's excess e = p_eff - p_free is a function F(s, f_IR),
whose zero-contour is the braking surface. Every single-variable story died; this is the
minimal surviving structure.

MODELS (weighted least squares, weights 1/err^2 on e):
  M0: e = a                      (nothing)
  M1a: e = a + b*f_IR            (coherence only)
  M1b: e = a + b*ln(s)           (proximity only)
  M2: e = a + b*f_IR + c*ln(s)   (the two-coordinate law)
RULE: M2 is claimable only if it beats BOTH M1's by delta-chi2 > 6 (2 sigma-ish for one
extra parameter). Otherwise report the best single variable, or noise.

Also: the r-collapse (law strength vs f_IR alone) tested on the 5 measured r points, and
the killer prediction evaluated at (beta=2.4, L=16).
"""
import numpy as np

# ---- window data: (beta, midL, p_eff, err, p_free) ----
WIN = [
    (2.2, 10, -2.703, 0.507, -1.929), (2.2, 12, -2.242, 0.613, -1.952),
    (2.2, 14, -1.446, 0.794, -1.965), (2.2, 16, -1.819, 0.965, -1.973),
    (2.4, 10, -2.029, 0.621, -1.929), (2.4, 12, -2.943, 0.819, -1.952),
    (2.4, 14, -4.456, 1.370, -1.965), (2.4, 16, -2.061, 1.373, -1.973),
    (2.6, 10, -2.082, 0.223, -1.929), (2.6, 12, -1.771, 0.319, -1.952),
    (2.6, 14, -1.871, 0.356, -1.965), (2.6, 16, -1.591, 0.223, -1.973),
]
# ---- per-(beta,L) grids from runs 8+9: lam and f_IR ----
LAM = {(2.2, 8): .04037, (2.2, 10): .02105, (2.2, 12): .01374, (2.2, 14): .009939,
       (2.2, 16): .009045, (2.2, 18): .006236,
       (2.4, 8): .09339, (2.4, 10): .05531, (2.4, 12): .04319, (2.4, 14): .01883,
       (2.4, 16): .01265, (2.4, 18): .01127,
       (2.6, 8): .22137, (2.6, 10): .12099, (2.6, 12): .09535, (2.6, 14): .06905,
       (2.6, 16): .05667, (2.6, 18): .04636}
FIR = {(2.2, 8): .3806, (2.2, 10): .3829, (2.2, 12): .3839, (2.2, 14): .3008,
       (2.2, 16): .3564, (2.2, 18): .2802,
       (2.4, 8): .687, (2.4, 10): .622, (2.4, 12): .560, (2.4, 14): .464,
       (2.4, 16): .444, (2.4, 18): .409,
       (2.6, 8): .916, (2.6, 10): .829, (2.6, 12): .868, (2.6, 14): .883,
       (2.6, 16): .869, (2.6, 18): .860}


def lam_free(L):
    return 4 * np.sin(np.pi / L) ** 2


rows = []
for beta, mL, pe, err, pf in WIN:
    s = LAM[(beta, mL)] / lam_free(mL)
    fir = FIR[(beta, mL)]
    rows.append((beta, mL, np.log(s), fir, pe - pf, err))
rows = np.array(rows)
b_, L_, lns, fir, e, err = rows.T
W = 1 / err ** 2

print("PHASE PORTRAIT (12 windows):")
print(f"{'beta':>5} {'midL':>5} {'ln s':>7} {'f_IR':>6} {'excess':>8} {'err':>6}")
for r in rows:
    print(f"{r[0]:>5.1f} {int(r[1]):>5} {r[2]:>7.2f} {r[3]:>6.2f} {r[4]:>8.2f} {r[5]:>6.2f}")


def wfit(X, y, w):
    Xw = X * np.sqrt(w[:, None]); yw = y * np.sqrt(w)
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    chi2 = float(((y - X @ beta) ** 2 * w).sum())
    return beta, chi2


one = np.ones_like(e)
models = {
    "M0 const":        np.column_stack([one]),
    "M1a f_IR":        np.column_stack([one, fir]),
    "M1b ln s":        np.column_stack([one, lns]),
    "M2 f_IR + ln s":  np.column_stack([one, fir, lns]),
}
print("\nMODEL COMPARISON (12 points):")
res = {}
for name, X in models.items():
    beta_c, chi2 = wfit(X, e, W)
    res[name] = (beta_c, chi2)
    dof = len(e) - X.shape[1]
    print(f"  {name:<15} chi2 = {chi2:6.2f}  /dof = {chi2/dof:5.2f}   "
          f"coeffs = {np.array2string(beta_c, precision=3)}")

d_a = res["M1a f_IR"][1] - res["M2 f_IR + ln s"][1]
d_b = res["M1b ln s"][1] - res["M2 f_IR + ln s"][1]
print(f"\n  M2 vs coherence-only: delta chi2 = {d_a:.2f}")
print(f"  M2 vs proximity-only: delta chi2 = {d_b:.2f}")
print("  RULE: M2 claimable only if BOTH > 6.")
a2, b2, c2 = res["M2 f_IR + ln s"][0]
print(f"\n  M2 surface: excess = {a2:+.2f} {b2:+.2f}*f_IR {c2:+.2f}*ln(s)")
print(f"  braking contour (excess = 0): f_IR = {-a2/b2:+.2f} {-c2/b2:+.2f}*ln(s)")

# ---- the r-collapse: law strength vs f_IR alone ----
print("\nR-COLLAPSE (law strength vs coherence):")
RPTS = [(0.383, +0.059), (0.622, -0.502), (0.622, -0.437), (0.930, -0.511),
        (0.860, -0.511)]   # (f_IR at the r-measurement point, partial r)
for f, r in RPTS:
    print(f"  f_IR = {f:.2f}  ->  r = {r:+.3f}")
print("  activation near f_IR ~ 0.5, saturation |r| ~ 0.51")
print("\nKILLER PREDICTION at (beta=2.4, L=16): f_IR = 0.444 -> BELOW activation ->")
print("  the law should be OFF (|r| < 0.2) at a coupling where it is strong at L=10.")
print("  One n~48 GPU run decides the whole two-coordinate formulation.")
