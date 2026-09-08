import glob
import os
import re

import numpy as np
import pandas as pd
import pyvista as pv

CASE = os.path.dirname(os.path.abspath(__file__))
Rp_1, Rp_2 = 11.2e-6, 5.6e-6
Lx = 39.4e-6
Ly = 60.0e-6
X_A, X_B = Rp_1, Lx - Rp_2  # 11.2, 33.8 um

n_steps = 201
p_max = 800.0e6


def D_at(pts, D, x, y):
    dist = np.sqrt((pts[:, 0] - x) ** 2 + (pts[:, 1] - y) ** 2)
    i = np.argmin(dist)
    return float(D[i])


def p_crit_at_D(df, target=0.5):
    below = df[df["D_max"] < target]
    above = df[df["D_max"] >= target]
    if len(below) == 0 or len(above) == 0:
        return None, None, None
    i0, i1 = below.index[-1], above.index[0]
    p0, D0 = df.loc[i0, "pressure_MPa"], df.loc[i0, "D_max"]
    p1, D1 = df.loc[i1, "pressure_MPa"], df.loc[i1, "D_max"]
    return p0 + (target - D0) * (p1 - p0) / (D1 - D0), int(df.loc[i0, "step"]), int(df.loc[i1, "step"])


files = sorted(
    glob.glob(os.path.join(CASE, "output", "simulation_*.vtu")),
    key=lambda f: int(re.search(r"(\d+)", os.path.basename(f)).group(1)),
)

rows = []
for f in files:
    step = int(re.search(r"(\d+)", os.path.basename(f)).group(1))
    pressure = -p_max * step / (n_steps - 1)
    g = pv.read(f)
    pts = g.points
    D = g.point_data["Damage"]
    D_A = D_at(pts, D, X_A, 0.0)
    D_B = D_at(pts, D, X_B, 0.0)
    rows.append({"step": step, "pressure_MPa": pressure / 1e6, "D_max": D.max(), "D_A": D_A, "D_B": D_B})

df = pd.DataFrame(rows).sort_values("step").reset_index(drop=True)
df.to_csv(os.path.join(CASE, "damage_history.csv"), index=False)

p_crit, s0, s1 = p_crit_at_D(df, 0.5)
print(f"p_crit (D_max=0.5) = {p_crit:.3f} MPa (step {s0}->{s1})")

# which tip nucleates first / D_A vs D_B at nucleation
row_nuc_idx = (df["D_max"] - 0.5).abs().idxmin()
step_nuc = int(df.loc[row_nuc_idx, "step"])
D_A_nuc, D_B_nuc = df.loc[row_nuc_idx, "D_A"], df.loc[row_nuc_idx, "D_B"]
print(f"At nucleation frame (step {step_nuc}): D_A={D_A_nuc:.4f}  D_B={D_B_nuc:.4f}")

# D=0.9 frame
idx_09 = (df["D_max"] - 0.9).abs().idxmin()
step_09 = int(df.loc[idx_09, "step"])
D_A_09, D_B_09 = df.loc[idx_09, "D_A"], df.loc[idx_09, "D_B"]
print(f"At D_max~0.9 frame (step {step_09}): D_A={D_A_09:.4f}  D_B={D_B_09:.4f}  D_max_actual={df.loc[idx_09,'D_max']:.4f}")

g09 = pv.read(os.path.join(CASE, "output", f"simulation_{step_09:04d}.vtu"))
pts09 = g09.points
D09 = g09.point_data["Damage"]
mask_lig = (np.abs(pts09[:, 1]) < 1e-7) & (pts09[:, 0] >= Rp_1 - 1e-7) & (pts09[:, 0] <= Lx - Rp_2 + 1e-7)
frac_lig = (D09[mask_lig] > 0.5).mean() if mask_lig.sum() else np.nan
mask_dmg = D09 > 0.5
y_max_dmg = pts09[mask_dmg, 1].max() if mask_dmg.sum() else 0.0
print(f"D=0.9 frame: ligament D>0.5 fraction = {frac_lig*100:.2f}%  "
      f"vertical extent D>0.5 = {y_max_dmg*1e6:.4f} um ({y_max_dmg/Ly*100:.2f}% Ly)")

summary = {
    "p_crit_MPa": p_crit,
    "step_nucleation": step_nuc,
    "D_A_at_nucleation": D_A_nuc,
    "D_B_at_nucleation": D_B_nuc,
    "step_D09": step_09,
    "D_A_at_D09": D_A_09,
    "D_B_at_D09": D_B_09,
    "frac_ligament_D_gt_0.5_at_D09": frac_lig,
    "y_max_D_gt_0.5_um_at_D09": y_max_dmg * 1e6,
    "y_max_over_Ly_at_D09": y_max_dmg / Ly,
}
pd.DataFrame([summary]).to_csv(os.path.join(CASE, "summary.csv"), index=False)
print("\nWrote damage_history.csv and summary.csv")
