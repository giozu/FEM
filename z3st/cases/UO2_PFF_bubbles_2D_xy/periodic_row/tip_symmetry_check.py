import glob
import os
import re

import numpy as np
import pandas as pd
import pyvista as pv

CASE = os.path.dirname(os.path.abspath(__file__))
Rp = 11.2e-6
Lx = 45.0e-6
X_LEFT, X_RIGHT = Rp, Lx - Rp  # 11.2, 33.8 um

n_steps = 201
p_max = 800.0e6

meshes = {
    "h=0.30um (coarse)": "Lx_45_h030",
    "h=0.15um (reference)": "Lx_45",
    "h=0.075um (fine)": "Lx_45_h0075",
}

nucleation_steps = {"h=0.30um (coarse)": 106, "h=0.15um (reference)": 106, "h=0.075um (fine)": 105}


def D_at(pts, D, x, y):
    dist = np.sqrt((pts[:, 0] - x) ** 2 + (pts[:, 1] - y) ** 2)
    i = np.argmin(dist)
    return float(D[i])


all_curves = {}
for label, name in meshes.items():
    case_dir = os.path.join(CASE, name)
    files = sorted(
        glob.glob(os.path.join(case_dir, "output", "simulation_*.vtu")),
        key=lambda f: int(re.search(r"(\d+)", os.path.basename(f)).group(1)),
    )
    rows = []
    for f in files:
        step = int(re.search(r"(\d+)", os.path.basename(f)).group(1))
        pressure = -p_max * step / (n_steps - 1)
        g = pv.read(f)
        pts = g.points
        D = g.point_data["Damage"]
        D_left = D_at(pts, D, X_LEFT, 0.0)
        D_right = D_at(pts, D, X_RIGHT, 0.0)
        rows.append({"step": step, "pressure_MPa": pressure / 1e6, "D_left": D_left, "D_right": D_right})
    df = pd.DataFrame(rows).sort_values("step").reset_index(drop=True)
    df.to_csv(os.path.join(case_dir, "tip_symmetry.csv"), index=False)
    all_curves[label] = df
    print(f"{label}: wrote {len(df)} rows")

print()
print("=== 1. D at both tips, at the nucleation frame ===")
for label in meshes:
    df = all_curves[label]
    step_nuc = nucleation_steps[label]
    row = df[df["step"] == step_nuc].iloc[0]
    Dl, Dr = row["D_left"], row["D_right"]
    rel_diff = abs(Dl - Dr) / max(Dl, Dr) * 100 if max(Dl, Dr) > 0 else 0.0
    print(f"{label}: step {step_nuc}, p={row['pressure_MPa']:.1f} MPa -> "
          f"D_left={Dl:.4f}  D_right={Dr:.4f}  rel_diff={rel_diff:.2f}%")

out_all = pd.concat(
    [df.assign(mesh=label) for label, df in all_curves.items()], ignore_index=True
)
out_all.to_csv(os.path.join(CASE, "tip_symmetry_all_meshes.csv"), index=False)
