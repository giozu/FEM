import glob
import os
import re

import numpy as np
import pandas as pd
import pyvista as pv

CASE = os.path.dirname(os.path.abspath(__file__))
Rp = 11.2e-6
Lx_um = 45
Lx = Lx_um * 1e-6
Ly = 60.0e-6
n_steps = 201
p_max = 800.0e6

cases = {
    "h_cavity=0.30um (coarse)": ("Lx_45_h030", 2306),
    "h_cavity=0.15um (reference)": ("Lx_45", 5510),
    "h_cavity=0.075um (fine)": ("Lx_45_h0075", 15816),
}


def load_history(case_dir):
    files = sorted(
        glob.glob(os.path.join(case_dir, "output", "simulation_*.vtu")),
        key=lambda f: int(re.search(r"(\d+)", os.path.basename(f)).group(1)),
    )
    rows = []
    for f in files:
        step = int(re.search(r"(\d+)", os.path.basename(f)).group(1))
        pressure = -p_max * step / (n_steps - 1)
        g = pv.read(f)
        D = g.point_data["Damage"]
        rows.append({"step": step, "pressure_MPa": pressure / 1e6, "D_max": D.max()})
    return pd.DataFrame(rows).sort_values("step").reset_index(drop=True)


def p_crit_at_D(df, target=0.5):
    below = df[df["D_max"] < target]
    above = df[df["D_max"] >= target]
    if len(below) == 0 or len(above) == 0:
        return None
    i0, i1 = below.index[-1], above.index[0]
    p0, D0 = df.loc[i0, "pressure_MPa"], df.loc[i0, "D_max"]
    p1, D1 = df.loc[i1, "pressure_MPa"], df.loc[i1, "D_max"]
    return p0 + (target - D0) * (p1 - p0) / (D1 - D0)


results = []
for label, (name, n_elem) in cases.items():
    case_dir = os.path.join(CASE, name)
    df = load_history(case_dir)
    df.to_csv(os.path.join(case_dir, "damage_vs_pressure.csv"), index=False)

    p_crit = p_crit_at_D(df, 0.5)

    # nucleation frame: D_max closest to 0.5
    idx_nuc = (df["D_max"] - 0.5).abs().idxmin()
    step_nuc = int(df.loc[idx_nuc, "step"])
    g_nuc = pv.read(os.path.join(case_dir, "output", f"simulation_{step_nuc:04d}.vtu"))
    pts = g_nuc.points
    D_nuc = g_nuc.point_data["Damage"]
    imax = np.argmax(D_nuc)
    x_nuc, y_nuc = pts[imax, 0], pts[imax, 1]

    # D=0.9 frame
    idx_09 = (df["D_max"] - 0.9).abs().idxmin()
    step_09 = int(df.loc[idx_09, "step"])
    g09 = pv.read(os.path.join(case_dir, "output", f"simulation_{step_09:04d}.vtu"))
    pts09 = g09.points
    D09 = g09.point_data["Damage"]

    mask_lig = (np.abs(pts09[:, 1]) < 1e-7) & (pts09[:, 0] >= Rp - 1e-7) & (pts09[:, 0] <= Lx - Rp + 1e-7)
    frac_lig = (D09[mask_lig] > 0.5).mean() if mask_lig.sum() else np.nan

    mask_dmg = D09 > 0.5
    y_max_dmg = pts09[mask_dmg, 1].max() if mask_dmg.sum() else 0.0

    results.append({
        "resolution": label,
        "n_elements": n_elem,
        "p_crit_MPa": p_crit,
        "step_nucleation": step_nuc,
        "D_max_at_nucleation": float(D_nuc[imax]),
        "x_nucleation_um": x_nuc * 1e6,
        "y_nucleation_um": y_nuc * 1e6,
        "step_D09": step_09,
        "D_max_actual_at_D09": float(D09.max()),
        "frac_ligament_D_gt_0.5_at_D09": frac_lig,
        "y_max_D_gt_0.5_um_at_D09": y_max_dmg * 1e6,
        "y_max_over_Ly_at_D09": y_max_dmg / Ly,
    })
    print(f"{label}: p_crit={p_crit:.2f} MPa | nucleation=({x_nuc*1e6:.3f},{y_nuc*1e6:.3f})um "
          f"| D09: lig_frac={frac_lig*100:.1f}% y_max={y_max_dmg*1e6:.3f}um")

out = pd.DataFrame(results)
out.to_csv(os.path.join(CASE, "mesh_convergence_Lx45.csv"), index=False)
print()
pd.set_option("display.width", 200)
print(out.to_string(index=False))
