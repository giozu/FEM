import glob
import os
import re

import numpy as np
import pandas as pd
import pyvista as pv

CASE = os.path.dirname(os.path.abspath(__file__))
Rp = 11.2e-6

n_steps = 201
p_max = 1200.0e6  # Pa magnitude, matches boundary_conditions.yaml


def cavity_area(vtu_path):
    g = pv.read(vtu_path)
    pts = g.points
    disp = g.point_data["Displacement"]

    surf = g.extract_surface()
    edges = surf.extract_feature_edges(boundary_edges=True, feature_edges=False,
                                        non_manifold_edges=False, manifold_edges=False)
    epts = edges.points
    Lx = Ly = 60.0e-6
    tol = 1e-8
    on_outer = (
        (np.abs(epts[:, 0] - 0.0) < tol) |
        (np.abs(epts[:, 1] - 0.0) < tol) |
        (np.abs(epts[:, 0] - Lx) < tol) |
        (np.abs(epts[:, 1] - Ly) < tol)
    )
    cav_pts_ref = epts[~on_outer]
    if len(cav_pts_ref) < 2:
        return np.nan

    from scipy.spatial import cKDTree
    tree = cKDTree(pts)
    _, idx = tree.query(cav_pts_ref)
    cav_disp = disp[idx]
    cav_pts_def = cav_pts_ref[:, :2] + cav_disp[:, :2]

    order = np.argsort(-cav_pts_def[:, 0])
    arc = cav_pts_def[order]
    poly = np.vstack([[0.0, 0.0], arc, [0.0, 0.0]])
    x, y = poly[:, 0], poly[:, 1]
    quarter_area = 0.5 * np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
    return 4.0 * quarter_area


files = sorted(
    glob.glob(os.path.join(CASE, "output", "simulation_*.vtu")),
    key=lambda f: int(re.search(r"(\d+)", os.path.basename(f)).group(1)),
)

rows = []
for f in files:
    step = int(re.search(r"(\d+)", os.path.basename(f)).group(1))
    pressure = -p_max * step / (n_steps - 1)
    g = pv.read(f)
    D = g.point_data["Damage"]
    area = cavity_area(f)
    rows.append({"step": step, "pressure_MPa": pressure / 1e6, "D_max": D.max(), "cavity_area_m2": area})

df = pd.DataFrame(rows).sort_values("step").reset_index(drop=True)
df.to_csv(os.path.join(CASE, "damage_area_vs_pressure.csv"), index=False)

# p_crit at D_max=0.5, interpolated
below = df[df["D_max"] < 0.5]
above = df[df["D_max"] >= 0.5]
i0, i1 = below.index[-1], above.index[0]
p0, D0 = df.loc[i0, "pressure_MPa"], df.loc[i0, "D_max"]
p1, D1 = df.loc[i1, "pressure_MPa"], df.loc[i1, "D_max"]
p_crit = p0 + (0.5 - D0) * (p1 - p0) / (D1 - D0)
print(f"p_crit (D_max=0.5) = {p_crit:.3f} MPa (between step {df.loc[i0,'step']} "
      f"[{p0:.2f} MPa, D={D0:.4f}] and step {df.loc[i1,'step']} [{p1:.2f} MPa, D={D1:.4f}])")

idx_closest = (df["D_max"] - 0.5).abs().idxmin()
step_closest = int(df.loc[idx_closest, "step"])
print(f"Frame closest to D_max=0.5: step {step_closest} "
      f"(D_max={df.loc[idx_closest,'D_max']:.4f}, p={df.loc[idx_closest,'pressure_MPa']:.2f} MPa)")

with open(os.path.join(CASE, "pcrit_area_summary.txt"), "w") as f:
    f.write(f"p_crit_MPa={p_crit}\nstep_closest_to_D05={step_closest}\n")

print(df.to_string(index=False))
