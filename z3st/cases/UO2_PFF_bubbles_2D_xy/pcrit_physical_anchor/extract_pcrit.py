import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import pyvista as pv

case_dir = sys.argv[1]
n_steps = int(sys.argv[2])
p_max = float(sys.argv[3])  # Pa magnitude

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

df = pd.DataFrame(rows).sort_values("step").reset_index(drop=True)
df.to_csv(os.path.join(case_dir, "damage_vs_pressure.csv"), index=False)

below = df[df["D_max"] < 0.5]
above = df[df["D_max"] >= 0.5]
name = os.path.basename(case_dir)
if len(below) == 0 or len(above) == 0:
    print(f"{name}: D_max never crosses 0.5 in this range (max observed {df['D_max'].max():.4f})")
else:
    i0, i1 = below.index[-1], above.index[0]
    p0, D0 = df.loc[i0, "pressure_MPa"], df.loc[i0, "D_max"]
    p1, D1 = df.loc[i1, "pressure_MPa"], df.loc[i1, "D_max"]
    p_crit = p0 + (0.5 - D0) * (p1 - p0) / (D1 - D0)
    print(f"{name}: p_crit (D_max=0.5) = {p_crit:.3f} MPa "
          f"(step {df.loc[i0,'step']} [{p0:.2f} MPa, D={D0:.4f}] -> "
          f"step {df.loc[i1,'step']} [{p1:.2f} MPa, D={D1:.4f}])")
