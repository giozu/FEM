import glob
import math
import os
import re

import numpy as np
import pandas as pd
import pyvista as pv

SWEEP_DIR = os.path.dirname(os.path.abspath(__file__))
Rp = 11.2e-6
lc = 2.0e-6

n_steps = 201
p_max = 800.0e6  # Pa magnitude

Lx_values_um = [200, 100, 60, 45, 36, 33, 30]


def p_crit_at_D(df, target=0.5):
    below = df[df["D_max"] < target]
    above = df[df["D_max"] >= target]
    if len(below) == 0 or len(above) == 0:
        return None, None, None
    i0 = below.index[-1]
    i1 = above.index[0]
    p0, D0 = df.loc[i0, "pressure_MPa"], df.loc[i0, "D_max"]
    p1, D1 = df.loc[i1, "pressure_MPa"], df.loc[i1, "D_max"]
    p_c = p0 + (target - D0) * (p1 - p0) / (D1 - D0)
    return p_c, int(df.loc[i0, "step"]), int(df.loc[i1, "step"])


rows = []
for Lx_um in Lx_values_um:
    name = f"Lx_{Lx_um}"
    case_dir = os.path.join(SWEEP_DIR, name)
    Lx = Lx_um * 1e-6
    ligament = Lx - 2 * Rp
    coverage = 2 * Rp / Lx

    files = sorted(
        glob.glob(os.path.join(case_dir, "output", "simulation_*.vtu")),
        key=lambda f: int(re.search(r"(\d+)", os.path.basename(f)).group(1)),
    )
    if not files:
        print(f"[WARN] no output for {name}")
        continue

    steps = []
    for f in files:
        step = int(re.search(r"(\d+)", os.path.basename(f)).group(1))
        pressure = -p_max * step / (n_steps - 1)
        g = pv.read(f)
        D = g.point_data["Damage"]
        steps.append({"step": step, "pressure_MPa": pressure / 1e6, "D_max": D.max()})
    df = pd.DataFrame(steps).sort_values("step").reset_index(drop=True)
    df.to_csv(os.path.join(case_dir, "damage_vs_pressure.csv"), index=False)

    p_crit, s0, s1 = p_crit_at_D(df, 0.5)
    final_Dmax = df["D_max"].iloc[-1]

    row = {
        "Lx_um": Lx_um,
        "ligament_um": ligament * 1e6,
        "ligament_sur_lc": ligament / lc,
        "p_crit_MPa": p_crit if p_crit is not None else np.nan,
        "couverture_lineique": coverage,
        # NOTE: no "p_crit / 571.3" ratio here on purpose -- 571.3 MPa (Lx=Ly=60um,
        # theta sweep) is NOT a domain-matched reference for this Lx-varying sweep
        # (see README: elastic concentration factor itself depends on Lx/Ly aspect
        # ratio). Comparing against it would be misleading.
        "status": f"OK (step {s0}->{s1})" if p_crit is not None else f"NO_CRACK final_Dmax={final_Dmax:.4f}",
    }
    rows.append(row)
    print(f"{name}: ligament={ligament*1e6:.2f}um  p_crit={p_crit if p_crit else 'N/A'}")

out = pd.DataFrame(rows)
out.to_csv(os.path.join(SWEEP_DIR, "Lx_sweep_results.csv"), index=False)
pd.set_option("display.width", 200)
print()
print(out.to_string(index=False))
