import glob
import math
import os
import re

import numpy as np
import pandas as pd
import pyvista as pv

SWEEP_DIR = os.path.dirname(os.path.abspath(__file__))
Rp = 11.2e-6
E = 358.0e9
lc = 2.0e-6
sigma_c = 421.5e6
K_IC = 1.0e6  # Pa.m^0.5
REF_168_7 = K_IC / math.sqrt(math.pi * Rp)

Gc_derived = (8.0 / 3.0) * lc * sigma_c**2 / E

CASES = {
    "theta_30": 30.0,
    "theta_40": 40.0,
    "theta_48p1": 48.1,
    "theta_60": 60.0,
    "theta_71": 71.0,
}

n_steps = 201
p_max = 1200.0e6  # Pa magnitude


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
for name, theta_deg in CASES.items():
    case_dir = os.path.join(SWEEP_DIR, name)
    files = sorted(
        glob.glob(os.path.join(case_dir, "output", "simulation_*.vtu")),
        key=lambda f: int(re.search(r"(\d+)", os.path.basename(f)).group(1)),
    )
    if not files:
        print(f"[WARN] no output found for {name}")
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

    ay = Rp * math.tan(math.radians(theta_deg) / 2.0)

    p_crit, s0, s1 = p_crit_at_D(df, 0.5)
    final_Dmax = df["D_max"].iloc[-1]

    if p_crit is None:
        print(f"[FLAG] {name} (theta={theta_deg}): D_max never reaches 0.5 in [0,-{p_max/1e6:.0f} MPa] "
              f"(final D_max = {final_Dmax:.4f}) -- NOT auto-extending, reporting as-is.")
        rows.append({
            "theta_deg": theta_deg, "ay_um": ay * 1e6, "Gc_derived_J_m2": Gc_derived,
            "p_crit_MPa": np.nan, "p_crit_over_168_7": np.nan,
            "status": f"NO_CRACK_final_Dmax={final_Dmax:.4f}",
        })
        continue

    rows.append({
        "theta_deg": theta_deg, "ay_um": ay * 1e6, "Gc_derived_J_m2": Gc_derived,
        "p_crit_MPa": p_crit, "p_crit_over_168_7": abs(p_crit) / (REF_168_7 / 1e6),
        "status": f"OK (step {s0}->{s1})",
    })

out = pd.DataFrame(rows)
out.to_csv(os.path.join(SWEEP_DIR, "theta_sweep_results.csv"), index=False)
pd.set_option("display.width", 160)
print(out.to_string(index=False))
print(f"\n168.7 MPa reference check: K_IC/sqrt(pi*Rp) = {REF_168_7/1e6:.4f} MPa")
print(f"Gc_derived = {Gc_derived:.6f} J/m2 (target 2.645 +/- 0.005)")
