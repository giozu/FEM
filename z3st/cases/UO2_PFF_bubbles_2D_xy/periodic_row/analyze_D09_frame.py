import os

import numpy as np
import pandas as pd
import pyvista as pv

CASE = os.path.dirname(os.path.abspath(__file__))
Rp = 11.2e-6
Ly = 60.0e-6

d09_steps = {200: 156, 100: 142, 60: 123, 45: 108, 36: 93, 33: 86, 30: 77}

rows = []
for Lx_um, step in d09_steps.items():
    Lx = Lx_um * 1e-6
    g = pv.read(os.path.join(CASE, f"Lx_{Lx_um}", "output", f"simulation_{step:04d}.vtu"))
    pts = g.points
    D = g.point_data["Damage"]

    # ligament fraction with D>0.5 (along y=0, x in [Rp, Lx-Rp])
    mask_lig = (np.abs(pts[:, 1]) < 1e-7) & (pts[:, 0] >= Rp - 1e-7) & (pts[:, 0] <= Lx - Rp + 1e-7)
    frac_lig_gt05 = (D[mask_lig] > 0.5).mean() if mask_lig.sum() else np.nan

    # max vertical extent of the D>0.5 region (over the WHOLE domain, not just ligament)
    mask_dmg = D > 0.5
    y_max_dmg = pts[mask_dmg, 1].max() if mask_dmg.sum() else 0.0

    rows.append({
        "Lx_um": Lx_um, "step_D09": step,
        "D_max_actual": D.max(),
        "frac_ligament_D_gt_0.5": frac_lig_gt05,
        "n_ligament_pts": int(mask_lig.sum()),
        "y_max_D_gt_0.5_um": y_max_dmg * 1e6,
        "y_max_over_Ly": y_max_dmg / Ly,
    })
    print(f"Lx={Lx_um}um: frac_lig(D>0.5)={frac_lig_gt05:.3f}  y_max(D>0.5)={y_max_dmg*1e6:.2f}um "
          f"({y_max_dmg/Ly*100:.1f}% of Ly)")

out = pd.DataFrame(rows)
out.to_csv(os.path.join(CASE, "D09_frame_analysis.csv"), index=False)
print()
print(out.to_string(index=False))
