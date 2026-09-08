import os

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

CASE = os.path.dirname(os.path.abspath(__file__))
Rp = 11.2e-6

Lx_values = [200, 100, 60, 45, 36, 33, 30]

nucleation_info = {
    200: (154, -616.0, 0.4875), 100: (140, -560.0, 0.5752), 60: (121, -484.0, 0.5597),
    45: (106, -424.0, 0.5592), 36: (91, -364.0, 0.5127), 33: (84, -336.0, 0.4923), 30: (75, -300.0, 0.4409),
}
ligament_um = {200: 177.6, 100: 77.6, 60: 37.6, 45: 22.6, 36: 13.6, 33: 10.6, 30: 7.6}

d09 = pd.read_csv(os.path.join(CASE, "D09_frame_analysis.csv")).set_index("Lx_um")

fig, axes = plt.subplots(7, 3, figsize=(15.5, 21))

col_titles = ["Nucleation (D_max~0.5)", "D_max~0.9 (propagation)", "Final frame (p=-800MPa, reference)"]

for row, Lx_um in enumerate(Lx_values):
    step_n, p_n, D_n = nucleation_info[Lx_um]
    lig = ligament_um[Lx_um]
    r = d09.loc[Lx_um]

    for col, tag in enumerate(["nucleation", "d09", "final"]):
        ax = axes[row, col]
        img = mpimg.imread(os.path.join(CASE, f"_panel_Lx{Lx_um}_{tag}.png"))
        ax.imshow(img)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_edgecolor("black")
            spine.set_linewidth(0.6)

        if row == 0:
            ax.set_title(col_titles[col], fontsize=11)

        if col == 0:
            ax.text(0.03, 0.94, f"step {step_n}\np={p_n:.0f} MPa\nD_max={D_n:.2f}",
                    transform=ax.transAxes, fontsize=7.5, va="top", color="white",
                    bbox=dict(facecolor="black", alpha=0.4, pad=1.5))
        elif col == 1:
            ax.text(0.03, 0.94,
                    f"step {int(r['step_D09'])}\nD_max={r['D_max_actual']:.2f}\n"
                    f"lig. D>0.5: {r['frac_ligament_D_gt_0.5']*100:.1f}%\n"
                    f"y_max(D>0.5): {r['y_max_D_gt_0.5_um']:.2f} um ({r['y_max_over_Ly']*100:.1f}% Ly)",
                    transform=ax.transAxes, fontsize=7.5, va="top", color="white",
                    bbox=dict(facecolor="black", alpha=0.4, pad=1.5))

    axes[row, 0].text(-0.30, 0.5, f"Lx={Lx_um} um\nligament={lig} um\n({lig/2.0:.1f} lc)",
                       fontsize=10, ha="center", va="center", transform=axes[row, 0].transAxes)

fig.suptitle("Periodic row: propagation morphology (column 2, D_max~0.9) vs. saturated state (column 3)\n"
             "same framing and 0-1 (Damage) scale for every thumbnail -- Lx = center-to-center spacing",
             fontsize=13, y=0.995)
fig.subplots_adjust(left=0.14, right=0.98, top=0.955, bottom=0.035, hspace=0.30, wspace=0.05)

sm = ScalarMappable(cmap="turbo", norm=Normalize(vmin=0, vmax=1))
cbar_ax = fig.add_axes([0.30, 0.01, 0.4, 0.006])
fig.colorbar(sm, cax=cbar_ax, orientation="horizontal", label="Damage D")

out = os.path.join(CASE, "Lx_sweep_damage_panel_3col.png")
fig.savefig(out, dpi=165)
print(f"Wrote {out}")
