import os

import matplotlib.pyplot as plt
import pandas as pd

CASE = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(CASE, "damage_history.csv")).sort_values("step")

fig, ax = plt.subplots(figsize=(8, 6))
p = -df["pressure_MPa"]
ax.plot(p, df["D_A"], "-", color="tab:red", lw=2, label="D at tip A (Rp_1=11.2um, x=11.2um)")
ax.plot(p, df["D_B"], "-", color="tab:blue", lw=2, label="D at tip B (Rp_2=5.6um, x=33.8um)")
ax.axhline(0.5, color="gray", ls=":", lw=1)
ax.axvline(476.0, color="tab:red", ls="--", lw=1, alpha=0.6, label="p_crit (D_A=0.5) = 476.0 MPa")
ax.axvline(422.5, color="black", ls="--", lw=1, alpha=0.6,
           label="monodisperse Lx_45 reference, p_crit = 422.5 MPa")

ax.set_xlabel("|pressure| (MPa)")
ax.set_ylabel("Damage D at tip")
ax.set_title("Polydisperse cell (Lx=39.4um, ligament=22.6um): which tip nucleates?\n"
             "D_A vs D_B vs applied pressure")
ax.legend(fontsize=9, loc="upper left")
ax.grid(True, alpha=0.4)
ax.set_ylim(-0.02, 1.02)

fig.tight_layout()
out = os.path.join(CASE, "tip_damage_A_vs_B.png")
fig.savefig(out, dpi=200)
print(f"Wrote {out}")
