import os

import matplotlib.pyplot as plt
import pandas as pd

CASE = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(CASE, "tip_symmetry_all_meshes.csv"))

styles = {
    "h=0.30um (coarse)": dict(color="tab:blue"),
    "h=0.15um (reference)": dict(color="tab:orange"),
    "h=0.075um (fine)": dict(color="tab:green"),
}

fig, ax = plt.subplots(figsize=(8, 6))
for label, sty in styles.items():
    sub = df[df["mesh"] == label].sort_values("pressure_MPa")
    p = -sub["pressure_MPa"]
    ax.plot(p, sub["D_left"], "-", label=f"{label}, left tip (x=11.2um)", **sty)
    ax.plot(p, sub["D_right"], "--", label=f"{label}, right tip (x=33.8um)", **sty)

ax.set_xlabel("|pressure| (MPa)")
ax.set_ylabel("Damage D at tip")
ax.set_title("Lx=45um: D at both tips vs. pressure, 3 meshes\n"
             "(6 curves -- overlapping if the tips are mechanically equivalent)")
ax.legend(fontsize=8, loc="upper left")
ax.grid(True, alpha=0.4)
ax.set_xlim(380, 460)
ax.set_ylim(-0.02, 1.02)

fig.tight_layout()
out = os.path.join(CASE, "tip_symmetry_check.png")
fig.savefig(out, dpi=200)
print(f"Wrote {out}")
