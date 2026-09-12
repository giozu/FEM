import os

import matplotlib.pyplot as plt
import pandas as pd

SWEEP_DIR = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(SWEEP_DIR, "Lx_sweep_results.csv")).sort_values("Lx_um")

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

ax = axes[0]
ax.plot(df["couverture_lineique"], -df["p_crit_MPa"], "o-", color="tab:blue")
ax.axhline(571.3, color="red", ls="--", label="571.3 MPa (Lx=Ly=60um ref -- NOT domain-matched, see caveat)")
for _, r in df.iterrows():
    ax.annotate(f"Lx={int(r['Lx_um'])}", (r["couverture_lineique"], -r["p_crit_MPa"]),
                textcoords="offset points", xytext=(6, 4), fontsize=8)
ax.set_xlabel(r"Linear coverage $2 R_p / L_x$ (row/strip arrangement, NOT White's $F_c$)")
ax.set_ylabel(r"$|p_{crit}|$ at $D_{max}=0.5$ (MPa)")
ax.set_title("Critical pressure vs. row coverage")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.4)

ax = axes[1]
ax.plot(df["ligament_sur_lc"], -df["p_crit_MPa"], "s-", color="tab:green")
ax.axhline(571.3, color="red", ls="--")
for _, r in df.iterrows():
    ax.annotate(f"Lx={int(r['Lx_um'])}", (r["ligament_sur_lc"], -r["p_crit_MPa"]),
                textcoords="offset points", xytext=(6, 4), fontsize=8)
ax.set_xlabel(r"Ligament / $l_c$")
ax.set_ylabel(r"$|p_{crit}|$ at $D_{max}=0.5$ (MPa)")
ax.set_title("Critical pressure vs. normalized ligament")
ax.grid(True, alpha=0.4)

fig.suptitle("Periodic row of lenticular cavities: p_crit vs. spacing\n"
             "(all points: same Ly=60um, Clamp_x on xmin+xmax; Lx = center-to-center bubble spacing)")
fig.tight_layout()
out = os.path.join(SWEEP_DIR, "Lx_sweep_pcrit.png")
fig.savefig(out, dpi=200)
print(f"Wrote {out}")
