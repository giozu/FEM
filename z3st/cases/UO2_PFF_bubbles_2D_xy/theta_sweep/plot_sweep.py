import math
import os

import matplotlib.pyplot as plt
import pandas as pd

SWEEP_DIR = os.path.dirname(os.path.abspath(__file__))
Rp = 11.2e-6

df = pd.read_csv(os.path.join(SWEEP_DIR, "theta_sweep_results.csv")).sort_values("theta_deg")
# NOTE (correction): h_cavity was scaled per theta assuming rho=ay^2/Rp is
# a tip curvature radius to resolve. That premise is wrong -- the cavity
# tip is a wedge corner (two arcs meeting at 2*theta), not an ellipse; a
# corner has no curvature radius. The real (size-independent) constraint
# is h <= lc/4 = 0.5e-6, satisfied throughout. "rho_over_h" below is kept
# only as a resolution-doubling label for the N=15 vs N=30 mesh-refinement
# check (which remains a valid convergence test on its own), not as a
# curvature-resolution ratio.
df["rho_um"] = (df["ay_um"] * 1e-6) ** 2 / Rp * 1e6
df["rho_over_h"] = 15.0

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
ax.plot(df["theta_deg"], -df["p_crit_MPa"], "o-", color="tab:blue")
ax.axvline(48.1, color="gray", linestyle=":", linewidth=1)
ax.axhline(576, color="tab:red", linestyle="--", linewidth=1, label="576 MPa (rescaled anchor @ 48.1°)")
ax.set_xlabel(r"Semi-dihedral angle $\theta$ (deg)")
ax.set_ylabel(r"$|p_{crit}|$ at $D_{max}=0.5$ (MPa)")
ax.set_title("Critical pressure vs. dihedral angle")
ax.grid(True, alpha=0.4)
ax.legend()

ax = axes[1]
ax.plot(df["theta_deg"], df["p_crit_over_168_7"], "s-", color="tab:green")
ax.axhline(3.42, color="tab:red", linestyle="--", linewidth=1, label="3.42 (target @ 48.1°)")
ax.axvline(48.1, color="gray", linestyle=":", linewidth=1)
ax.set_xlabel(r"Semi-dihedral angle $\theta$ (deg)")
ax.set_ylabel(r"Blunting factor $p_{crit}/168.7$ MPa")
ax.set_title("Blunting factor vs. dihedral angle")
ax.grid(True, alpha=0.4)
ax.legend()

fig.tight_layout()
out1 = os.path.join(SWEEP_DIR, "theta_sweep_pcrit.png")
fig.savefig(out1, dpi=200)
print(f"Wrote {out1}")

# second figure: mesh-convergence check at the 48.1 deg anchor point
# (N = rho/h_cavity = 15, the value used for the whole sweep, vs N = 30)
fig2, ax2 = plt.subplots(figsize=(6, 5))
N_vals = [15, 30]
p_vals = [571.303464, 570.381]
ax2.plot(N_vals, p_vals, "o-", color="tab:purple")
for n, p in zip(N_vals, p_vals):
    ax2.annotate(f"{p:.1f} MPa", (n, p), textcoords="offset points", xytext=(8, 6))
ax2.set_xlabel(r"Tip-curvature mesh resolution $N = \rho / h_{cavity}$")
ax2.set_ylabel(r"$|p_{crit}|$ at $D_{max}=0.5$, $\theta=48.1$deg (MPa)")
ax2.set_title("Mesh convergence check at the 48.1 deg anchor\n(0.16% change from N=15 to N=30)")
ax2.set_xlim(10, 35)
ax2.grid(True, alpha=0.4)
fig2.tight_layout()
out2 = os.path.join(SWEEP_DIR, "theta_sweep_convergence_check.png")
fig2.savefig(out2, dpi=200)
print(f"Wrote {out2}")

df.to_csv(os.path.join(SWEEP_DIR, "theta_sweep_results.csv"), index=False)
print(df.to_string(index=False))
