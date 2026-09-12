#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Figures for verification/cohesive/bar_1D.

Renders ``output/cohesive_bar_1D.png`` from ``output/response.csv`` (written by
diagnostics.py) plus, when the run's VTU files are still present, the phase-field
profile. Four panels, each against its closed form in Vicentini et al. (2026):

  (a) structural response, Eq. (85) -- including the snap-back branch that
      displacement control cannot follow, which is why the computed curve drops
      vertically instead;
  (b) cohesive law, Eq. (78), the traction-separation relation the model
      delivers;
  (c) energies, with the fracture energy tending to Gc;
  (d) phase-field profile, Eq. (75), which is what fixes the regularisation
      length.

The dashed reference in (a)-(c) uses the fracture energy the discretisation
actually delivers, Gc (1 + h/2l) by Eq. (149); the continuum Gc is shown too, so
the gap between them is the discretisation bias rather than a modelling error.
"""

import os
import re
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
import yaml

from z3st.utils.plotstyle import apply as _apply_plotstyle

_apply_plotstyle()

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(CASE_DIR, "output")
CSV = os.path.join(OUT, "response.csv")
PNG = os.path.join(OUT, "cohesive_bar_1D.png")

# --.. ..- .-.. .-.. --- case parameters --.. ..- .-.. .-.. ---
with open(os.path.join(CASE_DIR, "geometry.yaml")) as f:
    L2 = float(yaml.safe_load(f)["Lx"])            # 2 L
with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    ell = float(yaml.safe_load(f)["models"]["cohesive"]["ell"])
with open(os.path.join(CASE_DIR, "material.yaml")) as f:
    mat = yaml.safe_load(f)
E, Gc, sigma_c = float(mat["E"]), float(mat["Gc"]), float(mat["p_c"])

with open(os.path.join(CASE_DIR, "mesh.geo")) as f:
    h = L2 / (int(re.search(r"nx\s*=\s*(\d+)", f.read()).group(1)) - 1)

U_peak = L2 * sigma_c / E
Gc_h = Gc * (1.0 + h / (2.0 * ell))               # Eq. (149), AT2
B = L2 * sigma_c**2 / (Gc * E)                    # continuum brittleness ratio
B_h = L2 * sigma_c**2 / (Gc_h * E)

data = np.genfromtxt(CSV, delimiter=",", names=True)
U = np.atleast_1d(data["U_t_m"]) / U_peak
sigma = np.atleast_1d(data["sigma_xx_Pa"]) / sigma_c
alpha = np.atleast_1d(data["alpha_max"])
E_el = np.atleast_1d(data["E_el_J"])
E_frac = np.atleast_1d(data["E_frac_J"])
cracked = alpha > 1e-6

fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))

# --.. ..- .-.. .-.. --- (a) structural response, Eq. (85) --.. ..- .-.. ---
ax = axes[0, 0]
a0 = np.linspace(0.0, 0.999, 2000)
for ratio, label, style in ((B, rf"Eq. (85), $G_c$", "-"),
                            (B_h, rf"Eq. (85), $G_c(1{{+}}h/2\ell)$", "--")):
    ax.plot(a0 / ((1 - a0) * ratio) + (1 - a0) ** 2, (1 - a0) ** 2,
            style, marker="", color="0.35" if style == "-" else "#D55E00", lw=1.5,
            label=label)
ax.plot(U, sigma, "-", marker="", lw=2.0, color="#0072B2", label="Z3ST")
ax.plot(U[cracked], sigma[cracked], marker="o", ms=3.5, ls="none",
        color="#0072B2", label="localized states")
ax.set_xlabel(r"$U_t\,E\,/\,(2L\,\sigma_c)$")
ax.set_ylabel(r"$\sigma\,/\,\sigma_c$")
ax.set_xlim(0, 5.2)
ax.set_ylim(0, 1.1)
ax.legend(loc="upper right", framealpha=0.95)

# The softening branch lives in the bottom few percent of the axis, where the
# agreement that the case actually verifies is invisible at full scale.
inset = ax.inset_axes([0.42, 0.20, 0.54, 0.44])
inset.plot(a0 / ((1 - a0) * B) + (1 - a0) ** 2, (1 - a0) ** 2,
           "-", marker="", color="0.35", lw=1.2)
inset.plot(a0 / ((1 - a0) * B_h) + (1 - a0) ** 2, (1 - a0) ** 2,
           "--", marker="", color="#D55E00", lw=1.2)
inset.plot(U[cracked], sigma[cracked], marker="o", ms=3.5, ls="none",
           color="#0072B2")
inset.set_xlim(0.9, 3.0)
inset.set_ylim(0, 0.06)
inset.tick_params(labelsize=9)
inset.grid(alpha=0.3)
inset.set_xlabel(r"softening branch, $B=%.2f$ ($G_c$) / $%.2f$ (discrete)"
                 % (B, B_h), fontsize=8, labelpad=1)
ax.set_title("(a) structural response")
ax.grid(alpha=0.3)

# --.. ..- .-.. .-.. --- (b) cohesive law, Eq. (78) --.. ..- .-.. .-.. ---
ax = axes[0, 1]
if "crack_opening_m" in data.dtype.names:
    delta = np.atleast_1d(data["crack_opening_m"])[cracked]
    d_ref = np.linspace(0.0, max(delta.max() * 1.05, 1e-12), 400)
    ax.plot(d_ref * sigma_c / Gc, 1.0 / (1.0 + sigma_c * d_ref / Gc) ** 2,
            "-", marker="", lw=1.5, color="0.35", label=r"Eq. (78), $G_c$")
    ax.plot(d_ref * sigma_c / Gc, 1.0 / (1.0 + sigma_c * d_ref / Gc_h) ** 2,
            "--", marker="", lw=1.5, color="#D55E00",
            label=r"Eq. (78), $G_c(1+h/2\ell)$")
    ax.plot(delta * sigma_c / Gc, sigma[cracked], marker="o", ms=4,
            ls="none", color="#0072B2", label="Z3ST")
    ax.set_xlabel(r"$[\![u]\!]\,\sigma_c\,/\,G_c$")
    ax.set_ylabel(r"$\hat{\sigma}\,/\,\sigma_c$")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)
else:
    ax.text(0.5, 0.5, "re-run the case:\nresponse.csv predates crack_opening_m",
            ha="center", va="center", transform=ax.transAxes)
    ax.set_axis_off()
ax.set_title("(b) cohesive law")

# --.. ..- .-.. .-.. --- (c) energies --.. ..- .-.. .-.. ---
ax = axes[1, 0]
ax.plot(U, E_el, "-", marker="", lw=1.8, label=r"$\Psi$ (elastic)")
ax.plot(U, E_frac, "-", marker="", lw=1.8, label=r"$\mathcal{D}_\ell$ (fracture)")
ax.axhline(Gc, ls="-", lw=1.2, color="0.35", label=r"$G_c$")
ax.axhline(Gc_h, ls="--", lw=1.2, color="#D55E00", label=r"$G_c(1+h/2\ell)$")
ax.set_xlabel(r"$U_t\,E\,/\,(2L\,\sigma_c)$")
ax.set_ylabel("energy (J)")
ax.set_xlim(0, 5.2)
ax.set_yscale("log")
ax.legend(loc="lower right")
ax.set_title("(c) energy balance")
ax.grid(alpha=0.3)

# --.. ..- .-.. .-.. --- (d) phase-field profile, Eq. (75) --.. ..- .-.. ---
ax = axes[1, 1]
vtus = sorted(glob(os.path.join(OUT, "fields_*.vtu")))
if vtus:
    import pyvista as pv

    mesh = pv.read(vtus[-1])
    x = np.asarray(mesh.points)[:, 0]
    d = np.asarray(mesh.point_data["Damage"])
    order = np.argsort(x)
    x, d = x[order], d[order]
    x0 = x[np.argmax(d)]
    ax.plot((x - x0) / ell, d, "-", marker="", lw=1.8, color="#0072B2",
            label="Z3ST")
    xi = np.linspace((x - x0).min(), (x - x0).max(), 800)
    ax.plot(xi / ell, d.max() * np.exp(-np.abs(xi) / ell), "--", marker="",
            lw=1.5, color="0.35",
            label=r"Eq. (75), $\alpha_0 e^{-|x|/\ell}$")
    ax.set_xlim(-8, 8)
    ax.set_xlabel(r"$(x - x_{\rm crack})\,/\,\ell$")
    ax.set_ylabel(r"$\alpha$")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)
else:
    ax.text(0.5, 0.5, "no VTU output found\n(run ./Allrun to populate)",
            ha="center", va="center", transform=ax.transAxes)
    ax.set_axis_off()
ax.set_title(r"(d) phase-field profile at the last step")

fig.tight_layout()
fig.savefig(PNG)
print(f"[INFO] wrote {PNG}")
