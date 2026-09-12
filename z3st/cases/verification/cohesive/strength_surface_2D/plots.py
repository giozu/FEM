#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Figures for the cohesive strength-surface cases.

Renders ``output/strength_surface.png`` from ``output/response.csv``. The same
script serves the r = 1, r = 2 and r = inf cases: it reads which one it is from
the case's own input.yaml and material.yaml.

Three panels:

  (a) the pressure-shear plane, where the strength surface is written. The
      analytic dS_0 is Eq. (110) with its p < 0 branch from Eq. (109); the
      loading path is a straight ray because the state stays homogeneous and
      elastic until it reaches the surface.
  (b) the same in the plane of the in-plane normal stresses, which is how the
      paper plots it (Fig. 22). The analytic locus here is obtained by sweeping
      the loading angle and, for each, scaling the elastic ray until it meets
      dS_0 -- no extra simulation is needed for it.
  (c) the measurement itself: the phase field lifting off zero as the ray
      crosses the surface, with the bracket the load stepping resolves it to.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import yaml

from z3st.utils.plotstyle import apply as _apply_plotstyle

_apply_plotstyle()

CASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(CASE_DIR, "output")
CSV = os.path.join(OUT, "response.csv")
PNG = os.path.join(OUT, "strength_surface.png")

with open(os.path.join(CASE_DIR, "input.yaml")) as f:
    r = str(yaml.safe_load(f)["models"]["cohesive"]["r_norm"])
with open(os.path.join(CASE_DIR, "material.yaml")) as f:
    mat = yaml.safe_load(f)
E, nu = float(mat["E"]), float(mat["nu"])
p_c, tau_c = float(mat["p_c"]), float(mat["tau_c"])
kappa, mu = E / (3 * (1 - 2 * nu)), E / (2 * (1 + nu))

LABEL = {"1": "rectangle", "2": "ellipse", "inf": "Drucker-Prager"}[r]


def surface_tau(p):
    """tau on dS_0 for a given p. Eq. (110) for p >= 0, Eq. (109) below it."""
    p = np.atleast_1d(np.asarray(p, dtype=float))
    out = np.full_like(p, tau_c)               # p < 0: tau = tau_c for every r
    pos = p >= 0
    x = np.clip(p[pos] / p_c, 0.0, 1.0)
    if r == "1":                               # rectangle: tau = tau_c up to p_c
        out[pos] = np.where(x <= 1.0, tau_c, np.nan)
    else:
        e = float(r) / (float(r) - 1.0) if r != "inf" else 1.0
        out[pos] = tau_c * (1.0 - x ** e) ** (1.0 / e)
    return out


def ray_scale(p_dir, tau_dir):
    """Factor that puts the ray (p_dir, tau_dir) onto dS_0."""
    if r == "1":
        return 1.0 / max(p_dir / p_c, tau_dir / tau_c)
    if r == "inf":
        return 1.0 / (p_dir / p_c + tau_dir / tau_c)
    e = float(r) / (float(r) - 1.0)
    return ((p_dir / p_c) ** e + (tau_dir / tau_c) ** e) ** (-1.0 / e)


data = np.genfromtxt(CSV, delimiter=",", names=True)
p, tau = np.atleast_1d(data["p_Pa"]), np.atleast_1d(data["tau_Pa"])
sxx, syy = np.atleast_1d(data["sigma_xx_Pa"]), np.atleast_1d(data["sigma_yy_Pa"])
alpha = np.atleast_1d(data["alpha_max"])

i = int(np.flatnonzero(alpha > 1e-6)[0])
M = 1e6                                        # Pa -> MPa for the axes

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

# --.. ..- .-.. .-.. --- (a) pressure-shear plane --.. ..- .-.. .-.. ---
ax = axes[0]
pp = np.linspace(-1.2 * p_c, p_c, 800)
ax.plot(pp / M, surface_tau(pp) / M, "-", marker="", lw=1.8, color="0.35",
        label=rf"$\partial S_0$, Eq. (110) [{LABEL}]")
if r == "1":                                   # close the corner of the rectangle
    ax.plot([p_c / M, p_c / M], [0, tau_c / M], "-", marker="", lw=1.8, color="0.35")
ax.plot(p / M, tau / M, "-", marker="", lw=1.6, color="#0072B2", label="loading path")
ax.plot([p[i - 1] / M, p[i] / M], [tau[i - 1] / M, tau[i] / M],
        marker="o", ms=5, ls="none", color="#D55E00", label="nucleation bracket")
ax.set_xlabel(r"$p$ (MPa)")
ax.set_ylabel(r"$\tau$ (MPa)")
ax.set_xlim(-1.2 * p_c / M, 1.25 * p_c / M)
ax.set_ylim(0, 1.25 * tau_c / M)
ax.legend(loc="lower left", fontsize=9)
ax.set_title("(a) pressure-shear plane")
ax.grid(alpha=0.3)

# --.. ..- .-.. .-.. --- (b) in-plane normal stresses --.. ..- .-.. .-.. ---
# Sweep the loading angle; for each, the elastic ray is scaled onto dS_0. The
# stress follows sigma = p I + 2 mu dev(eps), with eps_zz = 0 (plane strain).
ax = axes[1]
th = np.linspace(0, 2 * np.pi, 721)
exx, eyy = np.cos(th), np.sin(th)
tr = exx + eyy
p_dir = kappa * tr
dev_norm = np.sqrt(2 / 3 * (exx**2 + eyy**2 - exx * eyy))
tau_dir = 2 * mu * dev_norm
Sxx, Syy = [], []
for k in range(th.size):
    if p_dir[k] >= 0:
        s = ray_scale(p_dir[k], tau_dir[k])
    else:                                      # p < 0 branch: tau = tau_c
        s = tau_c / tau_dir[k]
    a = s
    P = p_dir[k] * a
    Sxx.append(P + 2 * mu * a * (exx[k] - tr[k] / 3))
    Syy.append(P + 2 * mu * a * (eyy[k] - tr[k] / 3))
ax.plot(np.array(Sxx) / M, np.array(Syy) / M, "-", marker="", lw=1.8,
        color="0.35", label=r"$\partial S_0$ (analytic)")
ax.plot(sxx / M, syy / M, "-", marker="", lw=1.6, color="#0072B2",
        label="loading path")
ax.plot([sxx[i - 1] / M, sxx[i] / M], [syy[i - 1] / M, syy[i] / M],
        marker="o", ms=5, ls="none", color="#D55E00", label="nucleation bracket")
ax.axhline(0, lw=0.6, color="0.7", marker="")
ax.axvline(0, lw=0.6, color="0.7", marker="")
ax.set_xlabel(r"$\sigma_{xx}$ (MPa)")
ax.set_ylabel(r"$\sigma_{yy}$ (MPa)")
ax.legend(loc="lower right", fontsize=9)
ax.set_title("(b) in-plane normal stresses")
ax.grid(alpha=0.3)

# --.. ..- .-.. .-.. --- (c) the measurement --.. ..- .-.. .-.. ---
ax = axes[2]
# Normalise the distance along the ray by where the ray meets dS_0, so that
# 1.0 is exactly the analytic onset and the curve can be read as an error.
k = ray_scale(p[-1], tau[-1])
ratio = np.hypot(p, tau) / np.hypot(p[-1] * k, tau[-1] * k)
ax.plot(ratio, alpha, "-", marker="", lw=1.8, color="#0072B2", label=r"$\max\alpha$")
ax.axvline(1.0, ls="--", lw=1.4, color="0.35", marker="",
           label=r"$\partial S_0$ (analytic)")
ax.axvspan(ratio[i - 1], ratio[i], color="#D55E00", alpha=0.25, label="bracket")
ax.set_xlabel(r"stress along the ray, normalised by $\partial S_0$")
ax.set_ylabel(r"$\max_\Omega \alpha$")
ax.set_xlim(0.95, min(ratio.max(), 1.06))
ax.set_ylim(bottom=0)
ax.legend(loc="upper left", fontsize=9)
ax.set_title("(c) onset of damage")
ax.grid(alpha=0.3)

fig.tight_layout()
fig.savefig(PNG)
print(f"[INFO] wrote {PNG}")
