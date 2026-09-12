#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Comparison figure across the three strength-surface cases.

Reads what the r = 1, r = 2 and r = inf cases already computed and draws them
together, which is the claim the paper is built around: the same model, the
same code path, three deliberately different strength surfaces.

Run it after the three cases (each ``./Allrun``); it starts no simulation of
its own. Writes ``strength_surfaces.png`` next to this script.

Left panel normalises by p_c and tau_c, so the three surfaces share a unit box
and differ in shape alone. Right panel keeps physical units, where they also
differ in size: the paper picks p_c = 8.9 MPa for r = 1 against 12.4 MPa for
the other two, so that every variant sits at the same ell / ell_ch = 1/4.
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import yaml

from z3st.utils.plotstyle import apply as _apply_plotstyle

_apply_plotstyle()

HERE = os.path.dirname(os.path.abspath(__file__))
PNG = os.path.join(HERE, "strength_surfaces.png")

CASES = [
    ("strength_surface_r1_2D", "1", "rectangle", "#0072B2"),
    ("strength_surface_2D", "2", "ellipse", "#D55E00"),
    ("strength_surface_rinf_2D", "inf", "Drucker-Prager", "#009E73"),
]
M = 1e6


def load(case):
    with open(os.path.join(HERE, case, "material.yaml")) as f:
        mat = yaml.safe_load(f)
    csv = os.path.join(HERE, case, "output", "response.csv")
    if not os.path.exists(csv):
        raise SystemExit(f"{case}: no output/response.csv -- run ./Allrun there first.")
    d = np.genfromtxt(csv, delimiter=",", names=True)
    return mat, d


def surface_tau(p, r, p_c, tau_c):
    """tau on dS_0: Eq. (110) for p >= 0, the tau = tau_c branch below it."""
    p = np.atleast_1d(np.asarray(p, dtype=float))
    out = np.full_like(p, tau_c)
    pos = p >= 0
    x = np.clip(p[pos] / p_c, 0.0, 1.0)
    if r == "1":
        out[pos] = np.where(x <= 1.0, tau_c, np.nan)
    else:
        e = float(r) / (float(r) - 1.0) if r != "inf" else 1.0
        out[pos] = tau_c * (1.0 - x ** e) ** (1.0 / e)
    return out


def ray_scale(p_dir, tau_dir, r, p_c, tau_c):
    if r == "1":
        return 1.0 / max(p_dir / p_c, tau_dir / tau_c)
    if r == "inf":
        return 1.0 / (p_dir / p_c + tau_dir / tau_c)
    e = float(r) / (float(r) - 1.0)
    return ((p_dir / p_c) ** e + (tau_dir / tau_c) ** e) ** (-1.0 / e)


fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2))

for case, r, label, colour in CASES:
    mat, d = load(case)
    E, nu = float(mat["E"]), float(mat["nu"])
    p_c, tau_c = float(mat["p_c"]), float(mat["tau_c"])
    kappa, mu = E / (3 * (1 - 2 * nu)), E / (2 * (1 + nu))

    p, tau = np.atleast_1d(d["p_Pa"]), np.atleast_1d(d["tau_Pa"])
    sxx, syy = np.atleast_1d(d["sigma_xx_Pa"]), np.atleast_1d(d["sigma_yy_Pa"])
    i = int(np.flatnonzero(np.atleast_1d(d["alpha_max"]) > 1e-6)[0])
    p_m, tau_m = 0.5 * (p[i - 1] + p[i]), 0.5 * (tau[i - 1] + tau[i])
    err = abs(1.0 / ray_scale(p_m, tau_m, r, p_c, tau_c) - 1.0)

    # --- left: normalised, shape only ---
    pp = np.linspace(-1.3 * p_c, p_c, 800)
    axes[0].plot(pp / p_c, surface_tau(pp, r, p_c, tau_c) / tau_c, "-", marker="",
                 lw=1.8, color=colour,
                 label=rf"$r = {r}$ ({label}), {err:.2%}")
    if r == "1":
        axes[0].plot([1, 1], [0, 1], "-", marker="", lw=1.8, color=colour)
    axes[0].plot(p_m / p_c, tau_m / tau_c, marker="o", ms=7, ls="none",
                 color=colour, mec="k", mew=0.7)

    # --- right: physical stresses ---
    th = np.linspace(0, 2 * np.pi, 721)
    exx, eyy = np.cos(th), np.sin(th)
    tr = exx + eyy
    p_dir = kappa * tr
    tau_dir = 2 * mu * np.sqrt(2 / 3 * (exx**2 + eyy**2 - exx * eyy))
    Sxx, Syy = [], []
    for k in range(th.size):
        a = (ray_scale(p_dir[k], tau_dir[k], r, p_c, tau_c) if p_dir[k] >= 0
             else tau_c / tau_dir[k])
        Sxx.append(p_dir[k] * a + 2 * mu * a * (exx[k] - tr[k] / 3))
        Syy.append(p_dir[k] * a + 2 * mu * a * (eyy[k] - tr[k] / 3))
    axes[1].plot(np.array(Sxx) / M, np.array(Syy) / M, "-", marker="", lw=1.8,
                 color=colour, label=rf"$r = {r}$, $p_c = {p_c/M:.1f}$ MPa")
    axes[1].plot(0.5 * (sxx[i - 1] + sxx[i]) / M, 0.5 * (syy[i - 1] + syy[i]) / M,
                 marker="o", ms=7, ls="none", color=colour, mec="k", mew=0.7)

axes[0].set_xlabel(r"$p\,/\,p_c$")
axes[0].set_ylabel(r"$\tau\,/\,\tau_c$")
axes[0].set_xlim(-1.3, 1.25)
axes[0].set_ylim(0, 1.1)
axes[0].legend(loc="lower left", fontsize=9, title="measured vs Eq. (110)",
               title_fontsize=9)
# For p < 0 every r gives tau = tau_c (Eq. 109), so the three curves coincide
# there and only the last drawn is visible.
axes[0].set_title("(a) shape, normalised")
axes[0].grid(alpha=0.3)

axes[1].axhline(0, lw=0.6, color="0.7", marker="")
axes[1].axvline(0, lw=0.6, color="0.7", marker="")
axes[1].set_xlabel(r"$\sigma_{xx}$ (MPa)")
axes[1].set_ylabel(r"$\sigma_{yy}$ (MPa)")
axes[1].legend(loc="lower right", fontsize=9)
axes[1].set_title("(b) physical stresses")
axes[1].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(PNG)
print(f"[INFO] wrote {PNG}")
