#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Non-regression for verification/cohesive/strength_surface_2D.

Checks that the model nucleates where its strength surface says it should.
The block is driven through a homogeneous stress state until the phase field
first becomes non-zero; that stress must lie on dS_0, which for the r-norm
family is Eq. (110) of Vicentini et al. (2026),

    (p / p_c)^(r/(r-1)) + (tau / tau_c)^(r/(r-1)) = 1,    p >= 0,

with the r = 1 and r = inf limits being the rectangle of Eq. (100) and the
Drucker-Prager line of Eq. (112). This is the property the whole paper is
built around -- a strength surface that can be prescribed rather than
inherited from an energy split -- and the 1D bar case does not test it at all.

Load stepping resolves the onset only to within one increment, so the measured
point is the midpoint of the bracket (last elastic step, first damaged step)
and the half-bracket is reported as the resolution. The check that the analytic
surface falls *inside* that bracket is the sharper of the two tests, since it
cannot be satisfied by a compensating error in the step size.
"""

import os

import numpy as np

from z3st.utils.non_regression import case_paths, error_metric, finish, load_case, metric, tracked

CASE_DIR, _VTU, OUT_JSON = case_paths(__file__)
CSV = os.path.join(CASE_DIR, "output", "response.csv")

# One load increment is 0.5% of the nucleation stress, so the midpoint of the
# bracket carries about half that. The tolerance leaves room for it.
TOLERANCE = 0.02
ALPHA_ONSET = 1.0e-6        # the paper's own threshold for "damage has appeared"

geom, inp, mat = load_case(CASE_DIR)
r = str(inp["models"]["cohesive"]["r_norm"])
p_c, tau_c = float(mat["p_c"]), float(mat["tau_c"])


def surface_residual(p, tau):
    """Eq. (110) residual: zero on the strength surface, negative inside."""
    if r == "1":                      # rectangle, Eq. (100)
        return max(p / p_c, tau / tau_c) - 1.0
    if r == "inf":                    # Drucker-Prager line, Eq. (112)
        return p / p_c + tau / tau_c - 1.0
    e = float(r) / (float(r) - 1.0)   # Lame curve
    return (p / p_c) ** e + (tau / tau_c) ** e - 1.0


data = np.genfromtxt(CSV, delimiter=",", names=True)
p = np.atleast_1d(data["p_Pa"])
tau = np.atleast_1d(data["tau_Pa"])
alpha = np.atleast_1d(data["alpha_max"])

damaged = np.flatnonzero(alpha > ALPHA_ONSET)
if damaged.size == 0:
    raise RuntimeError(
        f"The block never nucleated (max alpha = {alpha.max():.3e}); "
        f"the load range does not reach the strength surface."
    )
i = int(damaged[0])
if i == 0:
    raise RuntimeError("Damage is present at the first step; the ramp starts too high.")

p_lo, tau_lo = p[i - 1], tau[i - 1]     # last elastic state, inside dS_0
p_hi, tau_hi = p[i], tau[i]             # first damaged state, outside
p_mid, tau_mid = 0.5 * (p_lo + p_hi), 0.5 * (tau_lo + tau_hi)

res_lo, res_hi = surface_residual(p_lo, tau_lo), surface_residual(p_hi, tau_hi)
brackets = res_lo < 0.0 < res_hi

# How far the measured point sits from the surface, measured along its own ray:
# the scale factor s that puts (p, tau) / s onto dS_0. The residual decreases
# with s, so the bracket [0.5, 2.0] straddles the root for any sane result.
ray = np.hypot(p_mid, tau_mid)
lo, hi = 0.5, 2.0
if not (surface_residual(p_mid / lo, tau_mid / lo) > 0
        > surface_residual(p_mid / hi, tau_mid / hi)):
    raise RuntimeError(
        "The nucleation point is not within a factor of two of the strength "
        "surface; the model and the reference disagree about more than the "
        "load-step resolution."
    )
for _ in range(60):
    s = 0.5 * (lo + hi)
    if surface_residual(p_mid / s, tau_mid / s) < 0:
        hi = s
    else:
        lo = s
radial_error = abs(0.5 * (lo + hi) - 1.0)

half_bracket = 0.5 * abs(np.hypot(p_hi, tau_hi) - np.hypot(p_lo, tau_lo)) / ray

print(f"[INFO] r = {r}, p_c = {p_c:.3e} Pa, tau_c = {tau_c:.3e} Pa")
print(f"[INFO] nucleation bracketed at step {i}: "
      f"p = {p_lo/1e6:.4f} -> {p_hi/1e6:.4f} MPa, "
      f"tau = {tau_lo/1e6:.4f} -> {tau_hi/1e6:.4f} MPa")
print(f"[INFO] Eq. (110) residual: {res_lo:+.4e} (inside) -> {res_hi:+.4e} (outside)")
print(f"[INFO] radial distance from the surface = {radial_error:.4%}, "
      f"load-step resolution = {half_bracket:.4%}")

errors = {
    # The measured nucleation point must sit on the surface.
    "nucleation_on_surface_radial": error_metric(radial_error),
    # Sharper: the surface must fall strictly between the last elastic and the
    # first damaged state. Step size cannot compensate for a wrong surface.
    "surface_brackets_nucleation": error_metric(0.0 if brackets else 1.0),
    "p_nucleation_Pa": tracked(float(p_mid)),
    "tau_nucleation_Pa": tracked(float(tau_mid)),
    "nucleation_step": tracked(float(i)),
    # p and tau track each other by construction at this loading angle, which
    # is what makes the point discriminate between the r-norms.
    "p_over_tau": metric(float(p_mid / tau_mid), 1.0),
    "load_step_resolution": tracked(float(half_bracket)),
}

finish(errors, TOLERANCE, OUT_JSON, CASE_DIR)
