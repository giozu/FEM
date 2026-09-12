#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST: An open-source FEniCSx framework for thermo-mechanical analysis
# Author: Giovanni Zullo
# Version: 0.3.2 (2026)
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Diagnostics for verification/cohesive/strength_surface_2D.

One row per step in ``output/response.csv``: the volume-averaged stress in the
two forms the paper plots it (the pressure-shear pair and the in-plane normal
components) together with the maximum of the phase field.

Averaging over the block is exact while the state is homogeneous, which is the
whole of the pre-nucleation history -- and the pre-nucleation history is what
this case measures. Once damage localizes the average stops meaning anything,
which is why non-regression.py reads only the first row with alpha > 0.

The pressure and the shear stress are the pair the strength surface is written
in (paper Eq. 89-90): sigma = p I + tau e_dev, so p = tr(sigma)/3 and
tau = ||dev(sigma)||.
"""

import os

import dolfinx
import numpy as np
import ufl
from mpi4py import MPI

_CSV = os.path.join(os.path.dirname(__file__), "output", "response.csv")
_HEADER = "step,U_t_m,p_Pa,tau_Pa,sigma_xx_Pa,sigma_yy_Pa,alpha_max\n"

_run_started = False
_forms = None


def _build_forms(problem):
    """Volume-average forms for the stress components, built once."""
    material = problem.materials["block"]
    u, tr_eta, dev_eta = problem.split_state(problem.w)
    sigma = problem.sigma_cohesive(u, tr_eta, dev_eta, material)
    dx = problem.dx_tags[problem.label_map["block"]]

    p = ufl.tr(sigma) / 3.0
    dev = sigma - p * ufl.Identity(3)
    tau = ufl.sqrt(ufl.inner(dev, dev))

    one = dolfinx.fem.Constant(problem.mesh, dolfinx.default_scalar_type(1.0))
    return {
        "area": dolfinx.fem.form(one * dx),
        "p": dolfinx.fem.form(p * dx),
        "tau": dolfinx.fem.form(tau * dx),
        "sxx": dolfinx.fem.form(sigma[0, 0] * dx),
        "syy": dolfinx.fem.form(sigma[1, 1] * dx),
    }


def per_step(problem, step, t):
    global _run_started, _forms

    if _forms is None:
        _forms = _build_forms(problem)

    comm = problem.mesh.comm

    def avg(key, area):
        return comm.allreduce(dolfinx.fem.assemble_scalar(_forms[key]), op=MPI.SUM) / area

    area = comm.allreduce(dolfinx.fem.assemble_scalar(_forms["area"]), op=MPI.SUM)
    p, tau = avg("p", area), avg("tau", area)
    sxx, syy = avg("sxx", area), avg("syy", area)

    # Prescribed magnitude U_t: the BC constants hold their signed components.
    ux = uy = 0.0
    for bc_list in problem.dirichlet_mechanical.values():
        for bc in bc_list:
            if not (isinstance(bc, dict) and isinstance(bc.get("raw"), list)):
                continue
            v = float(np.atleast_1d(bc["const"].value)[0])
            if bc["id"] == problem.label_map["right"]:
                ux = v
            elif bc["id"] == problem.label_map["top"]:
                uy = v
    U_t = float(np.hypot(ux, uy))

    alpha_max = comm.allreduce(float(problem.D.x.array.max()), op=MPI.MAX)

    if comm.rank != 0:
        return

    os.makedirs(os.path.dirname(_CSV), exist_ok=True)
    with open(_CSV, "w" if not _run_started else "a") as f:
        if not _run_started:
            f.write(_HEADER)
        f.write(f"{step:d},{U_t:.10e},{p:.10e},{tau:.10e},"
                f"{sxx:.10e},{syy:.10e},{alpha_max:.10e}\n")
    _run_started = True
