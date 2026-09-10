# SPDX-License-Identifier: Apache-2.0
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
# Z3ST <-> OpenMC coupling — axial power form factor from a mesh tally
# Author: Giovanni Zullo
# --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. --- --.. ..- .-.. .-.. ---
"""
Extract the axial power shape along one fuel element from an OpenMC statepoint
and write it as a Z3ST material-card fragment.

The coupling is one-way, neutronics -> Z3ST, and the exchange is a file: Z3ST
never imports OpenMC. This script runs in the OpenMC environment, its output is
pasted into the fuel card, and the case then runs with no OpenMC dependency.

The fragment feeds the tabulated axial bus already in Z3ST::

    axial_profile: materials.fuel_profiles.tabulated_axial
    axial_table_z: [...]   # (m) elevations, z = 0 at the bottom of active fuel
    axial_table_f: [...]   # (-) relative power

``set_power`` normalises the composite form factor to nodal mean 1, so only the
shape of the table is used, not its scale.

What this script does NOT provide is the absolute rating. A Cartesian core mesh
tally gives no per-element power fraction: the mesh cells do not follow the
element boundaries, and a flux tally is not a power tally. The linear heat rate
stays the case input ``lhr:``, set from the operating power and the element
peaking factor. To get the rating from OpenMC instead, tally ``fission-q-prompt``
(or ``heating-local``) over the fuel cell of the element and normalise with the
core-wide ``heating-local``.

Usage::

    python axial_power.py statepoint.300.h5 --pin-x 0 --pin-y 8.2 \
        --z-active 8.81 46.91 -o fuel_axial.yaml

Coordinates and lengths on the command line are in cm, the OpenMC convention;
the emitted table is in m, the Z3ST convention.
"""

import argparse
import numpy as np


def mesh_column(tally, mesh, x_cm, y_cm, radius_cm=0.0):
    """Mean and standard deviation per axial mesh row, summed over the mesh
    columns that fall within ``radius_cm`` of (x_cm, y_cm).

    A radius of 0 selects the single column containing the point. Widening it
    averages neighbouring columns, which trades axial resolution of the element
    for statistics when the tally is noisy.

    Returns (z_centres_cm, mean, std_dev), each of length nz.
    """
    nx, ny, nz = mesh.dimension
    ll = np.asarray(mesh.lower_left, dtype=float)
    ur = np.asarray(mesh.upper_right, dtype=float)

    # Mesh-filter bins run with x fastest, then y, then z.
    axis = [i for i, f in enumerate(tally.filters)
            if type(f).__name__ == "MeshFilter"]
    if not axis:
        raise ValueError(f"tally '{tally.name}' has no MeshFilter")
    axis = axis[0]

    def collapse(value):
        data = tally.get_reshaped_data(value=value)
        if value == "std_dev":
            # Variances add, standard deviations do not.
            data = data ** 2
        keep = data
        for ax in sorted(range(keep.ndim), reverse=True):
            if ax != axis:
                keep = keep.sum(axis=ax)
        return keep.reshape(nz, ny, nx)

    mean = collapse("mean")
    var = collapse("std_dev")

    xc = ll[0] + (np.arange(nx) + 0.5) * (ur[0] - ll[0]) / nx
    yc = ll[1] + (np.arange(ny) + 0.5) * (ur[1] - ll[1]) / ny
    zc = ll[2] + (np.arange(nz) + 0.5) * (ur[2] - ll[2]) / nz

    if radius_cm > 0.0:
        dx = xc[None, :] - x_cm
        dy = yc[:, None] - y_cm
        sel = (dx ** 2 + dy ** 2) <= radius_cm ** 2
        if not sel.any():
            raise ValueError(
                f"no mesh column within {radius_cm} cm of ({x_cm}, {y_cm}); "
                f"mesh cell size is {(ur[0]-ll[0])/nx:.2f} x {(ur[1]-ll[1])/ny:.2f} cm"
            )
    else:
        sel = np.zeros((ny, nx), dtype=bool)
        sel[int(np.argmin(np.abs(yc - y_cm))),
            int(np.argmin(np.abs(xc - x_cm)))] = True

    col = mean[:, sel].sum(axis=1)
    col_sd = np.sqrt(var[:, sel].sum(axis=1))
    return zc, col, col_sd


def to_table(z_cm, f, sd, z0_cm, z1_cm):
    """Restrict to the active fuel, shift z to the bottom of active fuel, convert
    to m, and normalise the shape to mean 1 over the retained rows."""
    inside = (z_cm >= z0_cm) & (z_cm <= z1_cm)
    if inside.sum() < 2:
        raise ValueError(
            f"only {inside.sum()} mesh rows inside the active range "
            f"[{z0_cm}, {z1_cm}] cm; the mesh is too coarse axially"
        )
    z = (z_cm[inside] - z0_cm) / 100.0
    v = f[inside]
    s = sd[inside]
    mean = v.mean()
    if mean <= 0:
        raise ValueError("non-positive mean over the active range")
    return z, v / mean, s / v


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("statepoint")
    p.add_argument("--tally", default="Phi(r)",
                   help="name of the mesh tally to read (default: Phi(r))")
    p.add_argument("--pin-x", type=float, required=True, help="(cm)")
    p.add_argument("--pin-y", type=float, required=True, help="(cm)")
    p.add_argument("--radius", type=float, default=0.0,
                   help="(cm) average over columns within this radius of the "
                        "element centre; 0 takes the single containing column")
    p.add_argument("--z-active", type=float, nargs=2, required=True,
                   metavar=("Z0", "Z1"), help="(cm) active fuel range")
    p.add_argument("-o", "--output", default="fuel_axial.yaml")
    args = p.parse_args()

    import openmc

    sp = openmc.StatePoint(args.statepoint)
    tally = sp.get_tally(name=args.tally)
    mesh = next(f.mesh for f in tally.filters if type(f).__name__ == "MeshFilter")

    z_cm, col, col_sd = mesh_column(tally, mesh, args.pin_x, args.pin_y, args.radius)
    z, f, rel_sd = to_table(z_cm, col, col_sd, *args.z_active)

    peak = float(f.max())
    worst = float(np.max(rel_sd))
    lines = [
        "# Axial power form factor for a Z3ST fuel card.",
        "# Generated by z3st/coupling/openmc/axial_power.py — do not edit by hand.",
        f"# statepoint : {args.statepoint}",
        f"# tally      : {args.tally!r}, scores {tally.scores}",
        f"# k-effective: {sp.keff}",
        f"# element    : (x, y) = ({args.pin_x}, {args.pin_y}) cm, "
        f"radius {args.radius} cm",
        f"# active fuel: {args.z_active[0]} to {args.z_active[1]} cm "
        f"({len(z)} mesh rows)",
        f"# axial peaking factor: {peak:.3f}",
        f"# worst per-row relative standard deviation: {worst*100:.2f} %",
        "",
        "axial_profile: materials.fuel_profiles.tabulated_axial",
        "axial_table_z: [" + ", ".join(f"{v:.5f}" for v in z) + "]",
        "axial_table_f: [" + ", ".join(f"{v:.4f}" for v in f) + "]",
        "",
    ]
    with open(args.output, "w") as fh:
        fh.write("\n".join(lines))

    print(f"[axial_power] {len(z)} rows written to {args.output}")
    print(f"[axial_power] axial peaking factor {peak:.3f}, "
          f"worst row sigma {worst*100:.2f} %")
    if worst > 0.05:
        print("[axial_power] WARNING: the shape is poorly converged; "
              "raise the particle count or widen --radius")


if __name__ == "__main__":
    main()
