# Polydisperse periodic row (large + small bubble per period)

## Why the quarter-model breaks here

The monodisperse `periodic_row/` case cuts each bubble in half with
`Clamp_x` at its own center, relying on that center being a valid
mirror plane for the infinite row. That requires equal spacing on
both sides of every bubble — true for identical, evenly-spaced
bubbles, but false as soon as neighbors differ in size: forcing
`Clamp_x` at bubble 1's center would force the material on both sides
of it to be identical, which contradicts having a differently-sized
bubble 2 nearby.

## What this model does instead — and what "free-edge approximation" means

**Half-model only**: `Clamp_y` on `ymin` (the grain-boundary plane
stays a valid symmetry plane regardless of bubble size — each bubble
individually is still symmetric top/bottom) and **no `Clamp_x`
anywhere**. `xmin`/`xmax` are tagged physical curves but carry no BC —
they are left as free (natural, zero-traction) surfaces.

**Concretely, what this means and why it's an approximation**: a
*true* periodic model would tie `xmin` and `xmax` together (matching
displacement/traction across them) so that the elastic field "sees"
the next bubble in the row continuing beyond each edge — that is what
a periodic boundary condition (PBC) does, and it is what the
`periodic_row/` quarter-model achieved *implicitly* via `Clamp_x` at
each bubble's own symmetric center. Here, with no PBC support assumed
available and no valid mirror plane to exploit, `xmin`/`xmax` are left
literally free instead: the model becomes "one isolated pair of
bubbles embedded in a finite plate with free left/right edges," not
"one period of an infinite row." The two situations are NOT the same
problem: a free edge lets the plate locally relax/expand there, which
a real periodic neighbor bubble would not allow (see the
`periodic_row/README.md` finding that changing a boundary from
`Clamp_x` to free shifted `p_crit` by a comparable amount to the
periodicity effect itself — free vs. clamped edges are not
interchangeable here). With `half_B = 15e-6` clearance on each side
(same order of magnitude as `ligament_A`, not a large multiple of it),
this free-edge simplification is not obviously negligible. It has not
been quantified in this pass (geometry check only, per request — no BC
file, no run yet); a follow-up should either implement true PBCs or at
least bound the error by comparing against a much larger `half_B`.

## Layout (`mesh.geo`, current defaults)

Two full lentilles (not cut in x), each individually clipped by `y=0`
only, placed along the grain boundary:

```
|<- half_B ->|<- bubble1 (2*Rp_1) ->|<- ligament_A ->|<- bubble2 (2*Rp_2) ->|<- half_B ->|
x=0        15.0                  37.4              52.4                  63.6        78.6  (um)
```

- `Rp_1 = 11.2e-6`, `Rp_2 = 5.6e-6` (factor 2), same `theta_deg = 48.1` for both.
- `ligament_A = 15e-6` (direct gap, bubble1-to-bubble2) and
  `ligament_B = 30e-6` (wrap-around gap, split evenly as `half_B=15e-6`
  clearance on each free edge) — chosen distinct (factor 2, matching
  the bubble-size ratio) and both comfortably above the
  `eps=0.5e-6` bounding-box margin and the `lc=2e-6` regularization
  length (7.5·lc and 15·lc respectively).
- `Ly = 60e-6` unchanged.

Verified via the mesh log Printf: `ymin` splits into 3 physical
segments (before bubble1, `ligament_A`, after bubble2 — 1 curve each,
3 total), `cavity` is 2 curves (1 full arc per bubble), `xmin`/`xmax`/
`ymax` are 1 curve each — all as expected.

## Mesh resolution: no per-bubble curvature criterion (correction)

An earlier version of this README claimed the same `h_cavity` under-
resolves bubble 2 relative to bubble 1 because "`rho_2 = rho_1/2`
(tip curvature radius scales with `Rp`)". **That premise is wrong**:
`rho = ay^2/Rp` is the curvature radius of an *ellipse* of semi-axes
`Rp` and `ay`. Our cavity is a **lentille** — two circular arcs
crossing at a corner of opening angle `2*theta_deg = 96.2 deg`. A
corner has no curvature radius; it is a singular point (the elastic
stress field there scales as `r^-0.45`, a wedge/notch singularity, not
a curvature-controlled concentration). There is therefore no
"resolve `rho`" mesh criterion to apply per bubble, and no resolution
mismatch between bubble 1 and bubble 2 to worry about on those
grounds.

The actual mesh constraint is simpler and size-independent: refine
proportionally to the size of the defect, subject to
**`h <= lc/4 = 0.5e-6`** (the phase-field regularization length sets
the floor, not the cavity geometry). `h_cavity = 0.15e-6` satisfies
this comfortably for both bubbles regardless of `Rp_1`/`Rp_2`. The
Étape-3 mesh-convergence study (`periodic_row/mesh_convergence_Lx45.*`)
already demonstrates `p_crit` stable to <1% between `h=0.30e-6` and
`h=0.075e-6`, both well inside this constraint — no new convergence
check is implied by the two bubbles having different `Rp`.

(The same correction applies to the Étape-2 angle sweep, which used a
`rho/h=15` criterion under the same wrong ellipse-curvature premise —
see `theta_sweep`/`periodic_row` READMEs. Its results are unaffected:
`h` there was always well below `lc/4` too, and the direct
mesh-convergence check confirms `p_crit` convergence independently of
that flawed reasoning.)

## Point to note: the small bubble is near the model's resolution floor

`Rp_2 = 5.6e-6` should be compared to the AT1 characteristic length
`l_ch = (8/3)*lc = (8/3)*2e-6 = 5.33e-6`. `Rp_2` is only ~5% above
`l_ch` — the small bubble sits right at the edge of what this
regularized model can represent as a distinct "defect" rather than a
feature smeared into the diffuse process zone. **If the small bubble
turns out not to participate in cracking (e.g. all damage localizes
at the large bubble's tip and the ligament next to it), that would be
a genuine physical/model result of this size ratio — not a meshing
bug.** Worth watching for explicitly once a run is performed.
