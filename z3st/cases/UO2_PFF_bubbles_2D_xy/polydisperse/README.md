# Polydisperse periodic row (large + small bubble per period)

## Why a bubble-centred cell still works

The monodisperse `periodic_row/` case cuts each bubble in half with `Clamp_x` at its
own centre, relying on that centre being a mirror plane of the infinite row. For
identical, evenly spaced bubbles it is. For bubbles of different sizes it is not: the
spacing on the two sides of a given bubble is no longer equal, so its centre stops
being a symmetry plane of the row.

The fix is to alternate the sizes: A, B, A, B, ... A row built this way has a mirror
plane at the centre of every bubble whatever its size, because each lentille is
symmetric about its own vertical median and the two neighbours of any bubble are then
the same size at the same distance. The unit cell runs from the centre of an A to the
centre of a B, and `Clamp_x` on `xmin` and `xmax` is exact, not an approximation.

## Boundary conditions

`Clamp_y` on `ymin` (the grain-boundary plane), `Clamp_x` on `xmin` (A's median) and
on `xmax` (B's median), Neumann pressure on `cavity_A` and `cavity_B`. `xmax` is a
symmetry plane, not a free edge — leaving it free would turn the model into one
isolated pair of bubbles in a finite plate, a different problem with weaker
confinement. `periodic_row/README.md` records the size of that difference.

## Results

`L_x` is adjusted per configuration so the ligament stays at 22.6 um, which is the
`Lx_45` value from `periodic_row/`. The ratio-1.00 point is that monodisperse case
reused, not a separate run. With the ligament fixed, any change comes from the size
difference alone.

| R_p,A / R_p,B | R_p,B (um) | L_x (um) | p_crit at D_max = 0.5 (MPa) | A saturates -> B initiates (MPa) |
|---|---|---|---|---|
| 1.00 | 11.2 | 45.0 | 422.5 | 0 |
| 1.33 | 8.4  | 42.2 | 450.6 | 4 |
| 2.00 | 5.6  | 39.4 | 476.0 | 20 |

The large bubble initiates first and alone; the small one damages later and on its own
concentration, with 15.2 um of intact material still between them at ratio 2.00. A
mixed pair is therefore stronger than a uniform pair at equal ligament: replacing one
bubble with a smaller one raises p_crit by 12.7 % at ratio 2.00. The gain per unit of
ratio falls between the two intervals, which two effects explain equally well at three
points — B's shrinking cross-section, and B approaching `l_ch` so the phase field
resolves it less and less as a defect. Separating them would need l varied at fixed
geometry, which was not done.

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
