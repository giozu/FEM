# Ramp-rate independence check — physical UO2 anchor

## Why this directory exists

The p_crit ≈ 200-210 MPa result (physical calibration, `sigma_c=150 MPa`
from `materials/uo2_sigma_anchor.yaml`, not the K_IC-anchored variant used
elsewhere in this case) was originally produced in scratch space
during Étape 1 and never saved anywhere in the repo. This reruns it
properly, persisted, so the "ramp-rate independence" claim in the
report has a traceable source.

## Setup

Single isolated lenticular bubble, quarter-model — identical geometry
to `theta_sweep/theta_48p1/` (`theta_deg=48.1`, `Rp=11.2e-6`,
`Lx=Ly=60e-6`, `h_cavity=1.487e-7` i.e. the resolved theta=48.1 value,
`h <= lc/4` satisfied). AT1 + `amor` split, `hybrid_constraint: true`,
`lc=2.0e-6`. Material: `../../../../materials/uo2_sigma_anchor.yaml` (the
real physical `sigma_c=150 MPa`, `Gc=0.335 J/m²` auto-derived — no
calibration-anchor trick needed here, this is the physical value).
Pressure ramp `0 -> -300e6 Pa` on `cavity`, identical target range in
all three runs, only `n_steps` differs.

## Result

Two series. The first is what the internship measured, on the core code at tag
`baptiste-final-2026-09`. The second is the same three cases re-run on `develop`
(2026-09-12) with the same mesh and the same inputs — only `n_steps` differs between the
three, verified by diff.

| Run | n_steps | Pa/step | p_crit at tag | p_crit on develop |
|---|---|---|---|---|
| `n20`  | 21  | 15.0 MPa | -214.8 MPa | **-202.68 MPa** |
| `n150` | 151 | 2.0 MPa  | -203.0 MPa | **-192.30 MPa** |
| `n201` | 201 | 1.5 MPa  | -202.0 MPa | **-191.31 MPa** |

**The ramp-rate conclusion holds, and is sharper on develop.** `n150` and `n201` agree to
0.52 %, so the converged value is 191-192 MPa. `n20` sits 5.9 % above it. As the tag-side
run already showed, that gap is an interpolation artefact, not ramp-rate physics: the
crossing is bracketed by two points 15 MPa apart (D = 0.299 at -195 MPa, D = 0.692 at
-210 MPa) either side of a very sharp transition. A 20-step ramp locates p_crit to a few
per cent, not to the MPa.

**The absolute level moved because the Amor split was corrected.** All three shifted down
by 5.3-5.6 %, a uniform offset rather than scatter. The cause is commit `dbb6fae`
(2026-07-07): `psi_amor_split` used the Lame parameter `lambda` as the coefficient of the
volumetric term, and now uses Amor's n-dimensional bulk modulus `K_n = lambda + 2*mu/n`
(Amor, Marigo & Maurini 2009), which is what the split is defined with. At `E = 358 GPa`,
`nu = 0.23`, `dim = 2` that coefficient rises by 2.174x, and the measured `psi+` at the
tip rises by 1.115x — consistent with a tip state roughly 10 % volumetric.

Verified directly: setting `k_n = lam` in `damage_model.py` and re-running `n20` returns
-214.795 MPa against the -214.8 MPa recorded at the tag, with the same two bracketing
steps (D = 0.376 at -210 MPa, D = 0.765 at -225 MPa). The revert was temporary and is not
in the tree.

**The develop values are the correct ones.** The tag-side numbers were produced with a
volumetric coefficient that was too small by a factor 2.174, so they overestimate p_crit.

Consequence for the report: its stated band of 200-210 MPa across 20/150/201 steps is the
tag-side reading, and even there only `n20` and the two converged runs at 202-203 fall
inside it. On develop the converged value is 191-192 MPa, below the band. Any figure
quoted from this case must name the code state it was produced on.

The same caution applies to the rest of `UO2_PFF_bubbles_2D_xy`: every number in it was
produced at the tag, and none has been re-run on develop except this case. Nothing here
should be blessed as a gold until that is done.
