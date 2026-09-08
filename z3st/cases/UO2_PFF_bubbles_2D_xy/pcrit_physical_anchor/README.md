# Ramp-rate independence check — physical UO2 anchor

## Why this directory exists

The p_crit ≈ 200-210 MPa result (physical calibration, `sigma_c=150 MPa`
from `materials/uo2.yaml`, not the K_IC-anchored variant used
elsewhere in this case) was originally produced in scratch space
during Étape 1 and never saved anywhere in the repo. This reruns it
properly, persisted, so the "ramp-rate independence" claim in the
report has a traceable source.

## Setup

Single isolated lenticular bubble, quarter-model — identical geometry
to `theta_sweep/theta_48p1/` (`theta_deg=48.1`, `Rp=11.2e-6`,
`Lx=Ly=60e-6`, `h_cavity=1.487e-7` i.e. the resolved theta=48.1 value,
`h <= lc/4` satisfied). AT1 + `amor` split, `hybrid_constraint: true`,
`lc=2.0e-6`. Material: `../../../../materials/uo2.yaml` directly (the
real physical `sigma_c=150 MPa`, `Gc=0.335 J/m²` auto-derived — no
calibration-anchor trick needed here, this is the physical value).
Pressure ramp `0 -> -300e6 Pa` on `cavity`, identical target range in
all three runs, only `n_steps` differs.

## Result

| Run | n_steps | Pa/step | p_crit (D_max=0.5) |
|---|---|---|---|
| `n20` | 21 | 15.0 MPa | **-214.8 MPa** |
| `n150` | 151 | 2.0 MPa | **-203.0 MPa** |
| `n201` | 201 | 1.5 MPa | **-202.0 MPa** |

See `ramp_rate_independence.png`.

**The claimed "200-210 MPa across 20/150/201 steps" is mostly right
but not exact**: `n150` and `n201` agree to <1 MPa and sit inside the
window; `n20` lands at 214.8 MPa, ~2.3% above the stated upper bound.
Looking at the actual data points (not just the interpolated
crossing), `n20`'s trajectory tracks the finer curves closely up to
the point right before nucleation — the 214.8 MPa figure comes from
linearly interpolating between two points 15 MPa apart (D=0.38 at
-210 MPa, D=0.77 at -225 MPa) that straddle a very sharp transition.
That is a **resolution artifact of the interpolation, not evidence of
a real ramp-rate-dependent shift in the physics** — but it does mean
a 20-step ramp alone should not be quoted as precise to the MPa; the
converged value (150 vs 201 steps) is ≈ 202-203 MPa.

Suggested correction for the report: *"p_crit converges to ≈202-203
MPa at 150 and 201 steps (agreement <1 MPa); a coarse 20-step ramp
still locates it within ~5% (214.8 MPa) but is limited by
interpolation across a 15 MPa gap, not by a real ramp-rate
sensitivity — unlike the load-step sensitivity seen on case 15."*
