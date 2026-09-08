# Periodic row of lenticular cavities (Étape 3)

## Physical configuration

This case models an **infinite periodic row of identical lenticular
cavities** along a grain boundary (`y = 0`), spaced `Lx` apart
center-to-center. It is **not** a pair of isolated bubbles.

The quarter-symmetry unit cell exploits two mirror planes at once:

- `y = 0` (`ymin`): the grain-boundary plane itself, standard
  quarter-model symmetry, as in the single-bubble case.
- `x = 0` (`xmin`): the mid-plane of one bubble.
- `x = Lx` (`xmax`): **also a mirror plane**, not a free surface. It
  sits at the mid-plane of the *next* bubble in the row.

The unit cell therefore contains **two half-lentilles**: one centered
at `x = 0` (cut by `xmin`) and one centered at `x = Lx` (cut by
`xmax`). `Lx` **is** the center-to-center spacing between neighboring
bubbles in the unfolded periodic row (unfold across `xmin`: bubble
centered at 0; unfold across `xmax`: next bubble centered at `Lx`).
Mirroring this cell repeatedly across `xmin`/`xmax` reconstructs the
full infinite row with spacing `Lx`.

**Note on an earlier, corrected convention**: an intermediate version
of this case used a separate parameter `s` with `Lx = s/2`, on the
(wrong) assumption that `Lx` was only half the bubble spacing. `Lx`
alone already **is** the spacing; `s` was redundant and has been
removed from `mesh.geo` and all scripts. Directories that were
previously named `s_XXX` are now named `Lx_YYY` with `YYY = XXX/2`
(the numerical mesh/results are unchanged — this was a naming/labeling
fix only, no case was rerun). `ligament_um` values are unaffected by
this fix since they were always computed from the correct `Lx`.

**Consequence for boundary conditions**: `xmax` must carry `Clamp_x`
(zero normal displacement), exactly like `xmin`. It is no longer a
traction-free edge as it was in the single-bubble case — a free
`xmax` would silently turn this into a *finite pair of bubbles in a
finite plate* instead of an infinite periodic row, which is a
different physical problem with a different (weaker) confinement.

## Geometry parameters (`mesh.geo`)

| Parameter | Value | Notes |
|---|---|---|
| `Rp` | 11.2e-6 m | cavity projected radius, fixed |
| `theta_deg` | 48.1 deg | semi-dihedral angle, fixed |
| `Lx` | overridable, default 60e-6 m | center-to-center bubble spacing = quarter-cell width |
| `Ly` | 60e-6 m | unchanged |
| `h_cavity` | 0.15e-6 m | retained from the Étape 2 angle sweep; refinement proportional to defect size, subject to `h <= lc/4 = 0.5e-6` (see correction note below) |

The **ligament** (solid material between the two facing tips) is
`Lx - 2*Rp`. This must stay positive and, per the `eps = 0.5e-6`
margin used in the bounding-box curve selection, well above
`~1e-6 m` for the `cavity`/`ymin` physical groups to be selected
unambiguously.

**Correction on mesh-resolution criterion**: the Étape 2 angle sweep
scaled `h_cavity` per `theta_deg` to keep `rho/h_cavity = 15` constant,
on the premise that `rho = ay^2/Rp` is the tip curvature radius. That
premise is wrong: the cavity tip is a wedge (two circular arcs meeting
at opening angle `2*theta_deg`), not an ellipse — a corner has no
curvature radius, it is a singular point (`stress ~ r^-0.45`). There
is no "resolve `rho`" criterion. The only real constraint is
`h <= lc/4 = 0.5e-6` (phase-field regularization length), which
`h_cavity = 0.15e-6` satisfies comfortably regardless of `theta_deg`
or `Rp`. **Results are unaffected**: the direct mesh-convergence check
below (independent of the `rho/h` reasoning) already shows `p_crit`
stable to <1% between `h=0.30e-6` and `h=0.075e-6`. See
`polydisperse/README.md` for the fuller correction and its
consequence there.

## Known-invalid regime

For `Rp = 11.2e-6` fixed, `Lx < 2*Rp + margin` (roughly `Lx < 24e-6`)
makes the two half-lentilles overlap: the ligament length goes
negative, `BooleanDifference` still runs but merges the two cavities
into a single notch reaching all the way through the ligament (see
`mesh_s36_full.png`, generated for `Lx = 18e-6`,
`2*Rp = 22.4e-6 > Lx` — file kept from an earlier exploratory pass,
predates the `s`->`Lx` naming fix). In that regime `ymin` disappears
entirely (0 curves — verified in the mesh log) since there is no
longer any solid boundary left at `y = 0`. This is not a bug in the
script; it is the geometry correctly reporting that the requested
spacing is physically inadmissible for this bubble size.

## Non-regression target (isolated-bubble limit)

As `Lx -> infinity`, the periodic-row result must converge to the
Étape 1/2 isolated-bubble result (`p_crit ~ 571 MPa` at
`theta_deg=48.1`, from the angle sweep).

**Result of the `Lx=200e-6` run: -616.3 MPa, not ~571 MPa (+7.9%).**
Root cause, isolated via two extra diagnostic runs (single bubble, no
periodicity, `Lx=200e-6`): the elastic tip-stress concentration factor
itself drops from 9.58 (`Lx=Ly=60e-6`, the Étape 1/2 reference) to
5.97 (`Lx=200e-6`, `Ly=60e-6` unchanged) — a genuine finite-size
effect from growing `Lx` alone while `Ly` stays fixed, changing the
domain from a square to an elongated strip. It is **not** a bug in
the periodic geometry or in `Clamp_x` on `xmax` (a free-`xmax` control
at the same `Lx` gives -680.0 MPa, i.e. further from 571, not closer).
**The `571.3 MPa` reference is specific to the `Lx=Ly=60e-6` square
domain from Étape 1/2 and was never checked for domain-size
convergence in `Ly`.** The 7-point `Lx` sweep below is internally
consistent (same `Ly=60e-6` throughout) but should not be read as
converging to 571.3 MPa as `Lx -> infinity`.

## Mid-height damage band in the saturated state (final frame)

At `Lx <= 60e-6`, the *final* frame (`p=-800 MPa`, up to 2.7x
`p_crit`) shows a horizontal damage band partway up the domain
instead of a clean localized crack. **This band is absent at the
physically relevant `D_max~0.9` frame** (see
`Lx_sweep_damage_panel_3col.png`): vertical extent of the `D>0.5`
region there is `0.72-1.08 um` (~1.3% of `Ly`) across the *entire*
sweep, regardless of `Lx` or how much of the ligament is already
bridged (0.5% to 38.9%) — a flat, thin band at `y~0`, not diffuse
damage.

Whether the *saturated-state* band position depends on `Ly` was
tested by rerunning `Lx=30e-6` with `Ly=120e-6` (`Lx_30_Ly120/`,
everything else identical). Band height (`y` where `D` last crosses
`0.5`) at the final frame:

- `Ly=60e-6`: `y = 16.5 um` (27.5% of `Ly`)
- `Ly=120e-6`: `y = 13.3 um` (11.1% of `Ly`)

A pure `ymax` boundary effect would put the band at a fixed *fraction*
of `Ly` (height scaling with `Ly`: ~16.5 then ~33 um). A purely
intrinsic length scale would put it at a fixed *absolute* height
(~16.5 both times). The observed `-19%` shift (16.5 -> 13.3) sits
closer to the intrinsic-length-scale hypothesis but does not match it
exactly either. **Not conclusive — noted as-is, no firm conclusion
drawn.** This is moot for the physical result either way: it concerns
an artificially overloaded state (2.7x `p_crit`) with no direct
physical relevance; the actual crack-propagation morphology
(`D_max~0.9` frame) does not show this band at all.

## Mesh-convergence check (acceptance criterion): does the nucleation
## site / propagation morphology depend on mesh resolution?

`Lx=45e-6` (middle of the isolated/coalescence transition) rerun at
three `h_cavity`: 0.30um (coarse, 2306 elements), 0.15um (reference,
5510 elements), 0.075um (fine, 15816 elements) — same `Lx`, `Ly`,
`lc=2e-6`, 201 steps, ramp to `-800e6` Pa throughout. See
`mesh_convergence_Lx45.py` / `.csv`.

| Resolution | Elements | p_crit (MPa) | Nucleation (x,y) um | D=0.9: ligament D>0.5 | D=0.9: vertical extent |
|---|---|---|---|---|---|
| h=0.30um | 2306 | -425.30 | (33.8, 0.0) -- right tip | 7.8% | 0.713 um (1.19% Ly) |
| h=0.15um | 5510 | -422.50 | (11.2, 0.0) -- left tip | 7.2% | 0.910 um (1.52% Ly) |
| h=0.075um | 15816 | -421.04 | (11.2, 0.0) -- left tip | 9.1% | 0.961 um (1.60% Ly) |

`p_crit` converges to <1% over a 6.9x range in element count
(monotone: -425.3 -> -422.5 -> -421.0). Nucleation is at a cavity tip
in all three cases (never mid-ligament) — the coarse mesh happens to
pick the *right* tip and the other two the *left* tip, which is not a
discrepancy: both tips are exactly symmetric under this loading, so
which one reaches `D=0.5` first is decided by sub-percent numerical
noise, not physics. Propagation morphology at `D_max~0.9` (thin band,
sub-micron vertical extent, ~7-9% ligament bridged) is stable across
the resolution range. **Acceptance criterion met: nucleation site and
propagation morphology are mesh-independent.**
