# Z3ST ← OpenMC coupling

One-way, neutronics → Z3ST, exchanged as a file. **Z3ST never imports OpenMC**:
the scripts here run in the OpenMC environment, their output is a YAML fragment
for a Z3ST material card, and the case then runs with no OpenMC dependency.

## What is coupled

| Quantity | Status | Bus |
|---|---|---|
| Axial power form factor f(z) | available | `materials.fuel_profiles.tabulated_axial` |
| Linear heat rate (W/m) | not extracted, see below | case input `lhr:` |
| Radial form factor f(r) | not extracted, see below | `radial_profile:` |
| Fast fluence / dpa | not implemented | — |

## axial_power.py

```bash
conda activate openmc
python axial_power.py statepoint.300.h5 \
    --pin-x 0 --pin-y 8.2 --z-active 8.81 46.91 -o fuel_axial.yaml
```

Reads a Cartesian mesh tally, sums the columns at the element position, keeps
the rows inside the active fuel, shifts z to the bottom of active fuel, converts
cm → m and normalises the shape to mean 1. The header of the emitted file
records the statepoint, the tally, k-effective, the axial peaking factor and the
worst per-row statistical uncertainty — enough to trace a table back to the run
that produced it.

Paste the fragment into the fuel card. `set_power` normalises the composite form
factor to nodal mean 1, so only the shape is used.

## Two things this does not give you

**The absolute rating.** A Cartesian core mesh tally has no per-element power
fraction: the cells do not follow the element boundaries and a flux tally is not
a power tally. The rating stays the case input `lhr:`, set from the operating
power and the element peaking factor. To take it from OpenMC instead, tally
`fission-q-prompt` over the fuel cell of the element and normalise with the
core-wide `heating-local`.

**The radial shape inside the element.** The same mesh resolves the whole core in
cells wider than one element, so it cannot resolve flux depression across the
fuel meat. That needs a cylindrical mesh tally, or a cell tally on radial rings,
on the OpenMC side. Until then f(r) ≡ 1.

Both are limitations of the tally, not of the reader: they are fixed in the
OpenMC model, and this script reads whatever mesh tally it is given.
