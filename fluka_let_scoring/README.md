# FLUKA LET scoring routines

This folder contains FLUKA user-routine code used for LET-related scoring in the SOBP/LET benchmark workflow.

## Files

- `source_sampler.f`  
  Source routine for sampling the FLUKA beam/source phase-space.

- `fluka_let_scoring.f`  
  Scoring-only user-routine file containing the `FLUSCW` and `COMSCW` routines used to calculate and filter LET, fluence, dose, and fragment-related scoring quantities.

## Purpose of `fluka_let_scoring.f`

The file provides custom FLUKA scoring weights for LET-moment and particle-filtered scoring. It is intended to be combined with a FLUKA source routine and linked into a custom FLUKA executable.

The scoring logic includes:

- all-proton and primary-proton fluence handling,
- dose-averaged and track-averaged LET moment scoring,
- primary-proton LET moment scoring using `LTRACK .EQ. 1`,
- light-fragment LET scoring for deuterons, tritons, helium-3, and helium-4 using `GETLET`,
- lithium-6 and lithium-7 filtering using heavy-fragment information from `TRACKR`/`FHEAVY`,
- primary-proton dose filtering through the `PRDO` scorer key in `COMSCW`.

## Main scoring routines

### `FLUSCW`

`FLUSCW` applies track-length and fluence scoring weights. It is used for LET-moment scoring of transported particles and fragments.

The routine uses FLUKA scorer names, accessed through `TITUSB(JSCRNG)`, to decide which quantity is being scored. The scorer name is normalized with `TRIM(ADJUSTL(...))` before comparison.

Examples of scorer keys handled in `FLUSCW` include:

- `PHL1`, `PHL2` for proton LET moments,
- `PWL1`, `PWL2` for water-equivalent/local LET moments,
- `PRI_`, `PRL1`, `PRL2`, `PWR1`, `PWR2` for primary-proton fluence and LET moments,
- fragment scorer keys for deuterons, tritons, helium-3, helium-4, lithium-6, and lithium-7.

### `COMSCW`

`COMSCW` applies dose/energy-deposition filters.

The implemented dose filters include:

- `PRDO`: primary-proton dose only, requiring `IJ .EQ. 1` and `LTRACK .EQ. 1`,
- `LI6_`: lithium-6 dose only,
- `LI7_`: lithium-7 dose only.

The `PRDO` filter is different from a standard proton-only `AUXSCORE PROTON` dose scorer. A proton-only `AUXSCORE` includes both primary and secondary protons, whereas `PRDO` keeps only source-generation protons.

## Important limitations

The fragment filters identify transported fragment species, but they do not identify the production vertex, parent particle, target nucleus, or reaction channel. Those ancestry details would require additional FLUKA stack or history bookkeeping, such as through `STUPRF` or related tracking logic.

Lithium isotope LET is not obtained through `GETLET` in this implementation. Instead, lithium LET is reconstructed from transported heavy-fragment energy deposition and track-length information.

## Use in FLUKA

This file is not a complete standalone FLUKA executable. It provides the scoring routines and must be compiled/linked with the relevant FLUKA user-routine source code, including an appropriate source routine, to produce a custom FLUKA executable.
