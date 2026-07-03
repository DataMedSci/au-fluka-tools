# FLUKA LET scoring routines

## Required companion source routine

`fluka_let_scoring.f` is not a complete FLUKA executable by itself. It provides the `FLUSCW` and `COMSCW` scoring routines only.

To run the SOBP/LET benchmark input files in this repository, it must be compiled and linked together with the SOBP source routine:

`../fluka_sobp_source/source_sampler.f`

That source routine implements the FLUKA `SOURCE` entry point and reads the beam/source table referenced by the input cards, for example `sobp.dat`. Without a matching `SOURCE` routine, the benchmark input files cannot sample the intended proton field.


This folder contains FLUKA user-routine code for LET-moment, fluence-filter, and dose-filter scoring.

The main file is:

- `fluka_let_scoring.f`

It implements two FLUKA user routines:

- `FLUSCW`: user weighting for fluence-like and track-length-like estimators.
- `COMSCW`: user weighting for dose-like energy-deposition estimators.

These routines are intended for FLUKA simulations where LET-related quantities are reconstructed from scorer moments.

## Reference

This implementation follows the need for explicit LET definitions emphasized in:

Kalholm F, Grzanka L, Traneus E, Bassler N. A systematic review on the usage of averaged LET in radiation biology for particle therapy. Radiotherapy and Oncology. 2021 Aug 1;161:211-21.

When using or modifying these routines, always document:

- which particles are included,
- whether primary-only or all transported particles are scored,
- whether LET is evaluated in the local material or in water,
- whether the score is a first LET moment or a second LET-squared moment,
- how the final averaged LET quantity is reconstructed in post-processing.

## What this file is and is not

`fluka_let_scoring.f` is not a complete FLUKA input and not a complete executable by itself.

It provides scoring user routines that must be linked into a custom FLUKA executable together with the source routine used by the simulation.

A typical workflow is:

1. Write or prepare a normal FLUKA input file with `USRBIN` scorers.
2. Give selected `USRBIN` scorers one of the four-character names listed below.
3. Attach `FLUSCW` or `COMSCW` through `USERWEIG` and, where relevant, `AUXSCORE`.
4. Compile and link `fluka_let_scoring.f` with FLUKA.
5. Run FLUKA with the custom executable.
6. Post-process the scored moments into averaged LET quantities.

## Four-character scorer-key convention

The scorer identifiers are kept to four characters because the USRBIN/AUXSCORE workflow used here relies on four-character scorer keys.

The routine reads the active scorer name through `TITUSB(JSCRNG)` and stores it as `SCONAM`.

The first four characters of `SCONAM` decide which branch of `FLUSCW` or `COMSCW` is applied.

## Moment convention

The `L1` and `L2` suffixes mean:

- `L1`: first raw LET moment, meaning the contribution is weighted by LET.
- `L2`: second raw LET moment, meaning the contribution is weighted by LET².

Here, LET means linear energy transfer based on the electronic stopping power of the selected material. For local-material scorers, the selected material is the current FLUKA transport material. For water-reference scorers, the selected material is the configured water-equivalent material.

The returned units are:

- `L1`: keV/µm
- `L2`: (keV/µm)²

In post-processing, these moments can be combined with the corresponding unweighted estimator to reconstruct averaged LET quantities.

For example, in a track-length-like LET reconstruction:

- track-averaged LET is reconstructed from a first LET moment divided by the matching unweighted fluence or track-length estimator.

For a dose-like LET reconstruction using first and second LET moments:

- dose-averaged LET is reconstructed from the second LET moment divided by the first LET moment.

The exact normalization depends on the estimator and post-processing script, so the scorer pairings must be kept consistent.

## `GETLET` and electronic stopping power

`GETLET` is the central FLUKA helper used by this routine for protons, deuterons, tritons, helium-3, and helium-4.

In this implementation, `GETLET` is called as:

    LETW = GETLET(IJ, EKIN, PLA, ZERZER, MATLET)

where:

- `IJ` is the FLUKA particle identifier.
- `EKIN = -PLA` is the kinetic energy used for the stopping-power lookup.
- `MATLET` is the material in which the LET/stopping power is evaluated.

The routine treats the value returned by `GETLET` as a mass electronic stopping power. It is converted to linear LET with:

    LETLIN = RHO(MATLET) * LETW

This gives LET in keV/µm for the selected material.

This material choice is important:

- local-material scorer keys use `MATLET = MEDFLK(NREG,1)`;
- water-reference scorer keys use `MATLET = 30`.

Thus, the same transported particle can be scored with LET evaluated either in the local transport material or in the selected water-equivalent material.

Li-6 and Li-7 are exceptions in this implementation. The tested `GETLET` calls returned zero for transported lithium ions, so lithium LET is reconstructed from `TRACKR` energy-deposition and track-length quantities instead.

## Important FLUKA variables used

The routine uses FLUKA common-block variables provided by the included FLUKA header files.

Important variables include:

- `IJ`: FLUKA particle identifier passed to `FLUSCW` and `COMSCW`.
  - `IJ = 1`: proton
  - `IJ = -3`: deuteron
  - `IJ = -4`: triton
  - `IJ = -5`: helium-3
  - `IJ = -6`: helium-4 / alpha
- `LTRACK`: generation index for the currently transported particle.
  - `LTRACK = 1`: source-generation particle.
  - `LTRACK > 1`: secondary or later-generation particle.
- `MEDFLK(NREG,1)`: material index of the current FLUKA region.
- `RHO(MATLET)`: density of the material used for converting mass stopping power to linear LET.
- `TITUSB(JSCRNG)`: scorer name for the active scoring call.
- `TRACKR` / `FHEAVY`: heavy-fragment transport bookkeeping used here for Li-6 and Li-7.

## Material convention

For local-material LET scoring, the routine uses:

`MATLET = MEDFLK(NREG,1)`

This means LET is evaluated in the material of the current FLUKA region.

The current implementation restricts local-material LET scoring to material indices:

- `27`
- `28`
- `29`
- `30`

These numbers are geometry-specific material indices for the benchmark phantom/slab scoring media. If the material card order or geometry changes, these numbers must be checked and updated.

For water-reference LET scoring, the routine uses:

`MATLET = 30`

This assumes that material index `30` is the intended water-equivalent material in the benchmark geometry. This must also be checked if the material definitions are changed.

## Proton scorer keys handled by FLUSCW

These scorers are handled in `FLUSCW`.

| Key | Meaning | Particle selection | Material convention | Returned weight |
|---|---|---|---|---|
| `PAL1` | all-proton local-material LET moment | all transported protons | local material | LET |
| `PAL2` | all-proton local-material LET² moment | all transported protons | local material | LET² |
| `PAW1` | all-proton water-reference LET moment | all transported protons | water-reference material | LET |
| `PAW2` | all-proton water-reference LET² moment | all transported protons | water-reference material | LET² |
| `P1FL` | primary-proton fluence filter | protons with `LTRACK = 1` | not applicable | 1 or 0 |
| `P1L1` | primary-proton local-material LET moment | protons with `LTRACK = 1` | local material | LET |
| `P1L2` | primary-proton local-material LET² moment | protons with `LTRACK = 1` | local material | LET² |
| `P1W1` | primary-proton water-reference LET moment | protons with `LTRACK = 1` | water-reference material | LET |
| `P1W2` | primary-proton water-reference LET² moment | protons with `LTRACK = 1` | water-reference material | LET² |

“Primary proton” here means source-generation proton, implemented with `LTRACK .EQ. 1`.

“All protons” includes both source-generation protons and secondary or later-generation protons.

## Light-fragment scorer keys handled by FLUSCW

These scorers are handled in `FLUSCW`.

| Key | Meaning | FLUKA particle identifier | Returned weight |
|---|---|---:|---|
| `D2L1` | deuteron LET moment | `IJ = -3` | LET |
| `D2L2` | deuteron LET² moment | `IJ = -3` | LET² |
| `T3L1` | triton LET moment | `IJ = -4` | LET |
| `T3L2` | triton LET² moment | `IJ = -4` | LET² |
| `H3L1` | helium-3 LET moment | `IJ = -5` | LET |
| `H3L2` | helium-3 LET² moment | `IJ = -5` | LET² |
| `H4L1` | helium-4 / alpha LET moment | `IJ = -6` | LET |
| `H4L2` | helium-4 / alpha LET² moment | `IJ = -6` | LET² |

These light-fragment scorers use `GETLET`.

The routine obtains a mass stopping power from `GETLET`, then converts it to linear LET using:

`LETLIN = RHO(MATLET) * LETW`

## Lithium scorer keys handled by FLUSCW

Lithium isotope LET is not obtained through `GETLET` in this implementation. The tested `GETLET` calls returned zero for transported lithium ions, so Li-6 and Li-7 are handled through FLUKA heavy-fragment bookkeeping.

For lithium LET moments, the routine reconstructs LET from `TRACKR` quantities:

`LET = 100 * SUMD / SUMT`

where:

- `SUMD` is the summed energy deposition along the transported fragment step in GeV,
- `SUMT` is the summed track length in cm,
- the factor `100` converts GeV/cm to keV/µm.

| Key | Meaning | Isotope selection | Returned weight |
|---|---|---|---|
| `L6L1` | lithium-6 LET moment | Z = 3, A = 6 | LET |
| `L6L2` | lithium-6 LET² moment | Z = 3, A = 6 | LET² |
| `L7L1` | lithium-7 LET moment | Z = 3, A = 7 | LET |
| `L7L2` | lithium-7 LET² moment | Z = 3, A = 7 | LET² |
| `L6FL` | lithium-6 fluence filter | Z = 3, A = 6 | 1 or 0 |
| `L7FL` | lithium-7 fluence filter | Z = 3, A = 7 | 1 or 0 |

## Dose-filter scorer keys handled by COMSCW

These scorers are handled in `COMSCW`.

`COMSCW` is needed for dose-like energy-deposition estimators, because `FLUSCW` is not applied to dose scoring in the same way.

| Key | Meaning | Selection | Returned weight |
|---|---|---|---|
| `P1DO` | primary-proton dose filter | `IJ = 1` and `LTRACK = 1` | 1 or 0 |
| `L6DO` | lithium-6 dose filter | Z = 3, A = 6 | 1 or 0 |
| `L7DO` | lithium-7 dose filter | Z = 3, A = 7 | 1 or 0 |

`P1DO` is different from a standard `AUXSCORE PROTON` dose scorer.

A standard proton-only `AUXSCORE PROTON` includes both primary and secondary protons. `P1DO` keeps only source-generation protons.

## Minimal FLUSCW example

`FLUSCW` is used for fluence-like or track-length-like estimators. In practice, the FLUKA input must contain a `USRBIN` scorer with a four-character scorer key that this routine recognizes.

Schematic example for an all-proton local-material LET moment:

    USRBIN        ...    ENERGY      ...       PAL1
    USRBIN        ...       ...      ...

In this schematic example, the scorer key is `PAL1`. When FLUKA calls `FLUSCW` for that scorer, this routine:

1. checks that the transported particle is a proton;
2. uses the local material with `MATLET = MEDFLK(NREG,1)`;
3. calls `GETLET` for the proton electronic stopping power in that material;
4. converts the result to linear LET with `LETLIN = RHO(MATLET) * LETW`;
5. returns `LETLIN` as the scoring weight.

To reconstruct averaged LET, the `L1` and `L2` scorers must be paired with the corresponding unweighted scorer in post-processing.

Typical all-proton local-material set:

- unweighted all-proton fluence or track-length scorer;
- `PAL1`, all-proton local-material LET moment;
- `PAL2`, all-proton local-material LET² moment.

Typical primary-proton local-material set:

- `P1FL`, primary-proton fluence filter;
- `P1L1`, primary-proton local-material LET moment;
- `P1L2`, primary-proton local-material LET² moment.

Typical water-reference set:

- `PAW1` and `PAW2` for all transported protons evaluated in water;
- `P1W1` and `P1W2` for primary protons evaluated in water.

The exact `USRBIN` geometry, binning, estimator type, unit number, and output filename are defined in the FLUKA input file, not in this user routine.

## Minimal COMSCW example

`COMSCW` is used for dose-like energy-deposition estimators. It is needed when a `DOSE`-type scorer must reject some particles before the deposited-energy contribution is scored.

Schematic example for primary-proton dose:

    USRBIN        ...      DOSE      ...       P1DO
    USRBIN        ...       ...      ...

In this schematic example, the scorer key is `P1DO`. When FLUKA calls `COMSCW` for that scorer, this routine:

1. checks that the active call is for a dose-like scoring estimator;
2. checks that `IJ = 1`, meaning the transported particle is a proton;
3. checks that `LTRACK = 1`, meaning the proton is source-generation/primary;
4. returns `ONEONE` if both conditions are true;
5. returns `ZERZER` otherwise.

Therefore `P1DO` is a primary-proton dose filter.

This differs from a standard proton-only `AUXSCORE PROTON` dose scorer. A standard proton-only dose scorer includes both primary and secondary protons. `P1DO` keeps only source-generation protons.

The lithium dose filters work similarly:

- `L6DO` keeps only dose contributions from transported Li-6 fragments;
- `L7DO` keeps only dose contributions from transported Li-7 fragments.

## Compilation

This file must be compiled and linked with FLUKA's user-routine build tools.

Example:

    export FLUPRO=/path/to/fluka
    export FLUKADATA=$FLUPRO/data
    $FLUPRO/bin/ldpmqmd -o fluka_let_scoring_test fluka_let_scoring.f

The exact FLUKA path and wrapper name may differ between installations.

This repository has been compile-tested with literal FLUKA include-file names such as:

    INCLUDE 'dblprc.inc'

On the tested FLUKA installation, token-style includes such as `INCLUDE '(DBLPRC)'` did not compile with the local `ldpmqmd` wrapper.

## Practical checklist when adding or changing a scorer

When adding or changing a scorer key:

1. Keep the scorer key to four characters.
2. Add the key to the correct dispatch block in `FLUSCW` or `COMSCW`.
3. Document whether the scorer is fluence-like, track-length-like, or dose-like.
4. Document the particle selection.
5. Document whether the scorer includes all transported particles or only primary/source-generation particles.
6. Document whether LET is evaluated in the local material or in the water-reference material.
7. Document whether the returned quantity is LET, LET², or a filter weight.
8. Update any FLUKA input cards using the scorer key.
9. Update post-processing scripts that read the corresponding output files.
10. Compile-test the modified user routine.
11. Run at least one small FLUKA smoke test before using the scorer in production.
