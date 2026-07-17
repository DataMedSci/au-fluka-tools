# FLUKA LET scoring routines

FLUKA `FLUSCW` and `COMSCW` user-weighting routines for scoring LET moments, from
which averaged LET quantities (track-averaged and dose-averaged LET) are reconstructed
in post-processing.

These are **user routines**, not a standalone program. You link them into a custom FLUKA
executable, then activate them from your input file with a `USERWEIG` card.

Reference: Kalholm F, Grzanka L, Traneus E, Bassler N. *A systematic review on the usage
of averaged LET in radiation biology for particle therapy.* Radiotherapy and Oncology.
2021;161:211-21.

## Quick start

1. In your FLUKA input, give a `USRBIN` scorer a name whose **first four characters** are
   a key this routine recognises (e.g. `ALL1`). The name is read from `TITUSB(JSCRNG)`;
   only the first four characters select the branch, so you may append a suffix for your
   own bookkeeping (e.g. `ALL1_ZN`).
2. Activate user weighting with a `USERWEIG` card (`WHAT(3)=1` to call `FLUSCW`,
   `WHAT(6)=1` to call `COMSCW`).
3. Compile and link the routine into a custom FLUKA executable (see *Compilation*).
4. Run FLUKA with that executable.
5. Post-process the scored moments into averaged LET (see the example below).

## The main quantity: all-particle dose-averaged LET

The headline scorers are `ALL1` and `ALL2`: the first and second LET moments over **all
charged hadrons and ions**. Their LET is the true *local* energy-deposition LET,
reconstructed per step from `TRACKR` as `LET = 100 * ΣDTRACK / ΣTTRACK` (dE_dep /
step-length, GeV/cm → keV/µm). This works for every charged hadron and ion FLUKA
transports — protons, light ions, and heavy fragments — with no per-species handling and
no dependence on `GETLET` or an explicit material-index lookup. The LET still depends
physically on the local material, through the energy deposited per unit track length.

Neutral particles are skipped, and **electrons and positrons are excluded** even when EMF
transport is active: averaged-LET reporting conventionally covers the hadron/ion field,
and including deltas would tie track-averaged LET to the EMF transport threshold rather
than to the physics.

Score three co-located `USRBIN` bins over the same region: `ALL1`, `ALL2`, and `ALFL`.
Schematic input:

```text
* userweig: call FLUSCW (WHAT(3)=1)
USERWEIG          0.0       0.0       1.0                              &
*
* first LET moment   (weight = LET)             -> ALL1  (unit 22)
USRBIN           11.0  ALL-PART      -22.  <xmax ymax zmax bins...>    ALL1
USRBIN         <xmin ymin zmin> ...                                   &
* second LET moment  (weight = LET^2)           -> ALL2  (unit 23)
USRBIN           11.0  ALL-PART      -23.  <xmax ymax zmax bins...>    ALL2
USRBIN         <xmin ymin zmin> ...                                   &
* unweighted fluence, same particle set         -> ALFL  (unit 21)
USRBIN           11.0  ALL-PART      -21.  <xmax ymax zmax bins...>    ALFL
USRBIN         <xmin ymin zmin> ...                                   &
```

All three bins use the `ALL-PART` generalized particle; `FLUSCW` does the actual particle
selection, and multiplies each track segment by the returned weight (LET for `ALL1`, LET²
for `ALL2`, 1 for `ALFL`). Then, bin by bin:

- **dose-averaged LET**  `LETd = ALL2 / ALL1`
- **track-averaged LET** `LETt = ALL1 / ALFL`

(`ALL2/ALL1` needs no separate fluence bin, because a segment's dose contribution is
∝ length·LET, so the dose-weighted mean of LET is `Σ(ℓ·LET·LET)/Σ(ℓ·LET)`.)

> **Use `ALFL`, not a plain unweighted `ALL-PART` bin, as the track-average denominator.**
> A plain `ALL-PART` bin gets no `FLUSCW` filtering and so counts neutrons, photons and
> electrons — none of which contribute to the `ALL1` numerator. Dividing by it
> underestimates `LETt` (~3% at entrance in `tests/`, and unboundedly past the distal
> edge, where the neutral fluence is non-zero while `ALL1` is exactly zero). `ALFL`
> applies the same particle selection as `ALL1`/`ALL2`, so the ratio is consistent by
> construction. `LETd = ALL2/ALL1` is unaffected either way.

## Two LET definitions

This routine offers two different, complementary LET definitions — the distinction the
Kalholm review stresses:

- **Local energy-deposition LET** (`ALL1/ALL2`, and the lithium scorers): the actual
  dE/dx deposited along the step in the current material, from `TRACKR`. Universal (all
  charged hadrons and ions), always evaluated in the **local** material.
- **Electronic stopping-power LET** (`GETLET`): analytic mass stopping power converted to
  linear LET via `LETLIN = RHO(MAT) * GETLET(...)`. Can be evaluated in the **local**
  material or in **water**, but only for the light particles `GETLET` supports
  (p, d, t, ³He, ⁴He).

## Material convention

- **Local scorers** use `MATLET = MEDFLK(NREG,1)` — LET in whatever material the particle
  is currently in. The only guard is skipping vacuum / non-material regions
  (`MATLET ≤ 0` or `RHO ≤ 0`). Restricting *where* you score is the job of the `USRBIN`
  geometry, not of a material list. (The old hardcoded `27/28/29/30` filter has been
  removed.)
- **Water-reference scorers** use `MATLET = MWATER`, resolved **automatically** on the
  first scoring call — there are no hardcoded material numbers. `MWATER` is taken from
  FLUKA's built-in `MATQLT` (the "extra water material for Q(L) calculations" in
  `flkmat.inc`, present even when the input defines no explicit `WATER`); failing that,
  from the first material named `WATER` in `MATNAM`. The chosen index is written to the
  FLUKA output as `fluka_let_scoring: water MWATER = <n>` so you can verify it. If no
  water material can be found, the water-reference scorers return zero and a warning is
  printed.

## Moment convention

`L1`/`L2` (and `W1`/`W2` for water) suffixes mean:

- `1`: first raw LET moment — weighted by LET. Units keV/µm.
- `2`: second raw LET moment — weighted by LET². Units (keV/µm)².

## Scorer keys handled by `FLUSCW`

All-particle (recommended starting point):

| Key | Meaning | Selection | Material | Weight |
|---|---|---|---|---|
| `ALL1` | all-particle LET moment | charged hadrons + ions (no e±) | local | LET |
| `ALL2` | all-particle LET² moment | charged hadrons + ions (no e±) | local | LET² |
| `ALFL` | all-particle fluence (`ALL1`/`ALL2` denominator) | charged hadrons + ions (no e±) | — | 1 |
| `ALW1` | all-particle water LET moment | p, d, t, ³He, ⁴He | water | LET |
| `ALW2` | all-particle water LET² moment | p, d, t, ³He, ⁴He | water | LET² |
| `ALWF` | water fluence (`ALW1`/`ALW2` denominator) | p, d, t, ³He, ⁴He | — | 1 |

Protons:

| Key | Meaning | Selection | Material | Weight |
|---|---|---|---|---|
| `PAL1`/`PAL2` | all-proton LET / LET² | all protons | local | LET / LET² |
| `PAW1`/`PAW2` | all-proton LET / LET² | all protons | water | LET / LET² |
| `P1FL` | primary-proton fluence filter | protons, `LTRACK=1` | — | 1 or 0 |
| `P1L1`/`P1L2` | primary-proton LET / LET² | protons, `LTRACK=1` | local | LET / LET² |
| `P1W1`/`P1W2` | primary-proton LET / LET² | protons, `LTRACK=1` | water | LET / LET² |

Light fragments (via `GETLET`, `IJ` = FLUKA particle id):

| Key | Particle | `IJ` | Weight |
|---|---|---:|---|
| `D2L1`/`D2L2` | deuteron | -3 | LET / LET² |
| `T3L1`/`T3L2` | triton | -4 | LET / LET² |
| `H3L1`/`H3L2` | helium-3 | -5 | LET / LET² |
| `H4L1`/`H4L2` | helium-4 / α | -6 | LET / LET² |

Lithium (via `TRACKR`, because `GETLET` returns zero for transported Li):

| Key | Isotope | Weight |
|---|---|---|
| `L6L1`/`L6L2` | Li-6 (Z=3, A=6) | LET / LET² |
| `L7L1`/`L7L2` | Li-7 (Z=3, A=7) | LET / LET² |
| `L6FL`/`L7FL` | Li-6 / Li-7 fluence filter | 1 or 0 |

## Scorer keys handled by `COMSCW`

`COMSCW` weights dose-like (energy-deposition) estimators, where `FLUSCW` does not apply.

| Key | Meaning | Selection | Weight |
|---|---|---|---|
| `ALDD` | dirty dose (LET > 30 MeV cm²/g) | charged hadrons + ions (no e±) | 1 or 0 |
| `P1DO` | primary-proton dose filter | `JTRACK=1`, `LTRACK=1` | 1 or 0 |
| `L6DO` | Li-6 dose filter | Z=3, A=6 | 1 or 0 |
| `L7DO` | Li-7 dose filter | Z=3, A=7 | 1 or 0 |

`P1DO` keeps only source-generation protons — unlike a standard `AUXSCORE PROTON` dose
scorer, which includes secondary protons too.

> Identify the particle in `COMSCW` by `JTRACK`, **not** by the `IJ` argument. Unlike in
> `FLUSCW`, `COMSCW`'s `IJ` is the *generalized quantity being deposited* (208 = `ENERGY`
> for dose scoring), not the particle type — so `IJ .EQ. 1` is never true and silently
> scores zero everywhere. `COMSCW` is also called for point-like depositions whose
> `JTRACK` is a pseudo-particle id above the normal range (208 heavy recoil, 211 e/γ below
> threshold, 308 low-energy neutron kerma), which will run off the end of `ICHRGE(-6:64)`
> if not screened first.

## Dirty dose

Dirty dose is the dose deposited by particles whose LET exceeds a threshold — here
**30 MeV cm²/g** of unrestricted mass stopping power, equivalently **3 keV/µm in water**.
It answers "how much of this dose was delivered by high-LET particles?", and pairs
naturally with the LET scorers above.

Score it with `ALDD` plus an unfiltered bin of the same generalized particle:

```text
* dirty dose, LET judged in the local medium        (unit 24)
USRBIN           10.0      DOSE      -24.  <xmax ymax zmax bins...>    ALDD
USRBIN         <xmin ymin zmin> ...                                   &
* unfiltered dose over the same region              (unit 25)
USRBIN           10.0      DOSE      -25.  <xmax ymax zmax bins...>    DOSE
USRBIN         <xmin ymin zmin> ...                                   &
```

then `dirty fraction = ALDD / DOSE`, bin by bin.

**The material the LET threshold is judged in is not a separate key.** There are two axes
— the dose being scored, and the material the LET is judged in — but only the two matching
combinations mean anything, so `ALDD` reads the binning's own generalized particle
(`IDUSBN`) and follows it:

| `USRBIN` WHAT(2) | LET judged in | Route | Covers |
|---|---|---|---|
| `DOSE` (228) | local medium | `TRACKR` | all charged hadrons + ions, incl. heavy fragments |
| `DOSE-H2O` (252) | water | `GETLET` | p, d, t, ³He, ⁴He only |

One key serves both, and the meaningless cross combinations (dose-to-water thresholded on
medium LET, or vice versa) cannot be expressed. `ALDD` on any other binning is rejected
with a warning rather than scoring something arbitrary.

Caveats worth knowing before quoting a number:

- **`DOSE-H2O` misses heavy fragments.** `GETLET` has no water stopping power for Li and
  above, so they drop out. That bites harder here than for the LET moments, because
  fragments are exactly the high-LET component dirty dose is meant to capture. The `DOSE`
  variant has no such gap.
- **Point-like depositions are excluded**, because `TRACKR` carries no step data for them
  and their LET cannot be reconstructed. Two of these are genuinely high-LET and so are
  under-counted: heavy recoils (`JTRACK` 208) and low-energy neutron kerma (308).
- **Electrons and positrons are excluded**, as for `ALL1`/`ALL2`. Low-energy electrons can
  exceed the threshold, so this is a real choice: it keeps dirty dose a property of the
  hadron/ion field.

> **`DOSE-H2O` requires a `RAD-BIOL` card**, or FLUKA dies with a SIGFPE in
> `dedx/alphbt.f:382` (`alphbt` called with `rbealp=-1`, `rbebet=0`) via `score/usrsco.f`
> — no message about the real cause, and only in materials other than water, so a water
> phantom runs and a carbon one aborts. This is FLUKA's own behaviour: it reproduces with
> the stock `fluka` executable and no user routines at all. The α/β file content is *not*
> used by `DOSE-H2O` — a FLUKA developer
> [confirmed on the forum](https://fluka-forum.web.cern.ch/t/dose-to-water-issue/6128) that
> it is required only "for technical reasons — to be overcome at a later stage". Any valid
> file will do; `tests/letbio.dat` is a minimal one.

## A note on ancestry

`LTRACK` gives the generation (1 = source-generation), and the isotope filters classify
the *currently transported* particle. None of these record the production vertex, parent
particle, or reaction channel — that would require production-time tagging via `STUPRF`
or `MDSTCK`.

## Compilation

Compile the routine and link it into a custom FLUKA executable with the FLUKA build
tools. Compiled and linked with **FLUKA 4** (`fff` + `lfluka`):

```bash
export FLUPRO=/usr/local/fluka        # your FLUKA install
export PATH=$PATH:$FLUPRO/bin
fff fluka_let_scoring.f               # -> fluka_let_scoring.o
lfluka -m fluka -o flukalet fluka_let_scoring.o
```

The include files are referenced by their literal names (`INCLUDE 'dblprc.inc'`), which
is what `fff` expects. The exact tool names vary between FLUKA distributions (older
`ldpmqmd` wrappers work too).

## Testing

`tests/` contains a ready-to-run example, `plan01_field01_geoA_SOBPcent.inp`, that
exercises the scorer keys (proton, light-fragment, and lithium LET moments plus the
dose/fluence filters). It uses a `SOURCE` reading the accompanying `sobp.dat` spot list,
so the executable must also link the source sampler from `../fluka_sobp_source`:

```bash
fff ../fluka_sobp_source/source_sampler.f
lfluka -m fluka -o flukalet fluka_let_scoring.o source_sampler.o
cd tests && rfluka -e ../flukalet -N0 -M1 plan01_field01_geoA_SOBPcent
```

On the first scoring call the run prints `fluka_let_scoring: water MWATER = <n>` to the
`.out`, confirming the auto-detected water material.

## Adding or changing a scorer

1. Keep the key to four characters and add it to the right dispatch block
   (`FLUSCW` for fluence/track-length, `COMSCW` for dose).
2. Document: particle selection; all vs. primary-only; local vs. water; LET, LET², or
   filter weight.
3. Update the input cards and any post-processing that reads the output units.
4. Compile-test, then run a small smoke test before production use.
