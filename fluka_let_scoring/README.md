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
charged particles**. Their LET is the true *local* energy-deposition LET, reconstructed
per step from `TRACKR` as `LET = 100 * ΣDTRACK / ΣTTRACK` (dE_dep / step-length,
GeV/cm → keV/µm). This works for every charged particle FLUKA transports — protons,
light ions, and heavy fragments — with no per-species handling and no dependence on
`GETLET`.

Score three co-located `USRBIN` bins over the same region: an ordinary (unweighted)
track-length fluence, plus `ALL1` and `ALL2`. Schematic input:

```text
* userweig: call FLUSCW (WHAT(3)=1)
USERWEIG          0.0       0.0       1.0                              &
*
* unweighted all-particle track-length fluence  -> Phi   (unit 21)
USRBIN           11.0  ALL-PART      -21.  <xmax ymax zmax bins...>    PHI
USRBIN         <xmin ymin zmin> ...                                   &
* first LET moment   (weight = LET)             -> ALL1  (unit 22)
USRBIN           11.0  ALL-PART      -22.  <xmax ymax zmax bins...>    ALL1
USRBIN         <xmin ymin zmin> ...                                   &
* second LET moment  (weight = LET^2)           -> ALL2  (unit 23)
USRBIN           11.0  ALL-PART      -23.  <xmax ymax zmax bins...>    ALL2
USRBIN         <xmin ymin zmin> ...                                   &
```

`FLUSCW` multiplies each track segment by the returned weight (1 for the unweighted bin,
LET for `ALL1`, LET² for `ALL2`). Then, bin by bin:

- **dose-averaged LET**  `LETd = ALL2 / ALL1`
- **track-averaged LET** `LETt = ALL1 / Phi`

(`ALL2/ALL1` needs no separate fluence bin, because a segment's dose contribution is
∝ length·LET, so the dose-weighted mean of LET is `Σ(ℓ·LET·LET)/Σ(ℓ·LET)`.)

## Two LET definitions

This routine offers two different, complementary LET definitions — the distinction the
Kalholm review stresses:

- **Local energy-deposition LET** (`ALL1/ALL2`, and the lithium scorers): the actual
  dE/dx deposited along the step in the current material, from `TRACKR`. Universal (all
  charged particles), always evaluated in the **local** material.
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
| `ALL1` | all-particle LET moment | all charged particles | local | LET |
| `ALL2` | all-particle LET² moment | all charged particles | local | LET² |
| `ALW1` | all-particle water LET moment | p, d, t, ³He, ⁴He | water | LET |
| `ALW2` | all-particle water LET² moment | p, d, t, ³He, ⁴He | water | LET² |

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
| `P1DO` | primary-proton dose filter | `IJ=1`, `LTRACK=1` | 1 or 0 |
| `L6DO` | Li-6 dose filter | Z=3, A=6 | 1 or 0 |
| `L7DO` | Li-7 dose filter | Z=3, A=7 | 1 or 0 |

`P1DO` keeps only source-generation protons — unlike a standard `AUXSCORE PROTON` dose
scorer, which includes secondary protons too.

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
