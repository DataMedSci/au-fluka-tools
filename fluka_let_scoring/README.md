# LET, Qeff and Dirty Dose scoring routines

FLUKA `FLUSCW` and `COMSCW` user-weighting routines that score the LET and Qeff moments
from which track-averaged and dose-averaged radiation-quality quantities are reconstructed
in post-processing.

These are **user routines**, not a standalone program: you link them into a custom FLUKA
executable (see [Compilation](#compilation)) and activate them from your input with a
`USERWEIG` card.

## Quick start

1. Add co-located `USRBIN` bins whose **first four characters** are a key from the recipe
   table below (e.g. `ALW1` for LET to water). Only the first four characters select the
   branch, so a suffix for your own bookkeeping is fine (`ALW1_ZN`).
2. Activate user weighting with a `USERWEIG` card: `WHAT(3)=1` calls `FLUSCW`, `WHAT(6)=1`
   calls `COMSCW`.
3. Compile and link into a custom executable ([Compilation](#compilation)), then run FLUKA
   with it.
4. Divide the scored bins as the recipe table says ([Post-processing](#post-processing)).

## What you can score

Score LET **to water** (`ALW*`) unless you have a specific reason not to — water is the
reference medium RBE models are built on, and scoring the local medium by mistake is the
single most common error here. The one catch: the water route covers p, d, t, ³He and ⁴He
only, so when you need the high-LET tail from Li and heavier fragments, score the local
medium (`ALL*`) instead — it covers every charged hadron and ion FLUKA transports.

Score the bins named in each ratio, **co-located over the same region**, run, then divide
bin by bin. Every key scores charged hadrons and ions; e± and neutrals are excluded (see
[Implementation details](#implementation-details)).

**LET** — to water by default:

| Quantity | To water (default) | In local medium (adds heavy fragments) |
|---|---|---|
| Dose-averaged LET (LETd) | `ALW2 / ALW1` | `ALL2 / ALL1` |
| Track-averaged LET (LETt) | `ALW1 / ALWF` | `ALL1 / ALFL` |

**Qeff** — material-independent; the same keys work in any medium:

| Quantity | Score & divide |
|---|---|
| Dose-averaged Qeff | `ALQD / DOSE` |
| Track-averaged Qeff | `ALQ1 / ALQF` |

**Dirty dose** — the dose fraction above a radiation-quality threshold:

| Threshold | Score & divide | Notes |
|---|---|---|
| LET > 30 MeV cm²/g, in water | `ALDD` on `DOSE-H2O` / `DOSE-H2O` | light ions only (p/d/t/³He/⁴He); needs a `RAD-BIOL` card |
| LET > 30 MeV cm²/g, in medium | `ALDD` on `DOSE` / `DOSE` | covers heavy fragments |
| Qeff > 30 | `ALDQ` / `DOSE` | material-independent |

Two rules keep the ratios meaningful:

- **Use the matching fluence key** as a track-average denominator — `ALWF` for water,
  `ALFL` for medium, `ALQF` for Qeff — *not* a plain `ALL-PART` bin. The filtered key
  applies the same particle selection as its numerator; a plain bin also counts neutrons,
  photons and electrons that never enter the numerator, underestimating the track average —
  and diverging past the distal edge, where neutral fluence is non-zero but the numerator is
  exactly zero. The dose averages (`ALW2/ALW1`, `ALL2/ALL1`, `ALQD/DOSE`) are unaffected.
- **Never mix routes across a ratio.** Numerator and denominator must share a particle set
  — `ALW1/ALWF`, never `ALW1/ALFL` (water numerator over a medium denominator).

For proton-only and per-species keys, see the [scorer key reference](#scorer-key-reference).

### Example input

All bins use the `ALL-PART` generalized particle; the routines do the real particle
selection and weighting internally. Here is the LET-to-water moment pair plus its fluence
bin, giving both dose- and track-averaged LET to water:

```text
* activate FLUSCW (WHAT(3)=1)
USERWEIG          0.0       0.0       1.0                              &
*
* 1st water-LET moment (weight LET)    -> ALW1 (unit 22)
USRBIN           11.0  ALL-PART      -22.  <xmax ymax zmax  bins...>    ALW1
USRBIN         <xmin ymin zmin> ...                                    &
* 2nd water-LET moment (weight LET^2)  -> ALW2 (unit 23)
USRBIN           11.0  ALL-PART      -23.  <xmax ymax zmax  bins...>    ALW2
USRBIN         <xmin ymin zmin> ...                                    &
* unweighted fluence, same set         -> ALWF (unit 21)
USRBIN           11.0  ALL-PART      -21.  <xmax ymax zmax  bins...>    ALWF
USRBIN         <xmin ymin zmin> ...                                    &
```

Swap `ALW*` for `ALL*` to score the local medium instead. For dirty dose or dose-averaged
Qeff, add `WHAT(6)=1` to `USERWEIG` (this calls `COMSCW`) and add a plain `DOSE` (or
`DOSE-H2O`) bin as the denominator.

## Post-processing

Merge the per-cycle binary files and convert them with
[pymchelper](https://datamedsci.github.io/pymchelper/)'s `convertmc`, as the
repository-root [README](../README.md) describes:

```bash
pip install pymchelper
convertmc image --many "*_fort.2*"     # image, txt, plotdata, ...
```

Then form the recipe-table ratios on the converted per-unit arrays — e.g.
`LETd = (unit 23) / (unit 22)`, `LETt = (unit 22) / (unit 21)`, dirty fraction =
`ALDD unit / DOSE unit`. Dose-averaged LET needs no separate fluence bin: a segment's dose
is ∝ length·LET, so `ALW2/ALW1` (or `ALL2/ALL1` in the medium) is already the dose-weighted
mean of LET.

## Scorer key reference

Which particles a key sees follows from its LET route, not its name (see
[Implementation details](#implementation-details)). e± and neutrals are excluded from every
key. Moment suffix `1`/`2` (and `W1`/`W2` for water) = first / second raw LET moment,
weighted by LET / LET²; units keV/µm and (keV/µm)². Material applies only to the
LET-weighted keys, where a stopping power is evaluated; for fluence and Qeff it is N/A.

### Handled by `FLUSCW` (fluence / track-length weighting)

| Key | Weight | Scores | Material |
|---|---|---|---|
| `ALL1` / `ALL2` | LET / LET² | all charged hadrons + ions | local |
| `ALFL` | 1 | all charged hadrons + ions (denominator for `ALL1`/`ALL2`) | *N/A* |
| `ALW1` / `ALW2` | LET / LET² | p, d, t, ³He, ⁴He | water |
| `ALWF` | 1 | p, d, t, ³He, ⁴He (denominator for `ALW1`/`ALW2`) | *N/A* |
| `ALQ1` | Qeff | all charged hadrons + ions | *N/A* |
| `ALQF` | 1 | all charged hadrons + ions (denominator for `ALQ1`) | *N/A* |
| `PAL1` / `PAL2` | LET / LET² | all protons | local |
| `PAW1` / `PAW2` | LET / LET² | all protons | water |
| `P1FL` | 1 / 0 | primary protons (`LTRACK=1`) | *N/A* |
| `P1L1` / `P1L2` | LET / LET² | primary protons | local |
| `P1W1` / `P1W2` | LET / LET² | primary protons | water |
| `D2L1` / `D2L2` | LET / LET² | deuterons | local |
| `T3L1` / `T3L2` | LET / LET² | tritons | local |
| `H3L1` / `H3L2` | LET / LET² | ³He | local |
| `H4L1` / `H4L2` | LET / LET² | ⁴He / α | local |
| `L6L1` / `L6L2` | LET / LET² | Li-6 (Z=3, A=6) | local |
| `L7L1` / `L7L2` | LET / LET² | Li-7 (Z=3, A=7) | local |
| `L6FL` / `L7FL` | 1 / 0 | Li-6 / Li-7 fluence | *N/A* |

### Handled by `COMSCW` (dose weighting)

| Key | Weight | Scores |
|---|---|---|
| `ALDD` | 1 / 0 | dirty dose, LET > 30 MeV cm²/g (LET local on `DOSE`, water on `DOSE-H2O`) |
| `ALQD` | Qeff | dose weighted by Qeff |
| `ALDQ` | 1 / 0 | dirty dose, Qeff > 30 |
| `P1DO` | 1 / 0 | primary-proton dose (`JTRACK=1`, `LTRACK=1`) |
| `L6DO` / `L7DO` | 1 / 0 | Li-6 / Li-7 dose |

`P1DO` keeps only source-generation protons — unlike `AUXSCORE PROTON`, which also counts
secondary protons. `ALDD` follows its own binning: `DOSE` judges LET in the local medium
(all fragments), `DOSE-H2O` in water (p/d/t/³He/⁴He only); any other binning is rejected
with a warning.

## Compilation

Compile the routine and link it into a custom FLUKA executable with the FLUKA build tools.
Compiled and linked with **FLUKA 4** (`fff` + `lfluka`):

```bash
export FLUPRO=/usr/local/fluka        # your FLUKA install
export PATH=$PATH:$FLUPRO/bin
fff fluka_let_scoring.f               # -> fluka_let_scoring.o
lfluka -m fluka -o flukalet fluka_let_scoring.o
```

The include files are referenced by their literal names (`INCLUDE 'dblprc.inc'`), which is
what `fff` expects. The exact tool names vary between FLUKA distributions (older `ldpmqmd`
wrappers work too).

## Testing

`../tests/` contains a ready-to-run example, `plan01_field01_geoA_SOBPcent.inp`, that
exercises the scorer keys (proton, light-fragment, and lithium LET moments, the Qeff
moments, and the dose/fluence filters). It uses a `SOURCE` reading the accompanying
`sobp.dat` spot list, so the executable must also link the source sampler from
`../fluka_source`:

```bash
fff ../fluka_source/source_sampler.f
lfluka -m fluka -o flukalet fluka_let_scoring.o source_sampler.o
cd ../tests && rfluka -e ../flukalet -N0 -M1 plan01_field01_geoA_SOBPcent
```

On the first scoring call the run prints `fluka_let_scoring: water MWATER = <n>` to the
`.out`, confirming the auto-detected water material.

## Adding or changing a scorer

1. Keep the key to four characters and add it to the right dispatch block (`FLUSCW` for
   fluence/track-length, `COMSCW` for dose).
2. Document: particle selection; all vs. primary-only; local vs. water; LET, LET², or
   filter weight.
3. Update the input cards and any post-processing that reads the output units.
4. Compile-test, then run a small smoke test before production use.

## Implementation details

Everything below is *why*, not *how* — skip it unless a number surprises you.

### The LET reported: unrestricted

Every LET key reports **unrestricted LET (LET_∞)** — full electronic stopping power,
including the energy handed to δ-rays. State this when you report a number: the Kalholm
review warns the restricted/unrestricted choice routinely goes unstated.

Because the LET already contains the δ-ray energy, **e± and neutrals are excluded from
every key**, even under EMF transport. FLUKA still transports the δ-rays — that is what
puts the dose in the right place — but scoring them *again* as LET carriers in their own
right would count that same energy twice.

### Two LET routes: `GETLET` and `TRACKR`

No single route covers the whole field, so the code dispatches over two:

| Route | Used for | Why |
|---|---|---|
| `GETLET` | p, d, t, ³He, ⁴He | Unrestricted by construction (the restriction-energy argument is passed as zero). Verified: 5.2 MeV cm²/g for a 160 MeV proton in water, matching NIST PSTAR's unrestricted value. |
| `TRACKR` (`ΣDTRACK/ΣTTRACK`) | Li and heavier fragments | `GETLET` returns exactly zero for these: FLUKA transports every heavy ion under one generic code (`JTRACK = -39`) with the real Z/A off in `FHEAVY`, and `GETLET`'s argument list has nowhere to accept them. |

`TRACKR` measures energy *deposited*, which in general is *restricted* at the δ-ray
production threshold (100 keV under `PRECISION`). For heavy fragments that distinction is
void: a δ exceeds 100 keV only when β²γ² > 0.098 (roughly 45 MeV/u), and fragments in a
proton field are far slower (typically < 1 MeV/u), so no δ is ever split off and the
deposited LET *is* the unrestricted LET. Both routes return the same quantity, and the seam
between them is invisible in the output.

**Limit of that argument.** It rests on *speed*, and the `TRACKR` branch takes everything
that is not p/d/t/³He/⁴He — heavy fragments, and also π±/K±/µ± where a field is energetic
enough to make them. Anything routed there that is fast enough to produce a δ above the
threshold (β²γ² > 0.098) is scored as restricted LET, not unrestricted. Sound for proton
therapy, where fragments are < 1 MeV/u and no pions are produced below ~290 MeV; re-derive
before trusting it for a 400 MeV/u carbon beam (T_max ≈ 800 keV) or any high-energy field.

The split is not academic. Before it, when `ALL1`/`ALL2` used `TRACKR` for everything,
entrance track-averaged LET came out **7% low** — measured 0.932 of the unrestricted water
reference, against 0.931 predicted by Bethe for Δ=100 keV vs T_max=378 keV at 160 MeV. With
the two routes it is 0.985, the remainder being the genuine solid-water-vs-water material
difference.

### Material: local vs water

- **Local keys** use `MATLET = MEDFLK(NREG,1)` — LET in whatever material the particle is
  currently in. The only guard skips vacuum / non-material regions (`MATLET ≤ 0` or
  `RHO ≤ 0`). Restrict *where* you score with the `USRBIN` geometry, not a material list.
- **Water keys** use `MATLET = MWATER`, resolved **automatically** on the first scoring
  call — no hardcoded material numbers. `MWATER` is taken from FLUKA's built-in `MATQLT`
  (the "extra water material for Q(L) calculations" in `flkmat.inc`, present even when the
  input defines no explicit `WATER`); failing that, from the first material named `WATER`
  in `MATNAM`. The chosen index is written to the `.out` as
  `fluka_let_scoring: water MWATER = <n>` so you can verify it. If no water material is
  found, the water keys return zero and a warning is printed.

### Dirty-dose thresholds

Dirty dose is the dose deposited by particles above a radiation-quality threshold. There
are two dirty-dose scorers, on two different quantities, each a single named `PARAMETER` in
`COMSCW` — change one line and recompile:

| Scorer | Threshold quantity | Parameter | Value |
|---|---|---|---|
| `ALDD` | unrestricted LET (mass stopping power) | `DDTHRE` | 30 MeV cm²/g = 3 keV/µm in water |
| `ALDQ` | Qeff = z_eff²/β² (dimensionless) | `DQTHRE` | 30 |

Both are "30" **by construction, not by coincidence**, and they are not the same 30 —
different quantities in different units. `DQTHRE` was deliberately anchored to `DDTHRE`: a
proton at the LET threshold (30 MeV cm²/g, T ≈ 16.9 MeV, β ≈ 0.188) has Qeff ≈ 28.4,
rounded up to 30. So for protons the two definitions cut at roughly the same place; for
heavier ions they diverge, because Qeff's `z_eff²` and LET's stopping power are different
functions of energy. On a 160 MeV proton SOBP the two dirty-dose fractions track closely by
depth (both ≈ 6% at the entrance, ≈ 97% at the Bragg peak), differing only by Bethe's
slowly-varying log term, which Qeff's pure `1/β²` scaling does not carry. Read Heuchel et
al. (2024) before changing `DDTHRE`.

`ALDD` reads its binning's generalized particle (`IDUSBN`) and follows it — the material the
LET is judged in is not a separate key:

| `USRBIN` WHAT(2) | LET judged in | Route | Covers |
|---|---|---|---|
| `DOSE` (228) | local medium | `GETLET` + `TRACKR` | all charged hadrons + ions, incl. heavy fragments |
| `DOSE-H2O` (252) | water | `GETLET` | p, d, t, ³He, ⁴He only |

The two meaningless cross combinations (dose-to-water thresholded on medium LET, or vice
versa) cannot be expressed; any other binning is rejected with a warning.

Caveats before quoting a dirty-dose number:

- **`DOSE-H2O` misses heavy fragments** — `GETLET` has no water stopping power for Li and
  above, so they drop out. That bites hardest here, because fragments are exactly the
  high-LET component dirty dose is meant to capture. The `DOSE` variant has no such gap.
- **Point-like depositions are excluded** — `TRACKR` carries no step data for them, so
  their LET cannot be reconstructed. Two are genuinely high-LET and thus under-counted:
  heavy recoils (`JTRACK` 208) and low-energy neutron kerma (308).
- **e± excluded**, as everywhere, to avoid double-counting the δ-ray energy already inside
  the primary's LET.

> **`DOSE-H2O` requires a `RAD-BIOL` card**, or FLUKA dies with a SIGFPE in
> `dedx/alphbt.f:382` (`alphbt` called with `rbealp=-1`, `rbebet=0`) via `score/usrsco.f`
> — no message about the real cause, and only in materials other than water, so a water
> phantom runs and a carbon one aborts. This is FLUKA's own behaviour: it reproduces with
> the stock `fluka` executable and no user routines. The α/β file content is *not* used by
> `DOSE-H2O` — a FLUKA developer
> [confirmed on the forum](https://fluka-forum.web.cern.ch/t/dose-to-water-issue/6128) that
> it is required only "for technical reasons — to be overcome at a later stage". Any valid
> file will do; `../tests/letbio.dat` is a minimal one.

### Qeff

Qeff is a radiation-quality metric used the same way as LET (as an RBE proxy, and to
threshold dirty dose) but computed completely differently — a closed-form function of the
particle's own charge and speed alone, with **no material and no stopping-power table**:

```text
beta  = v/c
z_eff = z * ( 1 - exp( -125 * beta * |z|^(-2/3) ) )     (Barkas)
Qeff  = z_eff^2 / beta^2
```

That material-independence buys two things the LET keys lack:

- **One code path for every charged particle**, fragments included — no `GETLET` coverage
  gap and no `TRACKR` restricted/unrestricted question (see `QEFCAL` in the source). β comes
  from `PTRACK/ETRACK` (no rest-mass lookup, which is awkward for fragments); z comes from
  `ICHRGE`, or from `FHEAVY` (`ICHEAV`) for heavier fragments under `JTRACK = -39` — the
  same fork the Li keys use.
- **`ALQD`/`ALDQ` need no binning fork.** Unlike `ALDD`, which must branch on `DOSE` vs
  `DOSE-H2O` because the material changes the LET, Qeff is identical on any dose-like
  binning.

The LET dose-average trick (`ALL2/ALL1`, no dose bin needed) does **not** carry over to
Qeff: it worked only because a segment's dose *is* length × LET, so a second LET moment
reproduces the dose weighting for free. Qeff carries no such relationship to the energy
deposited, so dose-averaged Qeff has to weight the real dose directly, in `COMSCW`
(`ALQD`). Kalholm et al. (2023, 2024) show Qeff and related energy/quality metrics improve
RBE prediction over LET alone for proton therapy; this implementation scores only the
moments and is not tied to a particular RBE model.

### A note on ancestry

`LTRACK` gives the generation (1 = source-generation), and the isotope filters classify the
*currently transported* particle. None of these record the production vertex, parent
particle, or reaction channel — that would require production-time tagging via `STUPRF` or
`MDSTCK`.

## References

**Averaged LET (the LET-moment scorers):**
> Kalholm F, Grzanka L, Traneus E, Bassler N. A systematic review on the usage of averaged LET in radiation biology for particle therapy. Radiotherapy and Oncology. 2021;161:211-221.

**Dirty dose (the `ALDD` scorer):**
> Heuchel L, Hahn C, Ödén J, Traneus E, Wulff J, Timmermann B, Bäumer C, Lühr A. The dirty and clean dose concept: towards creating proton therapy treatment plans with a photon-like dose response. Medical Physics. 2024;51(1):622-636.
  Introduces the concept, and discusses the choice of threshold.

**Dirty dose applied to RBE:**
> Kalholm F, Toma-Dasu I, Traneus E. "Dirty dose"-based proton variable RBE models — performance assessment on in vitro data. Medical Physics. 2025;52(2):1311-1322.

**Qeff (the `ALQ1`/`ALQF`/`ALQD`/`ALDQ` scorers):**
> Kalholm F, Grzanka L, Toma-Dasu I, Bassler N. Modeling RBE with other quantities than LET significantly improves prediction of in vitro cell survival for proton therapy. Medical Physics. 2023;50(1):651-659.

**Extended Qeff metrics:**
> Kalholm F, et al. Novel radiation quality metrics accounting for proton energy spectra for RBE proton models. Medical Physics. 2024;51(8):5773-5782.
