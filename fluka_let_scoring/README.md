# FLUKA LET and Qeff scoring routines

FLUKA `FLUSCW` and `COMSCW` user-weighting routines for scoring LET and Qeff moments, from
which averaged radiation-quality quantities (track-averaged and dose-averaged) are
reconstructed in post-processing.

These are **user routines**, not a standalone program. You link them into a custom FLUKA
executable, then activate them from your input file with a `USERWEIG` card.

References:

- Averaged LET (the LET-moment scorers): Kalholm F, Grzanka L, Traneus E, Bassler N.
  *A systematic review on the usage of averaged LET in radiation biology for particle
  therapy.* Radiotherapy and Oncology. 2021;161:211-21.
- Dirty dose (the `ALDD` scorer): Heuchel L, Hahn C, Ödén J, Traneus E, Wulff J,
  Timmermann B, Bäumer C, Lühr A. *The dirty and clean dose concept: towards creating
  proton therapy treatment plans with a photon-like dose response.* Medical Physics.
  2024;51(1):622-36. Introduces the concept, and discusses the choice of threshold.
- Dirty dose applied to RBE: Kalholm F, Toma-Dasu I, Traneus E. *'Dirty dose'-based proton
  variable RBE models — performance assessment on in vitro data.* Medical Physics.
  2025;52(2):1311-22.
- Qeff (the `ALQ1`/`ALQF`/`ALQD`/`ALDQ` scorers): Kalholm F, Grzanka L, Toma-Dasu I,
  Bassler N. *Modeling RBE with other quantities than LET significantly improves
  prediction of in vitro cell survival for proton therapy.* Medical Physics.
  2023;50(1):651-9.
- Kalholm F, et al. *Novel radiation quality metrics accounting for proton energy spectra
  for RBE proton models.* Medical Physics. 2024;51(8):5773-82.

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

The headline scorers are `ALL1` and `ALL2`: the first and second **unrestricted** LET
moments over **all charged hadrons and ions**, evaluated in the *local* material. They
cover every charged hadron and ion FLUKA transports — protons, light ions, and heavy
fragments alike — by dispatching internally over two routes (`GETLET` for p/d/t/³He/⁴He,
`TRACKR` for heavier fragments; see *The LET reported* below for why, and why both yield
the same unrestricted quantity).

Neutral particles are skipped, and **electrons and positrons are excluded** even when EMF
transport is active. This is not a matter of taste: LET here is *unrestricted* stopping
power, which already includes the energy handed to δ-rays. FLUKA transporting those δ-rays
is fine and desirable — it is how the dose ends up in the right place — but scoring them
*again* as LET carriers in their own right would count that same energy twice. Excluding
e± is what keeps the accounting closed. See Kalholm et al. on why the restricted /
unrestricted distinction has to be stated explicitly.

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

## The LET reported: unrestricted, in one of two materials

**Every scorer here reports unrestricted LET (LET_∞)** — the full electronic stopping
power, including the energy handed to δ-rays. The only axis that varies between keys is
the *material* the LET is evaluated in: the **local** medium, or **water**.

That matters because the restricted/unrestricted choice is exactly what the Kalholm review
warns goes unstated. State it when you report: these are unrestricted.

Getting there needs two implementation routes, because neither covers the whole field:

| Route | Used for | Why |
|---|---|---|
| `GETLET` | p, d, t, ³He, ⁴He | Unrestricted by construction (the restriction-energy argument is passed as zero). Verified: 5.2 MeV cm²/g for a 160 MeV proton in water, matching NIST PSTAR's unrestricted value. |
| `TRACKR` (`ΣDTRACK/ΣTTRACK`) | Li and heavier fragments | `GETLET` **cannot** serve these: FLUKA transports every heavy ion under one generic code (`JTRACK = -39`) with the real Z/A off in `FHEAVY`, and `GETLET`'s argument list has nowhere to accept them — it returns exactly zero for all of them. |

The `TRACKR` route measures energy *deposited*, which in general is *restricted* at the
δ-ray production threshold (100 keV under `PRECISION`). **For heavy fragments that
distinction is void**: a δ can only exceed 100 keV when β²γ² > 0.098, i.e. above roughly
45 MeV/u, and fragments in a proton field are far slower (typically < 1 MeV/u). No δ is
ever split off, so the deposited LET *is* the unrestricted LET. Both routes return the
same quantity, and the seam between them is invisible in the output.

> **Limit of that argument.** It rests on *speed*, and the `TRACKR` branch takes everything
> that is not p/d/t/³He/⁴He — heavy fragments, and also π±/K±/µ± where a field is energetic
> enough to make them. Anything routed there that is fast enough to produce a δ above the
> threshold (β²γ² > 0.098) is scored as restricted LET, not unrestricted. Sound for proton
> therapy, where fragments are < 1 MeV/u and no pions are produced below ~290 MeV; re-derive
> before trusting it for a 400 MeV/u carbon beam (T_max ≈ 800 keV) or any high-energy field.

The gap is not academic. Before the split, when `ALL1`/`ALL2` used `TRACKR` for
everything, entrance track-averaged LET came out **7% low** — measured 0.932 of the
unrestricted water reference, against 0.931 predicted by Bethe for Δ=100 keV vs
T_max=378 keV at 160 MeV. With the two routes it is 0.985, the remainder being the genuine
solid-water-vs-water material difference.

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

## Particle coverage at a glance

Which particles a key actually scores follows from its LET route, not from its name. The
`TRACKR` route reconstructs LET from the energy deposited per unit step, so it works for
anything charged; the `GETLET` route needs a tabulated stopping power and so is limited to
the five light species FLUKA supplies.

| Keys | Scores | Does **not** score | Route |
|---|---|---|---|
| `ALL1`, `ALL2`, `ALFL`, `ALDD` on `DOSE` | every charged particle FLUKA transports except e±: p, d, t, ³He, ⁴He, Li **and all heavier fragments**, plus π±/K±/µ± if present | e±, neutrals, point-like depositions (208, 211, 308) | `GETLET` local for p/d/t/³He/⁴He, `TRACKR` for the rest |
| `ALW1`, `ALW2`, `ALWF`, `ALDD` on `DOSE-H2O` | p, d, t, ³He, ⁴He | **Li and heavier fragments**, π±/K±/µ±, e±, neutrals | `GETLET`, water |
| `PAL1`, `PAL2` | protons, all generations | everything else | `GETLET`, local |
| `PAW1`, `PAW2` | protons, all generations | everything else | `GETLET`, water |
| `P1FL`, `P1L1`, `P1L2`, `P1W1`, `P1W2`, `P1DO` | protons with `LTRACK=1` (source generation) | secondary protons, everything else | `GETLET` |
| `D2L1`/`D2L2` | deuterons | everything else | `GETLET`, local |
| `T3L1`/`T3L2` | tritons | everything else | `GETLET`, local |
| `H3L1`/`H3L2` | ³He | everything else | `GETLET`, local |
| `H4L1`/`H4L2` | ⁴He / α | everything else | `GETLET`, local |
| `L6L1`, `L6L2`, `L6FL`, `L6DO` | Li-6 (Z=3, A=6) | everything else | `TRACKR` |
| `L7L1`, `L7L2`, `L7FL`, `L7DO` | Li-7 (Z=3, A=7) | everything else | `TRACKR` |
| `ALQ1`, `ALQF`, `ALQD`, `ALDQ` | every charged particle FLUKA transports except e±: p, d, t, ³He, ⁴He, Li **and all heavier fragments** | e±, neutrals, point-like depositions | `QEFCAL` (no material, no stopping-power table) |

Three consequences worth internalising:

- **The water-reference keys cannot see Li or heavier.** `GETLET` has no water stopping
  power for them, so every fragment above ⁴He silently drops out of `ALW1`/`ALW2`/`ALWF`
  and of `ALDD` on a `DOSE-H2O` binning. Only the `TRACKR` keys cover the full field. This
  matters most for dirty dose, where fragments are the high-LET component of interest.
- **e± are excluded everywhere by design**, even with EMF transport active. The LET is
  unrestricted, so the δ-ray energy is already inside the primary's LET; scoring the
  transported δ-rays again would double-count it. Transporting them is still correct and
  wanted — it is what puts the dose in the right place.
- **Never mix routes across a ratio.** Numerator and denominator must share a particle
  set: `ALL1/ALFL`, not `ALL1/ALWF` or `ALL1` over a plain `ALL-PART` bin.

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
| `ALQ1` | all-particle Qeff moment, track-length weighted | charged hadrons + ions (no e±) | *(none)* | Qeff |
| `ALQF` | all-particle fluence (`ALQ1` denominator) | charged hadrons + ions (no e±) | *(none)* | 1 |

See [Effective-charge radiation quality: Qeff](#effective-charge-radiation-quality-qeff)
below for `ALQ1`/`ALQF`, and for the `ALQD`/`ALDQ` pair in `COMSCW`.

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
| `ALQD` | dose weighted by Qeff | charged hadrons + ions (no e±) | Qeff |
| `ALDQ` | dirty dose (Qeff > 30) | charged hadrons + ions (no e±) | 1 or 0 |
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
naturally with the LET scorers above. Heuchel et al. (2024) introduce the dirty/clean dose
concept and its use in planning for a photon-like dose response; Kalholm et al. (2025)
assess dirty-dose-based variable-RBE models against in vitro data.

### The two hardcoded dirty-dose thresholds

There are two dirty-dose scorers, on two different quantities, each with its own hardcoded
threshold. Both are single named `PARAMETER`s in `COMSCW` — change one line and recompile.

| Scorer | Threshold quantity | Parameter | Value |
|---|---|---|---|
| `ALDD` | LET (unrestricted mass stopping power) | `DDTHRE` | **30 MeV cm²/g** = **3 keV/µm in water** |
| `ALDQ` | Qeff = z_eff²/β² (dimensionless) | `DQTHRE` | **30** |

The two values are both "30" **by construction, not by coincidence** — and they are not the
same 30. They are different quantities in different units. `DQTHRE` was deliberately
anchored to `DDTHRE`: a proton at the LET threshold (3 keV/µm in water, i.e. 30 MeV cm²/g,
T ≈ 16.9 MeV, β ≈ 0.188) has **Qeff ≈ 30** (28.4 measured, rounded up). So for protons the
two dirty-dose definitions cut at roughly the same place; for heavier ions they do not,
because Qeff's `z_eff²` and LET's stopping power are different functions of energy. See the
Qeff section below for the full derivation.

Heuchel et al. (2024) discuss the choice of LET threshold — read it before changing
`DDTHRE`.

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
| `DOSE` (228) | local medium | `GETLET` + `TRACKR` | all charged hadrons + ions, incl. heavy fragments |
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
- **Electrons and positrons are excluded**, as for `ALL1`/`ALL2`, and for the same reason:
  the threshold is applied to unrestricted LET, which already contains the δ-ray energy, so
  scoring the δ-rays as dirty-dose carriers in their own right would count it twice.

> **`DOSE-H2O` requires a `RAD-BIOL` card**, or FLUKA dies with a SIGFPE in
> `dedx/alphbt.f:382` (`alphbt` called with `rbealp=-1`, `rbebet=0`) via `score/usrsco.f`
> — no message about the real cause, and only in materials other than water, so a water
> phantom runs and a carbon one aborts. This is FLUKA's own behaviour: it reproduces with
> the stock `fluka` executable and no user routines at all. The α/β file content is *not*
> used by `DOSE-H2O` — a FLUKA developer
> [confirmed on the forum](https://fluka-forum.web.cern.ch/t/dose-to-water-issue/6128) that
> it is required only "for technical reasons — to be overcome at a later stage". Any valid
> file will do; `tests/letbio.dat` is a minimal one.

## Effective-charge radiation quality: Qeff

Qeff is a radiation-quality metric, not a LET: it is used the same way (as an RBE proxy,
and to threshold dirty dose), but computed completely differently.

**Qeff is material-independent and needs no stopping power.** It is a closed-form function
of the particle's own charge and speed alone — nothing about the medium enters it:

```text
beta   = v/c
z_eff  = z * ( 1 - exp( -125 * beta * |z|^(-2/3) ) )     (Barkas)
Qeff   = z_eff^2 / beta^2
```

That is the fundamental difference from every LET key in this file: LET is a stopping
power, so it depends on the material and (via `GETLET`) on a per-species table that only
covers five light ions. Qeff depends on neither.

Kalholm et al. (2023, 2024) show this and related energy/quality metrics improve RBE
prediction over LET alone for proton therapy; this implementation is not tied to a
particular RBE model, it just scores the moments needed to reconstruct Qeff and pairs it
with dirty dose.

Concretely, that material-independence buys two things none of the LET scorers have:

- **One code path for every charged particle**, fragments included, with no `GETLET`
  coverage gap and no `TRACKR` restricted/unrestricted question — see `QEFCAL` in the
  source. β comes from `PTRACK/ETRACK` (needs no rest-mass lookup, which is awkward for
  fragments); z comes from `ICHRGE` for ordinary particles and light ions, or from
  `FHEAVY` (`ICHEAV`) for heavier fragments transported under FLUKA's generic heavy-ion
  code `JTRACK = -39` — the same fork the Li scorers already use.
- **`ALDQ` needs no `IDUSBN` fork.** Compare with `ALDD`, which has to branch on whether
  the binning is `DOSE` or `DOSE-H2O` because the *material* the LET is judged in changes
  the answer. Qeff does not care what the binning scores: `ALQD` and `ALDQ` work
  identically on `DOSE`, `DOSE-H2O`, or any other dose-like generalized particle, with no
  warning and no restriction.

### Track-averaged and dose-averaged Qeff

Both flavours exist, as for LET, but unlike LET the trick used there (`ALL2/ALL1` for the
dose average, no dose bin needed) does **not** carry over: it worked because a segment's
dose *is* `length × LET`, so a second LET moment reproduces the dose weighting for free.
Qeff carries no such relationship to the energy actually deposited, so the dose-averaged
flavour has to weight the real dose directly, in `COMSCW`.

Score three co-located bins:

```text
* Qeff moment, weighted by track length          -> ALQ1  (unit 22)
USRBIN           11.0  ALL-PART      -22.  <xmax ymax zmax bins...>    ALQ1
USRBIN         <xmin ymin zmin> ...                                   &
* unweighted fluence, same particle set          -> ALQF  (unit 21)
USRBIN           11.0  ALL-PART      -21.  <xmax ymax zmax bins...>    ALQF
USRBIN         <xmin ymin zmin> ...                                   &
* dose weighted by Qeff                          -> ALQD  (unit 23)
USRBIN           10.0      DOSE      -23.  <xmax ymax zmax bins...>    ALQD
USRBIN         <xmin ymin zmin> ...                                   &
```

then, paired with an unfiltered `DOSE` bin over the same region:

- **track-averaged Qeff** `= ALQ1 / ALQF`
- **dose-averaged Qeff** `= ALQD / DOSE`

As always, `ALQF` — not a plain `ALL-PART` bin — is the track-average denominator, for the
same reason as `ALFL`: it applies exactly the same acceptance as `ALQ1`.

### Dirty dose by Qeff

`ALDQ` is the Qeff analogue of `ALDD`: dose from particles above a threshold, here on Qeff
rather than on LET.

```text
* dirty dose, Qeff > 30                          -> ALDQ  (unit 24)
USRBIN           10.0      DOSE      -24.  <xmax ymax zmax bins...>    ALDQ
USRBIN         <xmin ymin zmin> ...                                   &
```

then `dirty fraction = ALDQ / DOSE`, exactly as for `ALDD`.

**On the threshold, 30 (dimensionless):** Qeff and LET are different quantities with
different units, so their thresholds are independent choices — they cannot, in general,
be made to agree for every particle and energy. This one was anchored to a proton: at
30 MeV cm²/g of unrestricted mass stopping power (`ALDD`'s threshold) a proton sits at
T ≈ 16.9 MeV, β ≈ 0.188, where Qeff ≈ 28.4 (checked against this file's own `ALW1`/`ALWF`
via `GETLET`). `DQTHRE` is rounded up from that to 30 — slightly *stricter* for protons
than an exact match would be. The constant is a single named parameter in `COMSCW`
(`DQTHRE`); changing it means editing one line and recompiling.

As a consistency check, `ALDQ` and `ALDD` were run on the same 160 MeV proton SOBP: the
resulting dirty-dose fractions track each other closely by depth (e.g. both ≈ 6% at the
entrance, ≈ 97% at the Bragg peak), with small differences from Bethe's slowly-varying log
term, which Qeff's pure `1/β²` scaling does not carry. For a heavy fragment the two would
diverge much more, since Qeff's `z_eff²` dependence and LET's stopping-power `z²`
dependence are not the same function of energy.

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
exercises the scorer keys (proton, light-fragment, and lithium LET moments, the Qeff
moments, and the dose/fluence filters). It uses a `SOURCE` reading the accompanying
`sobp.dat` spot list,
so the executable must also link the source sampler from `../fluka_source`:

```bash
fff ../fluka_source/source_sampler.f
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
