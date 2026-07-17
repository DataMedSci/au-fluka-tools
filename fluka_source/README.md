## Spot-list Particle Source for Pencil Beam Scanning

A FLUKA user-defined `SOURCE` routine that samples primaries from a **spot list**: an
external table of beamlets, each with its own energy, position, size, divergence and
weight.

Nothing about the routine is specific to a spread-out Bragg peak — it samples whatever
spot list you give it, from a single pencil beam to a full scanned field. A SOBP is simply
the most common case, and is what the default file name (`sobp.dat`) and the example in
`../tests/` reflect: a set of beamlets whose energies and weights combine into a flat
depth-dose plateau.

### Features

- **Beamlet Configuration**: Multiple beamlets with individual relative weights, energies, and sizes
- **Energy Spread**: Optional energy spread specification per beamlet
- **Beam Geometry**: Choice between divergent beam from nozzle plane or point-like virtual source
- **Configurable Input**: Requires external data file (typically `sobp.dat`) for beam geometry and kinematics


### Compilation
Compile and link with the FLUKA 4 tools:

```bash
export FLUPRO=/usr/local/fluka      # your FLUKA installation
export PATH=$PATH:$FLUPRO/bin

fff source_sampler.f                # -> source_sampler.o
lfluka -m fluka -o fluka_source source_sampler.o
```

Use `ldpmqmd` in place of `lfluka` if you need the DPMJET/RQMD event generators
(heavy-ion projectiles); for protons `lfluka` is enough.

To use this source together with the LET scorers, link both objects into one
executable — see [`../fluka_let_scoring/README.md`](../fluka_let_scoring/README.md).

Place the `sobp.dat` file in the same directory as your FLUKA input file, then activate the custom source with a SOURCE card in your input file and run:

```bash
rfluka -N0 -M1 -e fluka_source your_input_file
```

**Note**: This implementation is based on the template from `$FLUPRO/src/user/source.f`.

### Invoking
In the FLUKA input file, you can specify the `SOURCE` card with a few arguments
```
*...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
SOURCE
```
- WHAT(1) : Nonzero flag enabling the virtual-source geometry based on `SADx` and `SADy`; `0` leaves that mode off. If not explicitly set, this routine will see it as `0`.
- WHAT(2) : Flag for debug info (`0` = off, nonzero = on)
- WHAT(3) : `SADx` - distance from the X-scanning magnet to the spotlist plane; must be positive when `WHAT(1)` is nonzero
- WHAT(4) : `SADy` - distance from the Y-scanning magnet to the spotlist plane; must be positive when `WHAT(1)` is nonzero
- WHAT(5) : (Not used)
- WHAT(6) : (Not used)
- SDUM : filename for the spotlist; `sobp.dat` if not set. **Maximum 8 characters** —
  FLUKA passes the SOURCE SDUM through `CHARACTER*8 SDUSOU`, so a longer name is
  truncated (`sobp.dat` is exactly 8). The file is read from the directory you launch
  `rfluka` in, not from the temporary working directory it creates.

### Spotlist format:
Spotlists can be generated from DICOM `RTPLAN` files with
[dicomexport](https://github.com/nbassler/dicomexport).

A plan alone is not enough: `RTPLAN` gives spot positions and monitor units, but not the
energy, spot size or divergence the machine actually delivers for them. dicomexport folds
in a **beam model** to supply those, and ships beam models for DCPT — so for that facility
the spot list can be produced straight from the plan, with no extra calibration work.

The FLUKA source sampler so far expects a file named `sobp.dat` with either 5,6,7,9 or 11 columns describing every spot:

```
5: ENERGY, XPOS, YPOS, FWHMxy, WEIGHT
6: ENERGY, XPOS, YPOS, FWHMX, FWHMY, WEIGHT
7: ENERGY, DE, XPOS, YPOS, FWHMX, FWHMY, WEIGHT
9: ENERGY, DE, XPOS, YPOS, FWHMX, FWHMY, DIVX, DIVY, WEIGHT
11: ENERGY, DE, XPOS, YPOS, FWHMX, FWHMY, DIVX, DIVY, CORX, CORY, WEIGHT
```

Where,
- `ENERGY`  particle energy in GeV/nucleon
- `DE`  particle energy spread (sigma) in GeV/nucleon
- `XPOS` beam spot center (X coordinate), in cm
- `YPOS` beam spot center (Y coordinate), in cm
- `FWHMX` beam spot size in X axis, in cm
- `FWHMY` beam spot size in Y axis, in cm
- `DIVX` beam spot angular divergence in X axis, in mrad
- `DIVY` beam spot angular divergence in Y axis, in mrad
- `CORX` correlation coefficient rho(x,tx) (dimensionless)
- `CORY` correlation coefficient rho(y,ty) (dimensionless)
- `PART` beamlet weight (relative, but recommend to use absolute primary particle numbers here for future features)

Headers with `#` will be skipped.
