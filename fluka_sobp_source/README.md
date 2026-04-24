## Particle Source for Pencil Beam Scanning in Hadrontherapy

A FLUKA user-defined source implementation for simulating pencil beam scanning in hadrontherapy applications. This source creates spread-out Bragg peak (SOBP) shapes in depth by simulating multiple beamlets with configurable properties.

### Features

- **Beamlet Configuration**: Multiple beamlets with individual relative weights, energies, and sizes
- **Energy Spread**: Optional energy spread specification per beamlet
- **Beam Geometry**: Choice between divergent beam from nozzle plane or point-like virtual source
- **Configurable Input**: Requires external data file (typically `sobp.dat`) for beam geometry and kinematics


### Compilation
Compile using the FLUKA utility:

```bash
ldpmqmd -oflukadpm_sobp source_sampler.f
```

Place the `sobp.dat` file in the same directory as your FLUKA input file, then activate the custom source with a SOURCE card in your input file and run:

```bash
rfluka -N0 -M1 -e flukadpm_sobp your_input_file
```

**Note**: This implementation is based on the template from `$FLUPRO/usermvax/source.f`.

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
- SDUM : filename for the spotlist; `sobp.dat` if not set.

### Spotlist format:
Spotlists can be generated from DICOM `RTPLAN` files together with a suitable beam model.
We here use the tool [dicomexport](https://github.com/nbassler/dicomexport) to generate the spotlists.

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
