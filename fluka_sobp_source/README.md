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
$FLUPRO/flutil/ldpm3qmd source_sampler.f -o flukadpm3_sobp
```

Place the `sobp.dat` file in the same directory as your FLUKA input file, then activate the custom source with a SOURCE card in your input file and run:

```bash
rfluka -N0 -M1 -e flukadpm3_sobp your_input_file
```

**Note**: This implementation is based on the template from `$FLUPRO/usermvax/source.f`.
