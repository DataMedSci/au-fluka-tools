[![CI](https://github.com/DataMedSci/au-fluka-tools/actions/workflows/ci.yml/badge.svg)](https://github.com/DataMedSci/au-fluka-tools/actions/workflows/ci.yml)

## AUFLUKATOOLS ##

Auxiliary tools for the Monte Carlo particle transport code
[FLUKA](https://www.fluka.org).

### FLUKA user routines ###

- [`fluka_let_scoring/`](fluka_let_scoring/) — `FLUSCW`/`COMSCW` routines for LET-moment
  scoring (track- and dose-averaged LET) and dirty dose. See its
  [README](fluka_let_scoring/README.md) for the scorer keys and post-processing.
- [`fluka_sobp_source/`](fluka_sobp_source/) — `SOURCE` routine sampling a spread-out
  Bragg peak from a spot list.

### Running FLUKA on a cluster ###

- `rcfluka.py` — parallelisation of FLUKA runs using the CONDOR job queuing system.
- `rtfluka.sh` — parallelisation on TORQUE.

### Post-processing USRBIN and friends ###

Use [pymchelper](https://datamedsci.github.io/pymchelper/). Its `convertmc` command parses,
merges and plots USRBIN, USRTRACK and USRBDX output:

```bash
pip install pymchelper
convertmc image --many "*.bnn"
```

The scripts this repository used to carry for that job (`usrbin2ascii.py`,
`usrtrack2ascii.py`, `usrbinmerge.py`, `flukaplot.py`, `dicom2vxl.py`, and the `usrbinav` /
`usrtrackav` / `usrbdxav` wrappers, together with a vendored copy of `flair`) have been
removed: they had to track FLUKA's output format by hand and had fallen behind, and
`convertmc` covers the same ground. See issue #8. They remain in the git history if needed.
