[![Build Status](https://travis-ci.org/grzanka/au-fluka-tools.svg?branch=master)](https://travis-ci.org/grzanka/au-fluka-tools)

## AUFLUKATOOLS ##
Provides auxiliary tools for the Monte Carlo particle transport code FLUKA (www.fluka.org).
Auflukatools features
- `usrbin2ascii.py` is an alternative script which formats binary usrbin output to a more plotting friendly ASCII file. Provide the usrbin filename as command argument.
- scripts for averaging multiple binary output files without user interaction,
- `rcfluka.py` for automatic parallelization of FLUKA runs on computer clusters using the CONDOR job queuing system.
- `rtfluka.sh` for parallelization on TORQUE

When running `usrbinav`, `usrtrackav` or `usrbdxav`, you must mention the unit numbers you wish to average over, i.e.

```bash
$ usrbinav 60 61 62
```
which will average all results from binary output unit 60, 61 and 62.

