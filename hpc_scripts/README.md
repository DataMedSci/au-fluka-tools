## HPC submission scripts ##

Helpers for spreading a FLUKA run across many cluster nodes. FLUKA itself is serial: you
parallelise it by running N independent copies with different random seeds and merging the
results afterwards. That is all these scripts do — set up the copies and submit them.

| Script | Queuing system |
|---|---|
| `rsfluka.sh` | SLURM — the usual choice on academic HPC clusters |
| `rcfluka.py` | HTCondor — high-throughput computing; CERN's batch system |

All paths target **FLUKA 4** (`$FLUPRO/bin`, `$FLUPRO/data`).

### rsfluka.sh (SLURM) ###

One array task per FLUKA cycle:

```bash
sbatch --array=1-20 rsfluka.sh
```

Defaults come from the environment, so nothing needs editing:

```bash
sbatch --array=1-20 \
       --export=ALL,INPUT=plan01_field01_geoA_SOBPcent,EXE=$PWD/flukalet,AUX="sobp.dat" \
       rsfluka.sh
```

- `INPUT` input file **without** the `.inp` suffix (default `example`)
- `EXE` user executable; leave unset to use the stock `fluka`
- `AUX` extra files the routines need at run time, e.g. the spot list

Each task runs in its own directory (`run_001`, `run_002`, …) and rewrites the
`RANDOMIZE` card with a seed taken from the array index. Two reasons that matters:

- **Separate directories** — tasks never share or overwrite each other's files.
- **Distinct seeds** — without them every task repeats the same history, and merging
  N identical results multiplies one answer rather than improving the statistics. The
  script refuses to run if the input has no `RANDOMIZE` card, rather than quietly
  producing N copies of the same thing.

Because the seed is the array index rather than a random number, a rerun reproduces the
same set of histories.

### rcfluka.py (HTCondor) ###

Compiles and links the user routines, writes a Condor submit file, and queues one job per
node, each with its own seed. Run it from the directory holding your input file, and give
the input **without** the `.inp` extension:

```bash
rcfluka.py -M20 plan01_field01_geoA_SOBPcent \
    -s ../fluka_let_scoring/fluka_let_scoring.f,../fluka_source/source_sampler.f \
    -f sobp.dat
```

- `-M` number of jobs to submit; `-N` start numbering at N
- `-s` user routines to compile and link in (comma-separated)
- `-f` extra files the routines need at run time — the spot list `sobp.dat` here
- `-l` linker, e.g. `ldpmqmd`, when the default is not right
- `-t` dry run: prepare everything but do not submit. Worth using first
- `-m` mail address to notify on completion

`rcfluka.py -h` lists the rest.

### Output must be BINARY, or you cannot merge ###

Each job writes its own output, and the whole point of running N of them is to combine
those into one result with better statistics. **That merge only works on binary
(unformatted) output.**

In FLUKA, `USRBIN` WHAT(3) is the logical output unit, and its *sign* picks the format:

- **negative** — unformatted binary. Use this.
- positive — formatted ASCII. Cannot be merged; you get N text files and no way to
  combine them properly.

```text
*         WHAT(1)   WHAT(2)   WHAT(3)   WHAT(4)   WHAT(5)   WHAT(6)   SDUM
USRBIN       10.0     228.0     -21.0      1.00      1.00     10.25   DOSE_ZN
                                 ^^^^^
                                 negative = binary = mergeable
```

The same applies to `USRTRACK`, `USRBDX` and the other estimators.

### Merging afterwards ###

Not supplied here. Use [pymchelper](https://datamedsci.github.io/pymchelper/), whose
`convertmc` command merges the per-job binary files and converts or plots the result:

```bash
pip install pymchelper
convertmc image --many "*_fort.21"
```

### Status ###

None of these are exercised by CI — there is no freely redistributable FLUKA to test
against, and no queuing system on a public runner. They have been updated to the FLUKA 4
layout (`$FLUPRO/bin`, `$FLUPRO/data`; FLUKA 3's `flutil/` no longer exists), and
`rsfluka.sh` has been tested end-to-end against a stub `rfluka`: run directories, seed
rewriting, argument handling and both guards. What has *not* been tested is `sbatch`
itself, so treat the SLURM script as unproven on a live cluster.

`rcfluka.py` dates from 2010 and is the less exercised of the two: its FLUKA 4 paths have
been corrected, but it has not been run against a live Condor pool. Its `-t` flag prepares
a submission without queuing it, which is the cheap way to check it.

A TORQUE/PBS script was dropped in favour of these two: TORQUE development wound down long
ago and sites have largely moved to SLURM.
