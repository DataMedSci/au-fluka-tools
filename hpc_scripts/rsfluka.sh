#!/bin/bash
#
# SLURM job-array submission for FLUKA 4.
#
# Submit one FLUKA cycle per array task:
#
#   sbatch --array=1-20 rsfluka.sh
#
# Override the defaults through the environment:
#
#   sbatch --array=1-20 \
#          --export=ALL,INPUT=plan01,EXE=$PWD/flukalet,AUX="sobp.dat" \
#          rsfluka.sh
#
# Each task runs in its own directory (run_001, run_002, ...) with its own
# random seed, so tasks never share or overwrite a file. Results are merged
# afterwards -- see README.md. Score to BINARY output (negative USRBIN unit)
# or the outputs cannot be merged.
#
#SBATCH --job-name=fluka
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

set -euo pipefail

INPUT="${INPUT:-example}"   # FLUKA input file, WITHOUT the .inp suffix
EXE="${EXE:-}"              # user executable; empty = stock fluka
AUX="${AUX:-}"              # extra runtime files, e.g. "sobp.dat"

: "${FLUPRO:?set FLUPRO to your FLUKA installation, e.g. /usr/local/fluka}"
i="${SLURM_ARRAY_TASK_ID:?submit with: sbatch --array=1-N rsfluka.sh}"

# Resolve before descending into the run directory, so relative values work.
[ -n "$EXE" ] && EXE="$(readlink -f -- "$EXE")"

run="$(printf 'run_%03d' "$i")"
mkdir -p "$run"
cp -- "${INPUT}.inp" "$run/"
if [ -n "$AUX" ]; then
    # shellcheck disable=SC2086  # AUX is deliberately word-split
    cp -- $AUX "$run/"
fi
cd "$run"

# Give every task its own random sequence, otherwise all N jobs repeat the
# same history and merging them only multiplies one result instead of
# improving the statistics. FLUKA takes the seed from RANDOMIZE WHAT(2);
# deriving it from the array index keeps the whole run reproducible.
if ! grep -q '^RANDOMIZ' "${INPUT}.inp"; then
    echo "error: no RANDOMIZE card in ${INPUT}.inp -- every task would use" >&2
    echo "       the same seed and produce identical results." >&2
    exit 1
fi
# Fixed-format card: name in cols 1-10, WHAT(1) in 11-20, WHAT(2) in 21-30.
seed_card="$(printf 'RANDOMIZ  %10.1f%10.1f' 1.0 "$i")"
sed -i "s|^RANDOMIZ.*|${seed_card}|" "${INPUT}.inp"

if [ -n "$EXE" ]; then
    "$FLUPRO/bin/rfluka" -N0 -M1 -e "$EXE" "$INPUT"
else
    "$FLUPRO/bin/rfluka" -N0 -M1 "$INPUT"
fi
