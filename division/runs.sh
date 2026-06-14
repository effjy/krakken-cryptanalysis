#!/usr/bin/env bash
# Krakken-2048 bit-based division-property integral search.
# One round count per invocation. Designed to be launched and left running.
#
#   ./runs.sh <rounds> [cube] [max_bits] [solver_time_limit]
#   ./runs.sh validate                 # run the correctness gates
#
# cube forms:
#   word:W        saturate all 64 bits of word W            (default word:0)
#   byte:W:B      saturate byte B (8 bits) of word W
#   bits:i,j,k    saturate an explicit list of global bit indices
#   dim:D[:seed]  saturate D pseudo-random bits (reproducible via seed)
#
# max_bits           test only the first N output bits (0 = all 2048)
# solver_time_limit  per-bit SCIP cap in seconds (0 = none). A bit that hits the
#                    cap is logged 'unknown' and treated conservatively as NOT
#                    balanced -- so partial / overnight results are always sound.
#
# Examples:
#   ./runs.sh 1                  # round 1, full word:0 cube, all bits
#   ./runs.sh 3 dim:30:0 0 300   # round 3, 30-bit cube, all bits, 300s/bit cap
#   ./runs.sh 2 word:0 64        # round 2, first 64 output bits (quick smoke)
#
# Output streams to the terminal AND to runs/round_<N>.log.
# A 12 GB virtual-memory cap guards against a runaway solve (override MEM_KB=...).
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
PY="${PY:-$HOME/venv/bin/python}"
MEM_KB="${MEM_KB:-12582912}"

# Build libkrakken.so + facets.pkl on first use (idempotent).
make -s PY="$PY" all

if [ "${1:-}" = "validate" ]; then
  ( ulimit -v "$MEM_KB"; make -s PY="$PY" validate )
  ( ulimit -v "$MEM_KB"; "$PY" cbdp.py )     # S-box selector gate (slow-ish)
  exit 0
fi

ROUNDS="${1:?usage: ./runs.sh <rounds> [cube] [max_bits] [solver_time_limit]}"
CUBE="${2:-word:0}"
MAXBITS="${3:-0}"
TLIM="${4:-0}"

mkdir -p runs
LOG="runs/round_${ROUNDS}.log"
echo "[runs] rounds=$ROUNDS cube=$CUBE max_bits=$MAXBITS solver_time_limit=$TLIM"
echo "[runs] S-box model: FACET (compact, exact)  | logging to $LOG  (vmem cap ${MEM_KB} KB)"

ARGS=(--rounds "$ROUNDS" --cube "$CUBE" --max-bits "$MAXBITS" \
      --solver-time-limit "$TLIM" --use-facets)

( ulimit -v "$MEM_KB"; "$PY" krakken_divprop.py "${ARGS[@]}" ) 2>&1 | tee "$LOG"
echo "[runs] done -> $LOG"
