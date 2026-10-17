#!/usr/bin/env bash
# Stage 1 parity on data that cannot live in this repository (run where the data are).
#
#   parity/stage1/run_external.sh NAME READS ANCHORS [EXPECTED_RECORDS] [-- OPTIONS...]
#
# Runs the frozen prototype (python/amfsplit_v0_3_1.py, needs python3 with edlib 1.3.9) and the Rust
# port (target/release/amfsplit) on the same reads and ITSxRust anchors, and requires every output
# file to be byte-identical. With EXPECTED_RECORDS, the Rust records.tsv must also equal that file
# (the v0.3.1 records stored in amf-novelty-benchmark). OPTIONS default to the options of the
# deposited copies: --max-err-frac 0.15 --only-anchored --open-ends.
#
# The three checks to run (paths as in amf-novelty-benchmark; B = that repository):
#   1. the 904 deposited culture sequences (adds PX215117.1, absent from this workspace):
#        run_external.sh stefani904 stefani_deposit.fasta $B/benchmarks/stefani/out/itsx/stefani.anchors.tsv
#   2. LongGloDB v1 test ASVs (reads from stage3/longglodb_overlap.py, never committed):
#        run_external.sh lgdb bench/longglodb/overlap.test_full.fasta \
#          $B/benchmarks/longglodb/stage1/itsx/lgdb_test.anchors.tsv \
#          $B/benchmarks/longglodb/stage1/windows/lgdb_test.records.tsv
#   3. EUKARYOME v2.1 Glomeromycota queries (euk21_glo_query.fasta, SHA-256 aaab407b...; anchors from
#      the local ITSxRust 0.3.0 run, hifi preset, tolerant thresholds as in run_stage1.sh):
#        run_external.sh euk euk21_glo_query.fasta bench/eukaryome/stage1/itsx/euk_test.anchors.tsv \
#          $B/benchmarks/eukaryome/stage1/euk_test.records.tsv
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
NAME=${1:?name}; READS=${2:?reads}; ANCHORS=${3:?anchors}; shift 3
EXPECTED=""
if [[ $# -gt 0 && "$1" != "--" ]]; then EXPECTED=$1; shift; fi
[[ $# -gt 0 && "$1" == "--" ]] && shift
OPTS=("$@")
[[ ${#OPTS[@]} -eq 0 ]] && OPTS=(--max-err-frac 0.15 --only-anchored --open-ends)
BIN=${AMFSPLIT_BIN:-$ROOT/target/release/amfsplit}
PANEL=${PANEL:-$ROOT/panels/amf_legacy_windows.tsv}
[[ -x "$BIN" ]] || { echo "build first: cargo build --release" >&2; exit 2; }

W=$(mktemp -d "${TMPDIR:-/tmp}/amfsplit-external-$NAME.XXXXXX")
[[ -n "${KEEP:-}" ]] || trap 'rm -rf "$W"' EXIT   # KEEP=1 keeps both sets of outputs for diffing
A="$ANCHORS"
if [[ "$A" == *.gz ]]; then gzip -dc "$A" > "$W/anchors.tsv"; A="$W/anchors.tsv"; fi   # the prototype reads plain text

python3 -I "$HERE/python/amfsplit_v0_3_1.py" -i "$READS" -a "$A" -p "$PANEL" -o "$W/py" "${OPTS[@]}" 2> "$W/py.log"
"$BIN" regions -i "$READS" -a "$ANCHORS" -p "$PANEL" -o "$W/rs" "${OPTS[@]}" 2> "$W/rs.log"

fail=0; n=0
for f in "$W"/py.*; do
  base=${f#"$W/py."}
  [[ "$base" == "log" ]] && continue
  n=$((n + 1))
  if ! cmp -s "$f" "$W/rs.$base"; then echo "DIFF $NAME: $base"; fail=1; fi
done
if [[ -n "$EXPECTED" ]]; then
  if cmp -s "$W/rs.records.tsv" "$EXPECTED"; then echo "ok   $NAME: records.tsv equals $EXPECTED"
  else echo "DIFF $NAME: records.tsv differs from $EXPECTED"; fail=1; fi
fi
if [[ $fail -eq 0 ]]; then
  echo "ok   $NAME: $n output files identical to the prototype ($(grep -c . "$W/rs.records.tsv") records lines)"
else
  echo "rerun with KEEP=1 to keep the prototype (py.*) and port (rs.*) outputs in $W"
fi
exit $fail
