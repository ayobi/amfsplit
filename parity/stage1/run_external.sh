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
# Every run ends with exactly one line starting "ok", "DIFF" or "FAIL"; it never fails silently.
# KEEP=1 keeps both sets of outputs (py.* and rs.*) and the logs for inspection.
#
# The three checks (B = your amf-novelty-benchmark clone):
#   1. the 904 deposited culture sequences:
#        run_external.sh stefani904 stefani904.fasta $B/benchmarks/stefani/out/itsx/stefani.anchors.tsv
#   2. LongGloDB v1 test ASVs (reads from stage3/longglodb_overlap.py):
#        run_external.sh lgdb overlap.test_full.fasta \
#          $B/benchmarks/longglodb/stage1/itsx/lgdb_test.anchors.tsv \
#          $B/benchmarks/longglodb/stage1/windows/lgdb_test.records.tsv
#   3. EUKARYOME v2.1 Glomeromycota queries (euk21_glo_query.fasta) with their ITSxRust anchors:
#        run_external.sh euk euk21_glo_query.fasta euk_test.anchors.tsv \
#          $B/benchmarks/eukaryome/stage1/euk_test.records.tsv
set -uo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
die() { echo "FAIL ${NAME:-?}: $*"; exit 1; }

if [[ $# -lt 3 ]]; then
  sed -n '2,10p' "$0" | sed 's/^# \{0,1\}//'
  exit 2
fi
NAME=$1; READS=$2; ANCHORS=$3; shift 3
EXPECTED=""
if [[ $# -gt 0 && "$1" != "--" ]]; then EXPECTED=$1; shift; fi
[[ $# -gt 0 && "$1" == "--" ]] && shift
OPTS=("$@")
[[ ${#OPTS[@]} -eq 0 ]] && OPTS=(--max-err-frac 0.15 --only-anchored --open-ends)
BIN=${AMFSPLIT_BIN:-$ROOT/target/release/amfsplit}
PANEL=${PANEL:-$ROOT/panels/amf_legacy_windows.tsv}

# inputs first: a wrong path is the commonest failure
[[ -f "$READS" ]] || die "reads not found: $READS"
[[ -f "$ANCHORS" ]] || die "anchors not found: $ANCHORS"
[[ -z "$EXPECTED" || -f "$EXPECTED" ]] || die "expected records not found: $EXPECTED"
[[ -f "$PANEL" ]] || die "panel not found: $PANEL"
[[ -x "$BIN" ]] || die "Rust binary not found: $BIN (run: cargo build --release)"
python3 -c "import edlib" 2> /dev/null || die "python3 has no edlib (activate the env that ran amfsplit.py)"

W=$(mktemp -d "${TMPDIR:-/tmp}/amfsplit-external-$NAME.XXXXXX")
[[ -n "${KEEP:-}" ]] || trap 'rm -rf "$W"' EXIT
A="$ANCHORS"
if [[ "$A" == *.gz ]]; then gzip -dc "$A" > "$W/anchors.tsv" || die "cannot decompress $ANCHORS"; A="$W/anchors.tsv"; fi

if ! python3 -I "$HERE/python/amfsplit_v0_3_1.py" -i "$READS" -a "$A" -p "$PANEL" -o "$W/py" "${OPTS[@]}" 2> "$W/py.log"; then
  echo "     prototype log (last lines):"; tail -n 5 "$W/py.log" | sed 's/^/     /'
  die "the prototype failed"
fi
if ! "$BIN" regions -i "$READS" -a "$ANCHORS" -p "$PANEL" -o "$W/rs" "${OPTS[@]}" 2> "$W/rs.log"; then
  echo "     amfsplit log (last lines):"; tail -n 5 "$W/rs.log" | sed 's/^/     /'
  die "the Rust port failed"
fi

counts=$(python3 -c "import json,sys; d=json.load(open(sys.argv[1])); print(f\"{d['reads']} reads, {d['reads_processed']} processed, {d['concatemers']['flagged']} concatemers\")" "$W/rs.summary.json")
fail=0; n=0
for f in "$W"/py.*; do
  base=${f#"$W/py."}
  [[ "$base" == "log" ]] && continue
  n=$((n + 1))
  if [[ ! -f "$W/rs.$base" ]]; then echo "     missing from the port: $base"; fail=1
  elif ! cmp -s "$f" "$W/rs.$base"; then echo "     differs: $base"; fail=1; fi
done
if [[ -n "$EXPECTED" ]]; then
  if cmp -s "$W/rs.records.tsv" "$EXPECTED"; then echo "     records.tsv equals $EXPECTED"
  else echo "     records.tsv differs from $EXPECTED"; fail=1; fi
fi
if [[ $fail -eq 0 ]]; then
  echo "ok   $NAME: $n output files identical to the prototype ($counts)"
else
  [[ -n "${KEEP:-}" ]] && echo "     outputs kept in $W" || echo "     rerun with KEEP=1 to keep the outputs"
  echo "DIFF $NAME ($counts)"
fi
exit $fail
