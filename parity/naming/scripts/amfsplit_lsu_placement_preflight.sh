#!/usr/bin/env bash
# Read-only environment check for the frozen Gate 0 LSU-placement rebuild.
set -u

printf '%s\n' '== AMFSplit LSU placement preflight =='
ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
printf 'Project: %s\n' "$ROOT"

# Explicit executable variables take precedence over PATH lookup.
MAFFT_BIN="${MAFFT_BIN:-$(command -v mafft 2>/dev/null || true)}"
RAxML_NG_BIN="${RAXML_NG_BIN:-$(command -v raxml-ng 2>/dev/null || true)}"
EPA_NG_BIN="${EPA_NG_BIN:-$(command -v epa-ng 2>/dev/null || true)}"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python 2>/dev/null || command -v python3 2>/dev/null || true)}"
EXPECTED_RAXML_SHA256='b5038949eb005d14c48efe47f1ddb4fa750e89dbafeb30c46e05c651bacba6d5'

show_tool() {
  local label="$1" path="$2" version_arg="$3"
  printf '\n--- %s ---\n' "$label"
  if [[ -n "$path" && -x "$path" ]]; then
    printf 'Path: %s\n' "$path"
    "$path" "$version_arg" 2>&1 | head -n 5 || true
    sha256sum "$path" || true
  else
    printf 'NOT FOUND or not executable: %s\n' "${path:-<unset>}"
  fi
}

show_tool 'MAFFT' "$MAFFT_BIN" '--version'
show_tool 'RAxML-NG' "$RAxML_NG_BIN" '--version'
show_tool 'EPA-ng' "$EPA_NG_BIN" '--version'

printf '\n--- Pinned RAxML-NG checksum ---\n'
if [[ -n "$RAxML_NG_BIN" && -f "$RAxML_NG_BIN" ]]; then
  actual="$(sha256sum "$RAxML_NG_BIN" | awk '{print $1}')"
  if [[ "$actual" == "$EXPECTED_RAXML_SHA256" ]]; then
    printf 'PASS: %s\n' "$actual"
  else
    printf 'FAIL: expected %s\nActual:   %s\n' "$EXPECTED_RAXML_SHA256" "$actual"
  fi
else
  printf 'Cannot check: RAXML_NG_BIN is not set to an existing file.\n'
fi

printf '\n--- Python / Biopython / Edlib ---\n'
if [[ -n "$PYTHON_BIN" ]]; then
  "$PYTHON_BIN" -c 'import sys, Bio, edlib; print("Executable:", sys.executable); print("Python:", sys.version.split()[0]); print("Biopython:", Bio.__version__); print("Edlib:", edlib.__version__ if hasattr(edlib, "__version__") else "import OK")' 2>&1 || true
else
  printf 'Python interpreter not found.\n'
fi

printf '\n--- V18 files ---\n'
if [[ -n "${V18_DIR:-}" ]]; then
  printf 'V18_DIR=%s\n' "$V18_DIR"
  for name in V18_LSUDB_052025_AMFONLY.fasta V18_LSUDB_052025.fasta Root_V18_LSUDB_052025.newick; do
    [[ -f "$V18_DIR/$name" ]] && printf 'FOUND: %s\n' "$name" || printf 'MISSING: %s\n' "$name"
  done
else
  printf 'V18_DIR is unset. Candidate AMF-only FASTA files under HOME:\n'
  find "$HOME" -type f -name 'V18_LSUDB_052025_AMFONLY.fasta' -print 2>/dev/null | sort
fi

printf '\n--- Committed query / golden hashes ---\n'
for rel in \
  parity/naming/inputs/benchmark/ref/fasta/LSU_LROR-FLR2.fasta \
  parity/naming/placement/lsu_placements.copies.tsv; do
  if [[ -f "$ROOT/$rel" ]]; then sha256sum "$ROOT/$rel"; else printf 'MISSING: %s\n' "$ROOT/$rel"; fi
done

printf '\nNote: MAFFT/EPA-ng hashes are printed for provenance recording. Their executable hashes did not match the archived provenance lock on the original machine, although the 20 Oct 2026 rebuild outputs matched the lock exactly.\n'
