# AMFSplit Gate 0: independent LSU-placement rebuild

This guide documents how to independently rebuild and verify the frozen LSU placement artifacts used by AMFSplit Gate 0. It does not modify naming rules or the Rust implementation.

## Current reproducibility result

On 20 October 2026, the independent rebuild passed against the archived lock. Verification checked the pinned V18 source commit, frozen placement script, input and intermediate artifact hashes, fitted reference tree and model, edge annotations, and all columns for all 903 query rows. The verifier reported:

```text
PLACEMENT REBUILD: PASS (all pinned artifacts and all 903 per-query summary rows match exactly)
```

The successful output directory on the original machine was `$HOME/Downloads/amfsplit_lsu_pinned_raxml_20261020_retry1`. Keep that directory and its build log as local evidence, outside the source tree.

## Locked inputs and reference

- V18 repository commit: `f4a0014336f49be6aea48e8374dae052360a88bd`.
- The V18 data directory must contain `V18_LSUDB_052025_AMFONLY.fasta`, `V18_LSUDB_052025.fasta`, and `Root_V18_LSUDB_052025.newick`.
- Query FASTA: `parity/naming/inputs/benchmark/ref/fasta/LSU_LROR-FLR2.fasta` (903 records; SHA-256 checked by the verifier).
- Frozen builder: `parity/naming/frozen_sources/m2b/lsu_placement.py`.
- Expected tool versions: MAFFT 7.526, RAxML-NG 1.2.2, EPA-ng 0.3.8; Python 3.12.14 and Biopython 1.88 were recorded for the reference environment. Edlib is also required by the frozen helper code.
- Required RAxML-NG executable SHA-256: `b5038949eb005d14c48efe47f1ddb4fa750e89dbafeb30c46e05c651bacba6d5`.

The successful run used the verified RAxML-NG executable from the official v1.2.2 release archive. Do not use an arbitrary `raxml-ng` from `PATH` solely because its version string says 1.2.2: the Conda executable previously tested had a different checksum and produced different fitted outputs.

**Tool-provenance note:** the installed MAFFT and EPA-ng reported the expected version strings and produced outputs that passed every locked artifact and per-query comparison. Their local executable hashes did not match the hashes recorded in `software_provenance.json`, so byte-identical executable provenance for those tools remains unresolved. Record actual versions and hashes when reproducing; do not claim their binary checksums match the lock without independently verifying that.

## 1. Set up paths and preflight

Run from the AMFSplit repository root with the intended Python environment activated. The original environment was named `amf-m4`; another machine may use a different name.

```bash
cd /path/to/amfsplit
conda activate amf-m4  # if this environment exists on this machine

# Set this to the directory containing the three V18 files above.
export V18_DIR="/path/to/2024_AMFPipeline_ASV"

# Override these variables as needed for your installed tools.
export MAFFT_BIN="$(command -v mafft)"
export EPA_NG_BIN="$(command -v epa-ng)"
export RAXML_NG_BIN="/path/to/pinned/raxml-ng"

bash parity/naming/scripts/amfsplit_lsu_placement_preflight.sh
```

Confirm the pinned RAxML-NG executable exists and has the required SHA-256 before running the builder:

```bash
printf '%s  %s\n' \
  'b5038949eb005d14c48efe47f1ddb4fa750e89dbafeb30c46e05c651bacba6d5' \
  "$RAXML_NG_BIN" | sha256sum -c -
```

`V18_DIR` must point to the directory containing the three files, not merely the repository top level. If the external source must be cloned, use a new directory so an existing checkout is not modified, then pin it to the specified commit:

```bash
V18_ROOT="$HOME/references/AMF-LSU-Database-and-Pipeline2-f4a001"
git clone https://github.com/c383d893/AMF-LSU-Database-and-Pipeline2.git "$V18_ROOT"
git -C "$V18_ROOT" checkout --detach f4a0014336f49be6aea48e8374dae052360a88bd
git -C "$V18_ROOT" rev-parse HEAD
export V18_DIR="$V18_ROOT/2024_AMFPipeline/2024_AMFPipeline_ASV"
```

The reported commit must equal `f4a0014336f49be6aea48e8374dae052360a88bd`.

## 2. Rebuild into a new output directory

Do not reuse a populated output directory. Keep generated files outside the Git checkout.

```bash
cd /path/to/amfsplit
: "${V18_DIR:?Set V18_DIR to the V18 data directory}"
: "${MAFFT_BIN:?Set MAFFT_BIN to the selected MAFFT executable}"
: "${RAXML_NG_BIN:?Set RAXML_NG_BIN to the checksum-verified RAxML-NG executable}"
: "${EPA_NG_BIN:?Set EPA_NG_BIN to the selected EPA-ng executable}"

QUERY="$PWD/parity/naming/inputs/benchmark/ref/fasta/LSU_LROR-FLR2.fasta"
OUT="$HOME/amfsplit_lsu_rebuild_$(date +%Y%m%d_%H%M%S)"

printf '%s\n' '=== Tool versions and checksums ==='
"$MAFFT_BIN" --version 2>&1 | head -n 3
"$RAXML_NG_BIN" --version 2>&1 | head -n 5
"$EPA_NG_BIN" --version 2>&1 | head -n 5
sha256sum "$MAFFT_BIN" "$RAXML_NG_BIN" "$EPA_NG_BIN"

PYTHONPATH="$PWD/parity/naming/frozen_sources${PYTHONPATH:+:$PYTHONPATH}" \
python parity/naming/frozen_sources/m2b/lsu_placement.py build \
  --v18 "$V18_DIR" \
  --queries "$QUERY" \
  --out "$OUT" \
  --mafft "$MAFFT_BIN" \
  --raxml-ng "$RAXML_NG_BIN" \
  --epa-ng "$EPA_NG_BIN" \
  --threads 8 \
  --mass 0.95

printf 'Rebuild output: %s\n' "$OUT"
```

Record actual tool paths, versions, checksums, and build log with local run notes. If the terminal may close, redirect output to a log or use `nohup`; do not start a second build while one is running.

## 3. Verify the rebuild

Use the verifier stored in the repository rather than a separately downloaded copy:

```bash
python parity/naming/scripts/verify_lsu_placement_rebuild.py \
  --project "$PWD" \
  --run "$OUT"
```

A full pass checks the pinned reference commit and frozen placement script, required file hashes, fitted reference tree and model, query alignment, edge annotations, the 903-row summary, 495 reference tips, 987 tree edges, 0.95 credible mass, 903/903 alignment-eligible queries, and exact per-cell agreement across every column.

If verification fails, it exits nonzero and can write `placement_row_differences.tsv` into the rebuild directory. Preserve outputs and logs while investigating. Do not change thresholds or automatically accept a mismatch as a placement tie.

## Interpretation boundary

This is a reproducibility test for the archived placement computation. It does not independently validate the taxonomy, prove that every species assignment is biologically correct, or constitute external biological validation. Likelihood-weight ratios (LWRs) are conditional on the chosen reference, alignment, and fitted model. Do not tune novelty thresholds during this gate.
