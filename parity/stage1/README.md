# Stage 1 parity

`amfsplit regions` is a port of the Stage 1 prototype, `amfsplit.py` v0.3.1 from
amf-novelty-benchmark (`stage1/amfsplit_v0.3.1/`). Its outputs must be byte-identical to the
prototype's: `records.tsv`, the 13–15 region and window FASTAs, and `summary.json`. The only
addition is `{out}.run.json`, which records versions, options and input checksums.

## What is here

| Path | What it holds |
|---|---|
| `python/amfsplit_v0_3_1.py` | The prototype, unchanged (SHA-256 `a7d79a63…`, the file in the benchmark repo at commit `02c0844`) |
| `inputs/stefani903/` | The 903 retained sequences of culture reference v4 (`reference/culture_ref_v4/fasta/full.fasta`, SHA-256 `db4e8c2f…`) and the ITSxRust anchors of the 904 deposited sequences (`benchmarks/stefani/out/itsx/stefani.anchors.tsv`, `49b57124…`), gzipped |
| `inputs/synthetic/` | Synthetic reads, anchors and panels written by `make_golden.py` (seed 7; reads A–N are those of the prototype's `tests/test_v03.py`) |
| `golden/<case>/` | The prototype's `records.tsv`, `summary.json` and `MANIFEST.tsv` (SHA-256 of every output); synthetic cases keep their FASTAs too |
| `cases.tsv` | Reads, anchors, panel and options of each case |
| `make_golden.py` | Rebuilds the inputs and golden files (Python 3 with edlib 1.3.9) |
| `check_parity.py` | Runs the Rust binary on every case, with a field-level diff on failure |
| `run_external.sh` | Prototype and port side by side on data that cannot live here |

## Cases

Twelve cases: nine synthetic, three on the 903 culture sequences.

The synthetic reads cover clean reads, primers missing or open at either end, a decoy primer in the
wrong region, tandem, inverted, partial and different-molecule concatemers, spurious repeats,
minus-strand reads with and without anchors, a strandless anchors row, lowercase input with N,
an empty confidence, a partial ITS chain, a read starting at the ITS, short and tiny reads, a
duplicate read ID, primers out of order, a too-short window and an empty read.

The cases run these under the prototype's flags: `--keep-concatemers`, `--no-region-limit`,
`--no-concat-screen`, `--keep-primers` without `--open-ends`, `--only-anchored`, region margins of
0 and 120, a FASTQ.gz input, every concatemer parameter, a panel with region columns and comment
lines, and no anchors at all.

The real cases use the options of the deposited copies (HiFi: `--max-err-frac 0.15 --only-anchored
--open-ends`), the ONT edit budget (0.20), and `--keep-primers --keep-concatemers`.

## Running

```bash
cargo test --release                        # includes the parity test (no Python needed)
python3 -I parity/stage1/check_parity.py    # the same, with a field diff on failure
```

## Not covered here, and how to cover it

- **PX215117.1.** It was dropped from culture reference v4 and GenBank is not reachable from this
  workspace, so `stefani903` has 903 of the 904 deposited sequences. Check 1 of
  `run_external.sh` runs all 904.
- **LongGloDB v1 and EUKARYOME v2.1.** Their sequences are not redistributed. Checks 2 and 3 of
  `run_external.sh` compare the port with the prototype and with the v0.3.1 records stored in
  amf-novelty-benchmark.
- **Empty reads.** With the concatemer screen on, the prototype fails with a TypeError on a
  zero-length read. The port stops with an explicit error at the same point. Case
  `syn_empty_read` covers an empty read with the screen off.
- **Lone CR line endings** (pre-OS X Mac files). Python's universal newlines treat a lone CR as a
  line break and the port does not. CRLF files behave the same in both.
