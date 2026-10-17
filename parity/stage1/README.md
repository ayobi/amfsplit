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

## External checks (`run_external.sh`)

Run on 16 October 2026 on Aaron's workstation: Ubuntu 24.04, g++ 13.3.0, rustc 1.97.1, Python edlib
1.3.9.post1. The anchors are those of the original runs. The expected records are the v0.3.1
files in amf-novelty-benchmark.

| Check | Reads | Processed | Concatemers | Result |
|---|---|---|---|---|
| `stefani904`: all deposited culture sequences, PX215117.1 from GenBank | 904 | 904 | 0 | 15 outputs identical to the prototype |
| `lgdb`: LongGloDB v1 test ASVs | 280 | 280 | 0 | 15 outputs identical; `records.tsv` equals `benchmarks/longglodb/stage1/windows/lgdb_test.records.tsv` |
| `euk`: EUKARYOME v2.1 Glomeromycota queries | 5,546 | 5,541 | 1 | 15 outputs identical; `records.tsv` equals `benchmarks/eukaryome/stage1/euk_test.records.tsv` |

The benchmark repo stores no v0.3.1 records for the 904 deposit (its `stefani.records.tsv` is a
v0.2 output), so that check compares the port with the prototype only.

## Not covered here, and how to cover it

- **PX215117.1, LongGloDB v1 and EUKARYOME v2.1.** They are not in this repository (one is absent
  from the reference, the others are not redistributed). They are covered by the external checks
  above.
- **Empty reads.** With the concatemer screen on, the prototype fails with a TypeError on a
  zero-length read. The port stops with an explicit error at the same point. Case
  `syn_empty_read` covers an empty read with the screen off.
- **Lone CR line endings** (pre-OS X Mac files). Python's universal newlines treat a lone CR as a
  line break and the port does not. CRLF files behave the same in both.
