# amfsplit

Long-read SSU–ITS–LSU metabarcoding of arbuscular mycorrhizal fungi: from accurate long reads to
names in the three systems surveys use (MaarjAM VT, UNITE SH, LSU species clade), a culture-linked
assignment, and how far that assignment can be trusted.

It grows out of O'Brien, Marín & Tedersoo, *Detecting unknown arbuscular mycorrhizal fungi from
linked long-read rRNA gene sequences* (preprint). The paper's code and data are in
[amf-novelty-benchmark](https://github.com/ayobi/amf-novelty-benchmark), release 0.3.1,
doi:10.5281/zenodo.23422616.

## Status

v0.1 in progress. v0.1 ports the paper's frozen method, and it ships only when it reproduces
release 0.3.1 assignment for assignment.

| Step | State |
|---|---|
| `regions`: Stage 1 (amfsplit.py v0.3.1) in Rust | Done. Byte-identical on the 12 parity cases, the 904 deposited culture sequences, the 280 LongGloDB ASVs and the 5,546 EUKARYOME reads (`parity/stage1`) |
| `validate` + `triage`: AMF or not, from the best hit of the LR0R–FLR2 window (validate_windows.py and triage.py v0.3.1) | Done. Byte-identical on 7 parity cases, vsearch included. New `--amf-reference`, and no default AMF rule, because the prototype's default outgroups do not fit V18 (`parity/triage`) |
| `gate0`: frozen naming method and per-copy parity | Gate 0a (nearest-reference LSU arm) and Gate 0b steps 1 and 2 (LSU species clades; species clades plus the ITS distance rule N6, the paper's method) are byte-identical in both culture-only holdouts (`amfsplit gate0`), and so are all 52 LSU-removal and genus rounds of amf-novelty-benchmark 0.3.2, 208 per-sequence tables (`--benchmark` or `AMFSPLIT_BENCHMARK` pointing at that checkout), and so is the LongGloDB test, conditions A and B under both VT analyses, 8 tables (`--lgdb` or `AMFSPLIT_LGDB`: a local folder with the test windows and evidence, which are never redistributed), and so is the EUKARYOME test, 8 tables over 5,511 reads (`--euk` or `AMFSPLIT_EUK`: a local folder with its four test windows). Gate 0b is complete: every per-call table behind the paper is reproduced byte for byte |
| `place`, `place-summary`: LSU windows to the placement table the naming arms read | `place` runs MAFFT `--addfragments --keeplength` and EPA-ng on the hash-checked frozen V18 reference with the frozen arguments and thread counts (MAFFT 4, EPA-ng 8), then the frozen `lsu_placement.py` summarize; on the 903 culture windows it compares the result with the frozen table (`1de103b6…`) cell by cell. A different `best_edge` is accepted only where that edge is tied for the top weight in the run: EPA-ng lists tied placements in run-dependent order, and no naming step reads `best_edge`. The summary alone is byte-identical on the 51 published rounds of the LSU-removal test (`amfsplit gate0 --policy placement --m4-runs`) |
| `evidence` (library; `gate0 --policy evidence`): VT, SH and LSU evidence from BLAST+ against MaarjAM, UNITE (Glomeromycota) and V18 AMF-only (frozen `arms.py` and `arm_sets.py`) | Rebuilds the frozen evidence of all 903 culture copies in all three arms from their windows and the three hash-pinned databases (BLAST+ 2.17.0 as in the frozen runs). The full UNITE release is streamed, never held in memory. `--evidence-arm` and `--evidence-sample` check one arm or an even sample of copies, since each query's hits depend only on the query and the database |
| `name`: production user-facing VT, SH and LSU naming workflow | Next; depends on Gate 0 and reference-search input wiring |
| `db`: references built from the user's downloads, checked by checksum | To do |

## Build and run

Needs a Rust toolchain (edition 2024) and a C++11 compiler (edlib is compiled from source).

```bash
cargo build --release
target/release/amfsplit regions -i reads.fq.gz -a reads.anchors.tsv -p panels/amf_legacy_windows.tsv \
    -o out/sample --only-anchored --open-ends
# AMF triage on the LR0R-FLR2 windows, against Delavaux V18 with its outgroups
target/release/amfsplit validate -w out/sample.LSU_LROR-FLR2.fasta -r V18_LSUDB_052025.fasta \
    -o out/sample.LSU_LROR-FLR2.vs_v18 --threads 8
target/release/amfsplit triage --hits out/sample.LSU_LROR-FLR2.vs_v18.hits.tsv --out out/sample \
    --records out/sample.records.tsv --fasta-prefix out/sample \
    --amf-reference V18_LSUDB_052025_AMFONLY.fasta --min-id 90      # 85 for ONT R9.4.1
cargo test --release
# Run the frozen naming reproducibility gate (requires parity/naming fixtures)
cargo run --release -- gate0
```

`-a` takes ITSxRust's `--anchors-tsv` for the same reads (ITSxRust 0.3.0, as in the paper). Options
and outputs are those of amfsplit.py v0.3.1:

- `{out}.records.tsv`
- one FASTA per ITS region and per primer window
- `{out}.summary.json`
- `{out}.run.json`, which records the tool version, options and SHA-256 of inputs and outputs

## Layout

| Path | What it holds |
|---|---|
| `crates/amfsplit` | The command-line tool |
| `crates/amfsplit-regions` | Stage 1: anchors, primer windows, concatemer screen |
| `crates/amfsplit-naming` | Frozen VT/SH/LSU culture-bridge and novelty logic plus the Gate 0 parity harness |
| `crates/amfsplit-edlib` | FFI to edlib, vendored from Python edlib 1.3.9.post1 so alignments equal the prototype's |
| `panels/` | The legacy primer-window panel used in the paper |
| `parity/stage1/`, `parity/triage/`, `parity/naming/` | Golden outputs of the prototypes, frozen evidence and per-copy naming parity checks |

## Licence

MIT (`LICENSE`). edlib, vendored in `crates/amfsplit-edlib/vendor/edlib`, is © Martin Šošić, MIT.
