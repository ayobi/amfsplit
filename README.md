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
| `validate` + `triage`: AMF or not, from the best hit of the LR0R–FLR2 window (validate_windows.py and triage.py v0.3.1) | Done. Byte-identical on 7 parity cases, vsearch included (`parity/triage`); the default outgroups do not fit V18 (see there) |
| `name`: VT, SH and LSU names, consensus, decision steps | Next |
| `db`: references built from the user's downloads, checked by checksum | To do |

## Build and run

Needs a Rust toolchain (edition 2024) and a C++11 compiler (edlib is compiled from source).

```bash
cargo build --release
target/release/amfsplit regions -i reads.fq.gz -a reads.anchors.tsv -p panels/amf_legacy_windows.tsv \
    -o out/sample --only-anchored --open-ends
cargo test --release
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
| `crates/amfsplit-edlib` | FFI to edlib, vendored from Python edlib 1.3.9.post1 so alignments equal the prototype's |
| `panels/` | The legacy primer-window panel used in the paper |
| `parity/stage1/`, `parity/triage/` | Golden outputs of the prototypes, and the parity checks |

## Licence

MIT (`LICENSE`). edlib, vendored in `crates/amfsplit-edlib/vendor/edlib`, is © Martin Šošić, MIT.
