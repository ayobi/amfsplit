# amfsplit

Long-read SSU–ITS–LSU metabarcoding of arbuscular mycorrhizal fungi: from accurate long reads to
names in the three systems surveys use (MaarjAM VT, UNITE SH, LSU species clade), a culture-linked
assignment, and how far that assignment can be trusted.

It grows out of O'Brien, Marín & Tedersoo, *Detecting unknown arbuscular mycorrhizal fungi from
linked long-read rRNA gene sequences* (preprint). The paper's code and data are in
[amf-novelty-benchmark](https://github.com/ayobi/amf-novelty-benchmark), release 0.3.1,
doi:10.5281/zenodo.23422616.

## Status

v0.1 reproduces the paper's method assignment for assignment: from the stored evidence and
placements of release 0.3.1 (Gates 0a and 0b), and from its own BLAST+ evidence and EPA-ng
placements (Gate 3). On the 280 LongGloDB test sequences, `amfsplit name` reproduces the published
calls end to end. Release packaging is in progress.

## Quick start: naming long reads

You need the tools at the versions the paper used: ITSxRust 0.3.0, BLAST+ 2.17.0, MAFFT 7.526 and
EPA-ng 0.3.8 (all on Bioconda; `packaging/bioconda/` holds the draft recipe for amfsplit itself). You also need three downloads:

- MaarjAM, the gDAT snapshot of July 2020 (3,081 sequences);
- the UNITE 10.0 general FASTA release of 19 February 2025 (`UNITE_public_19.02.2025.fasta.gz`);
- Delavaux's V18 LSU database, `2024_AMFPipeline/2024_AMFPipeline_ASV/V18_LSUDB_052025.fasta` at
  commit f4a0014 of [AMF-LSU-Database-and-Pipeline2](https://github.com/c383d893/AMF-LSU-Database-and-Pipeline2).

```bash
# 1. Once: check every file against its pin, build the VT and SH BLAST databases, and rebuild the
#    V18 reference alignment exactly as the frozen run did (MAFFT L-INS-i, checked by hash)
amfsplit db -o amfsplit_db --maarjam maarjam.fasta --unite UNITE_public_19.02.2025.fasta.gz \
    --v18-full V18_LSUDB_052025.fasta
# 2. Name: ITSxRust anchors, regions with the paper's Stage 1 options, then the naming
amfsplit name -i reads.fq.gz --hmm F.hmm --db amfsplit_db -o named/sample      # --preset ont for ONT
```

`--hmm` is ITSx's fungal profile (the paper's `F.hmm` has SHA-256 `a8c38c5c...`). With anchors from an
earlier ITSxRust run, pass `-a reads.anchors.tsv` instead. To run the steps separately, `amfsplit
regions` writes the windows and `amfsplit name -w <prefix>` names them.

`named/sample` then holds:

- `calls.tsv`: one row per sequence, with the status (known species, known species group, novel
  species with the genus known, novel lineage, boundary or divergent class), the placement, the
  novelty and artifact claims, and the reasons behind them;
- `rows.tsv`: every column of the frozen method's row that does not need the sequence's true culture;
- `skipped.tsv`: sequences that lack one of the four windows, with what they lack;
- `run.json`: inputs, database manifest and counts.

BLAST and EPA-ng results are kept in `named/sample/work/`, so a rerun only repeats the naming.

**Scope.** A sequence is named only if it has all four windows: the SSU flank, the full ITS, the
5.8S and the LSU LR0R-FLR2 window. That means long amplicons spanning SSU, ITS and LSU, as in the
paper. v0.1 is validated on PacBio HiFi reads and long-read ASVs; ONT consensus sequences are
accepted but untested. Calls are tied to the pinned references, and `db` refuses other files unless
told to accept them, in which case the results are not the paper's method.

**Data.** The culture reference and the frozen files ship with the tool: in
`share/amfsplit/naming` beside an installed binary, or in `parity/naming` in this checkout.
`AMFSPLIT_DATA` points elsewhere.

| Step | State |
|---|---|
| `regions`: Stage 1 (amfsplit.py v0.3.1) in Rust | Done. Byte-identical on the 12 parity cases, the 904 deposited culture sequences, the 280 LongGloDB ASVs and the 5,546 EUKARYOME reads (`parity/stage1`) |
| `validate` + `triage`: AMF or not, from the best hit of the LR0R–FLR2 window (validate_windows.py and triage.py v0.3.1) | Done. Byte-identical on 7 parity cases, vsearch included. New `--amf-reference`, and no default AMF rule, because the prototype's default outgroups do not fit V18 (`parity/triage`) |
| `gate0`: frozen naming method and per-copy parity | Gate 0a (nearest-reference LSU arm) and Gate 0b steps 1 and 2 (LSU species clades; species clades plus the ITS distance rule N6, the paper's method) are byte-identical in both culture-only holdouts (`amfsplit gate0`), and so are all 52 LSU-removal and genus rounds of amf-novelty-benchmark 0.3.2, 208 per-sequence tables (`--benchmark` or `AMFSPLIT_BENCHMARK` pointing at that checkout), and so is the LongGloDB test, conditions A and B under both VT analyses, 8 tables (`--lgdb` or `AMFSPLIT_LGDB`: a local folder with the test windows and evidence, which are never redistributed), and so is the EUKARYOME test, 8 tables over 5,511 reads (`--euk` or `AMFSPLIT_EUK`: a local folder with its four test windows). Gate 0b is complete: every per-call table behind the paper is reproduced byte for byte |
| `place`, `place-summary`: LSU windows to the placement table the naming arms read | `place` runs MAFFT `--addfragments --keeplength` and EPA-ng on the hash-checked frozen V18 reference with the frozen arguments and thread counts (MAFFT 4, EPA-ng 8), then the frozen `lsu_placement.py` summarize; on the 903 culture windows it compares the result with the frozen table (`1de103b6…`) cell by cell. A different `best_edge` is accepted only where that edge is tied for the top weight in the run: EPA-ng lists tied placements in run-dependent order, and no naming step reads `best_edge`. The summary alone is byte-identical on the 51 published rounds of the LSU-removal test (`amfsplit gate0 --policy placement --m4-runs`) |
| `evidence` (library; `gate0 --policy evidence`): VT, SH and LSU evidence from BLAST+ against MaarjAM, UNITE (Glomeromycota) and V18 AMF-only (frozen `arms.py` and `arm_sets.py`) | Rebuilds the frozen evidence of all 903 culture copies in all three arms from their windows and the three hash-pinned databases (BLAST+ 2.17.0 as in the frozen runs). The full UNITE release is streamed, never held in memory. `--evidence-arm` and `--evidence-sample` check one arm or an even sample of copies, since each query's hits depend only on the query and the database |
| `gate3`: the paper's calls from the tool's own evidence and placements | Rebuilds the VT, SH and LSU evidence of the 903 culture copies with BLAST+ (best units by the frozen first-in-BLAST-order rule) and their placements with MAFFT and EPA-ng, then names both culture-only holdouts with species clades plus N6 and compares them with the goldens. Built evidence and placements are cached in `--work`, so reruns repeat only the naming |
| `arm`: Protocol v3 rules, opt-in | R2, the inside gate, applies to any named-row table and refits its limits per round. It reproduces its development estimates exactly: EUKARYOME false known 13 to 1 of 2,706 (kept 152 to 149 of 199), LongGloDB 4 of 102 (kept 136 to 135 of 164), culture holdout kept 745 to 738 of 755, species holdout 18 of 803, LSU-removal rounds 53 to 48 of 673 (full-reference control 9 of 673). Each distance is computed once and shared by every round, so all of these run in seconds. R1a, unnamed cultures as taxa, reproduces its estimates (EUKARYOME false known 13 to 3, or 1 with R2; 9 of the 42 *Nanoglomus* calls inside unnamed cultures become U1) and changes no v4 design or LongGloDB call. R1b, its block variant, gives its first estimates: no change on v4 or LongGloDB; on EUKARYOME false known 13 to 3 (1 with R2) and kept 152 to 150 (148 with R2), every block coming from one of three LSU units held only by an unnamed culture. R5, the *Rhizoglomus* motif step, reproduces its EUKARYOME estimates (66 of 95 *Rhizoglomus* carriers called divergent by VT, 93 with the SSU fallback; 1 of 438 other carriers reclassified) and recognises more short-5.8S copies from their own sequence: under species holdout 29 of 43 frozen, 32 with R5 and 40 with the fallback, without changing any main-type call. R3, the provisional-genus flag G1, gives its first estimates: R3b flags 40 of 69 wrong genera in the 15 genus rounds with at most 1 false flag in 346 right-genus calls; R3a on any marker flags all 69 but also 19 of 487 right-genus calls in the species holdout (all *Septoglomus viscosum*, on SSU). `--rule r5s,r1a,r2,r3b` chains rules in protocol order (R5 or R5s, R1a or R1b, R2, then one R3 arm). `--score` reports each design's outcomes, frozen against arm, with the protocol's paired interval (95% cluster bootstrap over species, 10,000 replicates, fixed seed), and tests pin all eight R2 estimates. The other rules follow |
| `name`: name query sequences with the paper's method | Reads the four windows from `amfsplit regions` (ssu_flank, full_its, s58, LSU_LROR-FLR2), searches VT and SH with BLAST+ against the db, places the LSU windows with MAFFT and EPA-ng, and names each query as its own held-out culture with species clades plus N6, as the paper's external tests did. Writes `calls.tsv` (the call, its claims and reasons), `rows.tsv` (every column of the frozen row that does not need the query's true culture), `skipped.tsv` and `run.json`; BLAST and EPA-ng results are cached in `work/`. On the 280 LongGloDB test windows it reproduces the published condition A primary calls in every shared column |
| `db`: databases and tools for `name` | Checks MaarjAM, UNITE 10.0 and the frozen V18 reference folder against their pins, builds the VT and SH BLAST databases once (UNITE filtered to Glomeromycota and renamed as the frozen loader does), copies the V18 reference, and records each tool's path, version and hash in `manifest.json` |

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
