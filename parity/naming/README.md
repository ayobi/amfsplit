# Frozen M2b naming parity fixture

This directory freezes the **culture-only** input bundle used to check the Rust naming port against the M2b Python reference implementation.

- `inputs/benchmark/ref/`: curated Stefani culture reference v4 (903 copies, 151 cultures), plus the 13 reference FASTAs regenerated from the verified public culture-read inputs. Each FASTA hash matches the frozen reference manifest.
- `inputs/benchmark/evidence/`: three numeric per-marker identity evidence files (VT, SH, LSU), and the manifest. The payload and manifest SHA-256 values match the frozen M2b release exactly.
- `inputs/benchmark/baseline/`: the v4 top-hit tables needed to preserve each arm's frozen tie choices.
- `placement/`: frozen per-copy LSU placement annotations used as descriptive output columns; baseline naming uses the frozen culture bridge and arm evidence.
- `golden/`: archived outputs for physical-culture and physical-species holdouts, for the nearest-reference LSU arm (`baseline_*`, Gate 0a) and the LSU species-clade arm (`epa_species_clades_*`, Gate 0b step 1, copied from amf-novelty-benchmark `benchmarks/m4_external_holdout/frozen_m2b/comparison/`). These contain per-copy assignments, not sequence reads.
- `frozen_sources/m2b/`: the checksum-recorded Python bridge, reconciliation, novelty, and placement/reporting scripts plus the frozen protocol. This is the reference implementation for auditing the Rust port.

LongGloDB and EUKARYOME evidence payloads and sequences are intentionally absent. No external reference FASTA files are bundled; the external matches are represented only by the frozen numeric evidence.

Golden-output SHA-256:
- `baseline_physical_culture.copies.tsv`: `6cc5412acf147ba7644d767851872517f78551b372e614c3f07297227d539f34`
- `baseline_physical_species.copies.tsv`: `906420d70e9d8d24ab5b115ecd63e72451c396b4f1e7c60cfccf7b8d6a8a0800`
- `epa_species_clades_physical_culture.copies.tsv`: `a357b563a5def31932085dabe4697f3b49e0f586d56a6c91acefe8bb1cb7cfd4`
- `epa_species_clades_physical_species.copies.tsv`: `9b0d4ba422214a763cbf0c5949f8f7b2b163fdfcd8e2fd5ebf81a55db521a53a`
