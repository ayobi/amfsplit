# AMFSplit Gate 0: reference reproduction and Rust wiring status

**Date:** 19 October 2026

## Verified frozen evidence and reference

The uploaded evidence ZIP has exactly four members (`MANIFEST.json` and the VT, SH and LSU compressed evidence files). The four SHA-256 hashes match the frozen M2b lock exactly. The manifest identifies reference fingerprint `181938d695f3df40`.

All **16** outputs recorded in the frozen culture-reference manifest match their pinned SHA-256 hashes. The M2b source checksum lock matches its pinned hash and every source file listed in that lock matches. Each marker evidence payload has 903 queries and the expected arm/region/threshold settings.

### Evidence payloads

| File | Bytes | SHA-256 | Check |
|---|---:|---|---|
| `VT.evidence.json.gz` | 402,365 | `5bcd06453d32f8718f573368009b12d19b55c58c2712e439408d26e1a9217f56` | PASS |
| `SH.evidence.json.gz` | 1,500,239 | `de21b3fa16fe8df6ca43421260acfd912391ef513e6bf5ead86853be4b9e81f6` | PASS |
| `LSU.evidence.json.gz` | 276,076 | `1a259082f41ef44cdfca2d02b64dc7ef17c34473b693db56a292f5d35ecc9c4b` | PASS |

## Frozen Python reference reproduction

The original frozen M2b Python comparison was re-run with a small compatibility shim around the project’s **same vendored Edlib C++ source**. The shim covers only the `edlib.align` distance calls used in these scripts and was sanity-checked with known alignment examples. It is a reference-engine run, not a Rust compilation.

| Holdout | Output variant | Rows | Byte-identical to saved golden? | SHA-256 |
|---|---|---:|---|---|
| Physical-culture (148 folds) | Frozen baseline | 903 | PASS | `6cc5412acf147ba7644d767851872517f78551b372e614c3f07297227d539f34` |
| Physical-culture (148 folds) | LSU species-clade comparison | 903 | PASS | `a357b563a5def31932085dabe4697f3b49e0f586d56a6c91acefe8bb1cb7cfd4` |
| Species (50 folds) | Frozen baseline | 903 | PASS | `906420d70e9d8d24ab5b115ecd63e72451c396b4f1e7c60cfccf7b8d6a8a0800` |
| Species (50 folds) | LSU species-clade comparison | 903 | PASS | `9b0d4ba422214a763cbf0c5949f8f7b2b163fdfcd8e2fd5ebf81a55db521a53a` |

The species folds were run in batches because the combined process hit the execution limit; the last ten folds were run separately. Their outputs were merged by accession and compared byte-for-byte with the saved tables. All four 903-row tables were exact matches.

## Rust wiring patch

The project snapshot already contains a substantial `crates/amfsplit-naming` implementation, the frozen evidence and source locks, and a unit test named `frozen_m2b_assignment_rows_match_every_golden_cell`. The crate was dormant: it was absent from the root Cargo workspace and disconnected from the CLI.

The attached wiring patch makes five small changes:
- Adds `crates/amfsplit-naming` to the Cargo workspace and lockfile.
- Adds it as a dependency of the `amfsplit` executable.
- Exposes `amfsplit gate0`, defaulting to the checkout’s `parity/naming` fixtures; `--fixtures DIR` can select another fixture root.
- Leaves the production `amfsplit name` command explicitly as future work; `gate0` is a reproducibility command, not a claim that the user-facing naming workflow is finished.
- Updates README status and run instructions.

## Important remaining limitation

**Rust Gate 0 has not yet been declared passed.** This environment does not have `cargo` or `rustc`, so I could not compile the newly wired crate or run its Rust test. What is verified here is the input bundle, frozen source/reference integrity, and byte-exact output reproducibility of the Python reference. The Rust gate still must pass on the workstation.

## Commands to run locally

From the root of your `amfsplit` checkout, save `amfsplit_gate0_wiring.patch` there and run:

```bash
git apply --check amfsplit_gate0_wiring.patch
git apply amfsplit_gate0_wiring.patch
cargo test -p amfsplit-naming frozen_m2b_assignment_rows_match_every_golden_cell -- --nocapture
cargo run --release -- gate0
cargo test --workspace --release
```

If `git apply --check` passes, the patch matches the uploaded source snapshot. The decisive Gate 0 result is that the Rust command reports both holdout modes byte-identical with zero differing cells.

## Data boundary

The gate uses the curated culture-reference bundle, numeric marker evidence, and locked per-copy output. It does not require external raw read databases. No LongGloDB or EUKARYOME evidence payloads or sequences were added to the patch or report.
