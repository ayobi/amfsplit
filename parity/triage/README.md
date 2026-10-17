# Triage parity

`amfsplit validate` and `amfsplit triage` are ports of the Stage 1 prototypes
`validate_windows.py` and `triage.py` (amfsplit v0.3.1, in `python/`, unchanged; SHA-256
`ced7f223…` and `8a233d8e…`). Their outputs must be byte-identical:

- **validate:** `{out}.hits.tsv`, `{out}.summary.json` and stdout.
- **triage:** `{out}.triage.tsv`, `{out}.amf.ids`, `{out}.triage.summary.json`, stdout, and an
  `.amf.fasta` beside each Stage 1 FASTA.

## Inputs

`inputs/` is built by `make_golden.py` and gzipped. It starts from the prototype's Stage 1
outputs on the 903 culture sequences (`stefani903_hifi`), with 109 test windows appended to the
LR0R–FLR2 FASTA:

| Test windows | n | What they exercise |
|---|---|---|
| `MUT07_*` | 30 | Culture windows with 7% substitutions (AMF at 90%) |
| `MUT12_*` | 60 | Culture windows with 12% substitutions (mostly unresolved at 90%, AMF at 85%) |
| `RND_*` | 10 | Random sequence: a window, but no hit (no_hit) |
| `OUT_<V18 id>` | 9 | Non-AMF records of V18, copied exactly: Mortierella, Candida, Endogone, Rhizopus, Homo, Arabidopsis, Geosiphon (2), Schizosaccharomyces |

`records.tsv` adds one row per test window, plus five reads with no window (`NOWIN_*`).
`hits.tsv` holds the vsearch hits of `validate_windows.py` against V18 with outgroups
(`--id 0.70 --threads 1`).

V18 is not redistributed, so the `OUT_*` sequences are stored masked with N. Triage only copies
sequences, so the triage cases run on the masked files. The validate case restores the sequences
from V18 before running vsearch.

## Cases

| Case | Options |
|---|---|
| `validate` | `--id 0.70 --threads 1` (needs vsearch 2.27.0 and V18) |
| `triage_default` | `--records --fasta-prefix`, default outgroups, `--min-id 90` |
| `triage_amf_pattern` | `--amf-pattern _[A-Z][a-z]+aceae_` |
| `triage_ont_keep_unresolved` | `--min-id 85 --keep-unresolved` |
| `triage_hits_only` | hits only: no records, no FASTA prefix |
| `triage_custom_outgroups` | seven of the nine non-AMF genera as `--outgroups` |
| `triage_records_no_fasta` | records, no FASTA prefix (unhit windows become no_window, as in the prototype) |

## Running

```bash
cargo test --release --test parity_triage --test triage_amf_reference
# with the vsearch case and the real V18 AMF-only file:
AMFSPLIT_TEST_VSEARCH=/path/to/vsearch \
AMFSPLIT_TEST_V18=/path/to/V18_LSUDB_052025.fasta \
AMFSPLIT_TEST_V18_AMF=/path/to/V18_LSUDB_052025_AMFONLY.fasta \
    cargo test --release --test parity_triage --test triage_amf_reference
```

Without those variables, the two cases that need them print `SKIPPED` and pass. The test checks
the V18 files by SHA-256 and reports the vsearch version if a hit differs. Everything passed here on
17 October 2026, with vsearch v2.27.0 (linux x86_64 release) and V18 at commit `f4a0014`
(V18_LSUDB_052025.fasta `28a806a6…`, V18_LSUDB_052025_AMFONLY.fasta `e6794cbb…`).

## Where the port differs from the prototype: `--amf-reference`

The port has no default AMF rule (decided 17 October 2026). The caller gives one of three:

- `--amf-reference FASTA`, recommended: a target is AMF if and only if its ID is a record of
  this file. With Delavaux V18, search `V18_LSUDB_052025.fasta` and give
  `V18_LSUDB_052025_AMFONLY.fasta` here.
- `--amf-pattern REGEX`, as in the prototype.
- `--outgroups GENERA`, as in the prototype. Its v0.3.1 default
  (`Mortierella,Rhodotorula,Exophiala,Citrus`) now has to be given explicitly; the parity test
  does that.

`--amf-reference` cannot be combined with the other two. The summary's `params` records the file's
path, SHA-256 and record count; with the other two rules the summary is the prototype's.

Membership is strict: a target missing from the AMF file is non-AMF. Removing one record (the
most common best hit) moves 111 test reads to non_AMF. So the search reference and the AMF file
must come from the same release. The tool checks two cases:

- If no hit target is in the AMF file, triage stops with an error.
- If every hit target is in it, so no read could be called non_AMF, triage prints a note. That is
  the sign of a search against the AMF-only file itself.

The planned `db` step will store the pair together.

Tests (`tests/triage_amf_reference.rs`):

- On V18, the AMF-only file and the pattern `_[A-Z][a-z]+aceae_` select the same 460 of 495
  records. `--amf-reference` therefore has to reproduce the golden case `triage_amf_pattern`
  call for call, AMF-only FASTAs included. This is checked with a reference derived from the hit
  targets, and with the real V18 file when `AMFSPLIT_TEST_V18_AMF` is set.
- The other tests cover the error when no rule is given, the error on a reference from another
  release, the note when every target is AMF, and the option conflicts.

## Findings while porting

- **The default outgroups do not fit V18.** `triage.py`'s default, `--outgroups
  Mortierella,Rhodotorula,Exophiala,Citrus`, lists the outgroups of an older Delavaux release
  (`STAGE1.md` names `LRORFLR2.V10seqs_8.4.20.fasta`). V18 has 35 non-AMF records, and only
  Mortierella is on that list. In `triage_default`, 8 of the 9 non-AMF windows are therefore
  called AMF; `triage_amf_pattern` calls all 9 non-AMF.
  - The paper is unaffected: its tests used curated culture sequences, and the stored Stefani
    triage ran on the older release (1 of its 136 hit targets is a V18 name).
  - Decided: the port has no default rule, and `--amf-reference` is recommended (section above).
    With `--outgroups`, `amfsplit triage` still prints a note when a read is called AMF on a
    target that V18 labels non-AMF.
- **The vsearch hits order depends on threads.** With more than one thread, vsearch writes the
  same hit lines in another query order. Best hits, and so triage calls, are unchanged, because
  each query's hits stay together. The bytes of `hits.tsv`, and possibly the tie order of
  `top_targets` in the validation summary, are not. Parity is therefore checked with one thread.
- **Geosiphon.** *Geosiphon pyriformis* is a glomeromycete, but V18 keeps it out of its AMF-only
  file. The V18 pattern therefore calls it non-AMF.
