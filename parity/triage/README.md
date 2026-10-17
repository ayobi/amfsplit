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
cargo test --release --test parity_triage
AMFSPLIT_TEST_VSEARCH=/path/to/vsearch AMFSPLIT_TEST_V18=/path/to/V18_LSUDB_052025.fasta \
    cargo test --release --test parity_triage      # adds the vsearch case
```

Without those two variables, the validate case prints `SKIPPED` and passes. All seven cases
passed here on 17 October 2026, with vsearch v2.27.0 (linux x86_64 release) and V18 at commit
`f4a0014` (SHA-256 `28a806a6…`).

## Findings while porting

- **The default outgroups do not fit V18.** `triage.py`'s default, `--outgroups
  Mortierella,Rhodotorula,Exophiala,Citrus`, lists the outgroups of an older Delavaux release
  (`STAGE1.md` names `LRORFLR2.V10seqs_8.4.20.fasta`). V18 has 35 non-AMF records, and only
  Mortierella is on that list. In `triage_default`, 8 of the 9 non-AMF windows are therefore
  called AMF; `triage_amf_pattern` calls all 9 non-AMF.
  - The paper is unaffected: its tests used curated culture sequences, and the stored Stefani
    triage ran on the older release (1 of its 136 hit targets is a V18 name).
  - The port keeps the prototype's behaviour. When this happens, `amfsplit triage` prints a note
    on stderr. Whether to change the default is an open decision.
- **The vsearch hits order depends on threads.** With more than one thread, vsearch writes the
  same hit lines in another query order. Best hits, and so triage calls, are unchanged, because
  each query's hits stay together. The bytes of `hits.tsv`, and possibly the tie order of
  `top_targets` in the validation summary, are not. Parity is therefore checked with one thread.
- **Geosiphon.** *Geosiphon pyriformis* is a glomeromycete, but V18 keeps it out of its AMF-only
  file. The V18 pattern therefore calls it non-AMF.
