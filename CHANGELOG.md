# Changelog

## Unreleased

- Bioconda: vsearch 2.32.0 is a run dependency, for `validate`, and the recipe's test runs `validate`. The
  paper's hits were made with vsearch 2.27.0, which cannot be installed next to BLAST+ 2.17.0; 2.32.0 gives
  byte-identical outputs on the validate parity case (checked 5 November 2026).

## 0.1.0 (November 2026)

The paper's method, as in 0.1.0-rc.1, plus the fixes proposed since, as opt-in rules for Protocol v3 to
test. `crates/amfsplit-naming/src` and `crates/amfsplit/src` are byte for byte the arm code Protocol v3
locked on 4 November 2026 (tag protocol-v3-lock.1, doi:10.5281/zenodo.23962473). Every rule stays
opt-in until the protocol's confirmation run decides which become defaults.

- `arm`: Protocol v3's rules, with every rule's development estimates pinned in tests: R1a and R1b
  (unnamed cultures as taxa, inside and block), R2 (the inside gate), R3 (the provisional-genus flag:
  R3a, with its all-markers sensitivity, and R3b), R4 (within-species limits refitted on two collections,
  at the 97.5th percentile by default under decision 6), R5 (the *Rhizoglomus* motif step, by VT and with
  the SSU fallback), R6 (triage labels instead of filtering), R7 (the genus-clade fallback), R8 (the flag on
  novelty resting on one marker, and R8d) and S, the candidate default (R5 variant, R4, R1a, R2, R6, R7,
  R8, R3 variant). Rules chain in protocol order. `--score` compares each design with the published method
  using the protocol's paired cluster-bootstrap intervals, with the strata pooled over both external tests
  (conditions A and W), and `arm` takes several windows and placement tables in one run.
- `name --vt-spanned` (V): VT coverage over the part of the reference a sequence spans; reproduces the
  LongGloDB secondary analysis.
- `withhold`: an LSU-removal round rebuilt from scratch, as the frozen rounds were: Bio.Phylo's Newick
  prune and writer (pinned to its hashes on all 51 folds that remove V18 records), MAFFT L-INS-i,
  raxml-ng `--evaluate`, MAFFT `--addfragments`, EPA-ng and the placement summary. `--compare-with`
  checks each stage against an archived round, `--name` names the round's queries (stage 06), and
  `--reuse-reference` places other sequences in a finished round. Rebuilt this way, all 51 rounds give
  the published calls; three placement tables differ in the twelfth decimal of four sequences'
  likelihood weights.
- `withheld-k`: condition W, each K sequence of an external test named with its own species withheld, in
  its species' round; reproduces the published W tables (8 of 164 LongGloDB and 10 of 199 EUKARYOME K
  sequences called known).
- `adopt`: Protocol v3's section 8 applied rule by rule, S last, with its report; the dress rehearsal's
  decisions are pinned in tests.
- CI fetches amf-novelty-benchmark at 397a55b (release 0.3.2's rounds plus Protocol v3's inputs).

## 0.1.0-rc.1 (November 2026)

The paper's method (O'Brien, Marín & Tedersoo; benchmark release 0.3.1), from long reads to calls.

- `regions`, `validate`, `triage`: Stage 1, byte-identical to amfsplit.py 0.3.1.
- `evidence`: VT, SH and LSU evidence with BLAST+ against MaarjAM, UNITE 10.0 (Glomeromycota) and V18.
- `place`: MAFFT `--addfragments` and EPA-ng on the frozen V18 reference, and the placement summary.
- `db`: checks the databases against their pins, builds the BLAST databases and rebuilds the V18
  reference alignment from the release.
- `name`: names sequences with species clades plus the ITS distance rule (N6), from the four windows or
  straight from reads (ITSxRust anchors, then `regions`).
- `gate0`, `gate3`: reproducibility gates. Every per-call table of release 0.3.1 reproduces byte for
  byte (culture and species holdouts, the 52 LSU-removal and genus rounds, LongGloDB and EUKARYOME),
  and so do the paper's calls from the tool's own evidence and placements.

Known: EPA-ng lists tied placements in run-dependent order, so `best_edge` can differ between runs
when placements tie; no call depends on it, and `place` checks this against the frozen table.
