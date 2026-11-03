# Changelog

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
