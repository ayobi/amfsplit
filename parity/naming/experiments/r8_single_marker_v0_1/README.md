# AMFSplit R8 exploratory experiment (formerly N6)

> **Renamed, 1 November 2026,** from `n6_unplaced_v0_1`. Protocol v3 calls this rule R8, the single-marker novelty flag. N6 now means only the ITS distance rule, which is part of the paper's method and checked by `amfsplit gate0`. The files below keep their original wording.

This directory contains a read-only audit of the frozen M2b outputs, an experimental N6 evaluator, unit tests, and the result note. It does not modify or replace any frozen baseline files.

Run the tests from this directory:

```bash
python3 -m unittest -v test_unplaced_novelty_rule.py
```

See `N6_MULTI_UNMATCHED_v0.1_REPORT.md` for scope, observed results, limitations, and reproduction commands. `audit_novelty_calls.py` reads the baseline TSVs from an `amf-novelty-benchmark` checkout. `evaluate_unplaced_novelty_rule.py` consumes the audit TSV and writes results to a new output directory.
