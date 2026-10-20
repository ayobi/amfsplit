# AMFSplit N6 exploratory experiment

This directory contains a read-only audit of the frozen M2b outputs, an experimental N6 evaluator, unit tests, and the result note. It does not modify or replace any frozen baseline files.

Run the tests from this directory:

```bash
python3 -m unittest -v test_unplaced_novelty_rule.py
```

See `N6_MULTI_UNMATCHED_v0.1_REPORT.md` for scope, observed results, limitations, and reproduction commands. `audit_novelty_calls.py` reads the baseline TSVs from an `amf-novelty-benchmark` checkout. `evaluate_unplaced_novelty_rule.py` consumes the audit TSV and writes results to a new output directory.
