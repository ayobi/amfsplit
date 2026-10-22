#!/usr/bin/env python3
"""Golden per-copy tables in amf-novelty-benchmark that the Rust port must reproduce for Gate 0b
(the paper's method: LSU species clades with N6), with row counts and SHA-256 of the decompressed
content. Identifiers and derived numbers only; LongGloDB tables carry no sequences.

Usage: python3 v3_gate0b_goldens.py --repo /path/to/amf-novelty-benchmark [--out v3_gate0b_goldens.tsv]
"""
import argparse, csv, glob, gzip, hashlib, os, sys

GOLDENS = [
    ("species clades, culture rounds (148)", "benchmarks/m4_external_holdout/frozen_m2b/comparison/epa_species_clades_physical_culture.copies.tsv"),
    ("species clades, species rounds (50)", "benchmarks/m4_external_holdout/frozen_m2b/comparison/epa_species_clades_physical_species.copies.tsv"),
    ("LSU removal, V18 removed, without N6", "benchmarks/m4_external_holdout/aggregate_v2/species.v1.copies.tsv"),
    ("LSU removal, V18 removed, with N6", "benchmarks/m4_external_holdout/aggregate_v2/species.v2.copies.tsv"),
    ("LSU removal, V18 intact, without N6", "benchmarks/m4_external_holdout/aggregate_v2/species.control_v1.copies.tsv"),
    ("LSU removal, V18 intact, with N6", "benchmarks/m4_external_holdout/aggregate_v2/species.control_v2.copies.tsv"),
]

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--repo", required=True)
    ap.add_argument("--out", default="v3_gate0b_goldens.tsv")
    a = ap.parse_args()
    items = list(GOLDENS)
    for run, label in (("runs_euk", "EUKARYOME"), ("runs_lgdb", "LongGloDB")):
        for p in sorted(glob.glob(os.path.join(a.repo, "benchmarks", "m4_external_holdout", run, "06_scores", "*.copies.tsv*"))):
            items.append((f"{label}, {os.path.basename(p).split('.copies')[0]}", os.path.relpath(p, a.repo)))
    rows = []
    for what, rel in items:
        path = os.path.join(a.repo, rel)
        op = gzip.open if path.endswith(".gz") else open
        with op(path, "rb") as fh:
            data = fh.read()
        n = data.count(b"\n") - 1
        rows.append({"golden": what, "path": rel, "rows": n, "sha256_content": hashlib.sha256(data).hexdigest()})
        print(f"{rows[-1]['sha256_content'][:12]}  {n:5d}  {what}")
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("Not in the public release: per-copy N6 tables for the culture and species rounds (in the amfsplit "
          "repository's golden set) and the 15 genus rounds (export from the local M4 run).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
