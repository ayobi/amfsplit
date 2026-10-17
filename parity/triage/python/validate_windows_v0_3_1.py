#!/usr/bin/env python3
"""Validate amfsplit windows against the legacy reference database they should match.

For each window sequence, report its best vsearch hit (global alignment, both strands):
identity, query coverage and target coverage. A window cut at the right primer sites
should cover its reference region almost end to end (high qcov and tcov) and match
the expected taxon at high identity.

Usage:
  validate_windows.py -w results/bc1026.LSU_LROR-FLR2.fasta \
      -r LRORFLR2.V10seqs_8.4.20.fasta -o results/bc1026.LSU_LROR-FLR2.vs_ref
"""
import argparse
import json
import statistics
import subprocess
import sys
from collections import Counter


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-w", "--windows", required=True, help="window FASTA from amfsplit")
    ap.add_argument("-r", "--ref", required=True, help="legacy reference FASTA")
    ap.add_argument("-o", "--out", required=True, help="output prefix")
    ap.add_argument("--id", type=float, default=0.70, help="minimum identity to report")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    hits_path = f"{args.out}.hits.tsv"
    cmd = ["vsearch", "--usearch_global", args.windows, "--db", args.ref,
           "--id", str(args.id), "--strand", "both", "--maxaccepts", "16", "--maxrejects", "64",
           "--userout", hits_path, "--userfields", "query+target+id+qcov+tcov+alnlen",
           "--threads", str(args.threads), "--quiet"]
    subprocess.run(cmd, check=True)

    n_query = sum(1 for line in open(args.windows) if line.startswith(">"))
    best = {}
    for line in open(hits_path):
        q, t, ident, qcov, tcov, alnlen = line.rstrip("\n").split("\t")
        ident, qcov, tcov = float(ident), float(qcov), float(tcov)
        if q not in best or ident > best[q]["id"]:
            best[q] = {"target": t, "id": ident, "qcov": qcov, "tcov": tcov}

    def med(key):
        vals = [b[key] for b in best.values()]
        return round(statistics.median(vals), 1) if vals else None

    top = Counter(b["target"] for b in best.values()).most_common(5)
    summary = {
        "windows": n_query,
        "with_hit": len(best),
        "with_hit_pct": round(100 * len(best) / max(1, n_query), 1),
        "best_identity_median": med("id"),
        "query_coverage_median": med("qcov"),
        "target_coverage_median": med("tcov"),
        "identity_ge_97_pct": round(100 * sum(b["id"] >= 97 for b in best.values()) / max(1, n_query), 1),
        "top_targets": top,
    }
    with open(f"{args.out}.summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)
    json.dump(summary, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
