#!/usr/bin/env python3
"""Read-only audit of frozen AMF novelty decisions against per-marker reconciliation.

Reads the committed M2b baseline TSVs. Writes audit tables to a new output directory;
never modifies source/benchmark files.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

MODES = {
    "physical_culture": ("reconcile_culture.copies.tsv", "baseline_consistent_culture.copies.tsv"),
    "physical_species": ("reconcile_species.copies.tsv", "baseline_consistent_species.copies.tsv"),
}
MARKERS = ("VT", "SH", "LSU")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No data rows found in {path}")
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in columns})


def counts(rows: list[dict[str, object]], keyfn) -> dict[str, int]:
    c = Counter(str(keyfn(row) or "(blank)") for row in rows)
    return dict(sorted(c.items(), key=lambda kv: (-kv[1], kv[0])))


def require_unique(rows: list[dict[str, str]], path: Path, label: str) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    for row in rows:
        accession = row.get("accession", "")
        if not accession:
            raise ValueError(f"Missing accession in {path}")
        if accession in out:
            raise ValueError(f"Duplicate accession {accession} in {path}")
        out[accession] = row
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--benchmark-root", type=Path, required=True,
                    help="root of amf-novelty-benchmark Git checkout")
    ap.add_argument("--out", type=Path, required=True,
                    help="new output directory outside the source checkout is recommended")
    args = ap.parse_args()

    in_dir = args.benchmark_root / "benchmarks/m4_external_holdout/frozen_m2b/inputs/benchmark/baseline"
    if not in_dir.is_dir():
        ap.error(f"Expected frozen M2b baseline directory not found: {in_dir}")
    if args.out.exists():
        ap.error(f"Output directory already exists; refusing to overwrite: {args.out}")
    args.out.mkdir(parents=True)

    all_audit: list[dict[str, object]] = []
    unplaced_audit: list[dict[str, object]] = []
    summary: dict[str, object] = {
        "source_directory": str(in_dir.resolve()),
        "interpretation": (
            "Descriptive audit of frozen reference-holdout outputs only. A marker 'unmatched' means its "
            "recorded unit has no supporting training-culture holder under this reconciliation; 'no call' "
            "means no usable unit was assigned by that arm. Neither state alone proves biological novelty."
        ),
        "modes": {},
    }

    for mode, (rec_name, nov_name) in MODES.items():
        rec_path, nov_path = in_dir / rec_name, in_dir / nov_name
        rec_rows, nov_rows = read_tsv(rec_path), read_tsv(nov_path)
        rec = require_unique(rec_rows, rec_path, "reconciliation")
        nov = require_unique(nov_rows, nov_path, "novelty")
        if set(rec) != set(nov):
            missing_novelty = sorted(set(rec) - set(nov))[:5]
            missing_reconcile = sorted(set(nov) - set(rec))[:5]
            raise ValueError(f"Accession mismatch for {mode}: only in reconcile={missing_novelty}; only in novelty={missing_reconcile}")

        joined: list[dict[str, object]] = []
        for accession in sorted(rec):
            r, n = rec[accession], nov[accession]
            state_signature = "/".join(r.get(f"{m}_state", "") for m in MARKERS)
            um = [m for m in MARKERS if m in (r.get("unmatched_arms", "").split(";"))]
            no_call = [m for m in MARKERS if r.get(f"{m}_state") == "no call"]
            matched = [m for m in MARKERS if r.get(f"{m}_state") == "matched"]
            row: dict[str, object] = {
                "mode": mode,
                "accession": accession,
                "culture": r.get("culture", n.get("culture", "")),
                "organism": r.get("organism", n.get("organism", "")),
                "copy_class": r.get("copy_class", n.get("copy_class", "")),
                "div_kind": r.get("div_kind", n.get("div_kind", "")),
                "group": r.get("group", n.get("group", "")),
                "VT_state": r.get("VT_state", ""),
                "SH_state": r.get("SH_state", ""),
                "LSU_state": r.get("LSU_state", ""),
                "state_signature": state_signature,
                "matched_count": len(matched),
                "unmatched_count": len(um),
                "no_call_count": len(no_call),
                "rank": r.get("rank", ""),
                "call_species": r.get("call_species", ""),
                "call_genera": r.get("call_genera", ""),
                "unmatched_arms": r.get("unmatched_arms", ""),
                "status": n.get("status", ""),
                "novelty_claim": n.get("novelty_claim", ""),
                "reasons": n.get("reasons", ""),
                "flags": n.get("flags", ""),
            }
            joined.append(row)
            all_audit.append(row)
            if row["rank"] == "unplaced":
                unplaced_audit.append(row)

        unplaced = [x for x in joined if x["rank"] == "unplaced"]
        patterns = counts(unplaced, lambda x: x["state_signature"])
        status_rank = counts(joined, lambda x: f'{x["status"]} | {x["rank"]}')
        group_status = counts(joined, lambda x: f'{x["group"]} | {x["status"]}')
        mode_summary: dict[str, object] = {
            "n_rows": len(joined),
            "unique_accessions": len({x["accession"] for x in joined}),
            "group_counts": counts(joined, lambda x: x["group"]),
            "reconciliation_rank_counts": counts(joined, lambda x: x["rank"]),
            "novelty_status_counts": counts(joined, lambda x: x["status"]),
            "marker_state_counts": {m: counts(joined, lambda x, mm=m: x[f"{mm}_state"]) for m in MARKERS},
            "unplaced_n": len(unplaced),
            "unplaced_status_counts": counts(unplaced, lambda x: x["status"]),
            "unplaced_group_counts": counts(unplaced, lambda x: x["group"]),
            "unplaced_state_signatures": patterns,
            "unplaced_unmatched_count_distribution": counts(unplaced, lambda x: x["unmatched_count"]),
            "unplaced_no_call_count_distribution": counts(unplaced, lambda x: x["no_call_count"]),
            "unplaced_all_three_no_call": sum(x["no_call_count"] == 3 for x in unplaced),
            "unplaced_at_least_one_matched": sum(x["matched_count"] >= 1 for x in unplaced),
            "unplaced_one_unmatched_two_no_call": sum(x["unmatched_count"] == 1 and x["no_call_count"] == 2 for x in unplaced),
            "status_by_rank": status_rank,
            "group_by_status": group_status,
        }
        summary["modes"][mode] = mode_summary  # type: ignore[index]
        print(f"\n=== {mode} ===")
        print(f"Rows: {len(joined)}; unplaced: {len(unplaced)}; all-three-no-call among unplaced: {mode_summary['unplaced_all_three_no_call']}")
        print("Marker states:")
        for marker in MARKERS:
            print(f"  {marker}: {mode_summary['marker_state_counts'][marker]}")
        print("Unplaced marker-state patterns:")
        for pattern, n in patterns.items():
            print(f"  {pattern}: {n}")
        print("Unplaced novelty statuses:", mode_summary["unplaced_status_counts"])

    audit_columns = [
        "mode", "accession", "culture", "organism", "copy_class", "div_kind", "group",
        "VT_state", "SH_state", "LSU_state", "state_signature", "matched_count", "unmatched_count",
        "no_call_count", "rank", "call_species", "call_genera", "unmatched_arms", "status",
        "novelty_claim", "reasons", "flags",
    ]
    write_tsv(args.out / "novelty_decision_audit.tsv", all_audit, audit_columns)
    write_tsv(args.out / "unplaced_evidence_audit.tsv", unplaced_audit, audit_columns)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"\nWrote: {args.out / 'novelty_decision_audit.tsv'}")
    print(f"Wrote: {args.out / 'unplaced_evidence_audit.tsv'}")
    print(f"Wrote: {args.out / 'summary.json'}")
    print("Source tables were read only; no benchmark or frozen code was modified.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
