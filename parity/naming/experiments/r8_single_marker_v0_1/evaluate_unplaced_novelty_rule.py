#!/usr/bin/env python3
"""Run a small, read-only ablation of unplaced AMF novelty calls.

Input is the joined TSV produced by audit_novelty_calls.py. The frozen source
status is never overwritten: separate experimental_status / experimental_claim
columns are appended, with a summary of changed copies, cultures and species.

Rule N6_MULTI_UNMATCHED_v0.1 is a conservative evidence gate, not a claim of
biological novelty: only a legacy 'novel lineage' call with reconciliation rank
'unplaced', no matched arms, and at least two explicitly unmatched markers keeps
its experimental novel-lineage-candidate status. Other unplaced patterns become
'insufficient evidence'. All other legacy statuses are preserved.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

RULE_ID = "N6_MULTI_UNMATCHED_v0.1"
MARKERS = ("VT", "SH", "LSU")
REQUIRED = {
    "mode", "accession", "culture", "organism", "VT_state", "SH_state",
    "LSU_state", "matched_count", "unmatched_count", "no_call_count",
    "rank", "status", "novelty_claim",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header: {path}")
        missing = REQUIRED - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Input is not a novelty decision audit TSV; missing columns: {sorted(missing)}")
        return list(reader)


def experimental_decision(row: dict[str, str]) -> tuple[str, str, str]:
    """Return (experimental status, experimental novelty claim, explanation)."""
    baseline_status = row["status"]
    baseline_claim = row["novelty_claim"]
    if baseline_status != "novel lineage" or row["rank"] != "unplaced":
        return baseline_status, baseline_claim, "N6_NOT_APPLICABLE_BASELINE_PRESERVED"

    matched = int(row["matched_count"])
    unmatched = int(row["unmatched_count"])
    no_call = int(row["no_call_count"])
    if matched > 0:
        return (
            "insufficient evidence", "not_asserted",
            "N6_CONFLICT_UNPLACED_BUT_MATCHED_ARM_PRESENT",
        )
    if unmatched >= 2:
        return (
            "novel lineage", "candidate_only",
            "N6_AT_LEAST_TWO_UNMATCHED_MARKERS_NO_MATCHED_ARMS",
        )
    return (
        "insufficient evidence", "not_asserted",
        f"N6_TOO_FEW_UNMATCHED_MARKERS_unmatched={unmatched}_no_call={no_call}",
    )


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def counts(rows: list[dict[str, str]], key: str) -> dict[str, int]:
    return dict(sorted(Counter(r.get(key, "") or "(blank)" for r in rows).items(), key=lambda x: (-x[1], x[0])))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--audit", type=Path, required=True, help="novelty_decision_audit.tsv from the read-only audit")
    ap.add_argument("--out", type=Path, required=True, help="new output directory; must not already exist")
    a = ap.parse_args()
    if a.out.exists():
        ap.error(f"Output directory already exists; refusing to overwrite: {a.out}")

    source_rows = read_tsv(a.audit)
    by_mode: dict[str, list[dict[str, str]]] = {}
    seen: set[tuple[str, str]] = set()
    for row in source_rows:
        key = (row["mode"], row["accession"])
        if key in seen:
            raise ValueError(f"Duplicate mode/accession: {key}")
        seen.add(key)
        if row["VT_state"] not in {"matched", "unmatched", "no call"}:
            raise ValueError(f"Unexpected VT state {row['VT_state']!r} for {row['accession']}")
        if row["SH_state"] not in {"matched", "unmatched", "no call"}:
            raise ValueError(f"Unexpected SH state {row['SH_state']!r} for {row['accession']}")
        if row["LSU_state"] not in {"matched", "unmatched", "no call"}:
            raise ValueError(f"Unexpected LSU state {row['LSU_state']!r} for {row['accession']}")
        for k in ("matched_count", "unmatched_count", "no_call_count"):
            row[k] = str(int(row[k]))
        if sum(int(row[k]) for k in ("matched_count", "unmatched_count", "no_call_count")) != 3:
            raise ValueError(f"Marker counts do not sum to 3 for {row['accession']}")
        by_mode.setdefault(row["mode"], []).append(row)

    output: list[dict[str, str]] = []
    changes: list[dict[str, str]] = []
    summary: dict[str, object] = {
        "rule_id": RULE_ID,
        "input": str(a.audit.resolve()),
        "interpretation": (
            "Development-only ablation. Does not modify frozen source tables. It gates only existing novel-lineage calls: "
            "an unplaced call must have no matched arms and at least two explicitly unmatched markers to retain candidate-only status. "
            "This does not establish biological novelty, and these already-inspected holdout outputs are not independent validation."
        ),
        "modes": {},
    }

    for mode, rows in sorted(by_mode.items()):
        if len(rows) != 903 or len({r["accession"] for r in rows}) != 903:
            raise ValueError(f"Expected 903 unique copies for {mode}; observed {len(rows)} rows")
        exp_rows: list[dict[str, str]] = []
        for src in rows:
            r = dict(src)
            ex_status, ex_claim, reason = experimental_decision(r)
            r["baseline_status"] = src["status"]
            r["baseline_novelty_claim"] = src["novelty_claim"]
            r["experimental_status"] = ex_status
            r["experimental_novelty_claim"] = ex_claim
            r["experimental_rule_id"] = RULE_ID if reason.startswith("N6_") and reason != "N6_NOT_APPLICABLE_BASELINE_PRESERVED" else "BASELINE_PRESERVED"
            r["experimental_reason"] = reason
            exp_rows.append(r)
            output.append(r)
            if src["status"] != ex_status or src["novelty_claim"] != ex_claim:
                changes.append(r)

        unplaced = [r for r in exp_rows if r["rank"] == "unplaced" and r["baseline_status"] == "novel lineage"]
        new_insufficient = [r for r in changes if r["mode"] == mode and r["experimental_status"] == "insufficient evidence"]
        summary["modes"][mode] = {
            "rows": len(exp_rows),
            "baseline_status_counts": counts(exp_rows, "baseline_status"),
            "experimental_status_counts": counts(exp_rows, "experimental_status"),
            "baseline_candidate_claims": sum(r["baseline_novelty_claim"] == "candidate_only" for r in exp_rows),
            "experimental_candidate_claims": sum(r["experimental_novelty_claim"] == "candidate_only" for r in exp_rows),
            "legacy_unplaced_novel_lineage": len(unplaced),
            "unplaced_state_patterns": dict(sorted(Counter("/".join(r[f"{m}_state"] for m in MARKERS) for r in unplaced).items())),
            "changed_copy_count": sum(r in changes for r in exp_rows),
            "changed_accessions": [r["accession"] for r in new_insufficient],
            "changed_cultures": sorted({r["culture"] for r in new_insufficient}),
            "changed_organisms": sorted({r["organism"] for r in new_insufficient}),
        }

    if not by_mode:
        raise ValueError("No rows found")
    a.out.mkdir(parents=True)
    columns = list(source_rows[0].keys()) + [
        "baseline_status", "baseline_novelty_claim", "experimental_status",
        "experimental_novelty_claim", "experimental_rule_id", "experimental_reason",
    ]
    write_tsv(a.out / "novelty_rule_v0_1.tsv", output, columns)
    write_tsv(a.out / "changed_decisions.tsv", changes, columns)
    (a.out / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"Rule: {RULE_ID}")
    print("Scope: unplaced baseline 'novel lineage' calls only; all other statuses preserved.")
    for mode, data in summary["modes"].items():
        print(f"\n=== {mode} ===")
        print(f"Rows: {data['rows']}; baseline candidate claims: {data['baseline_candidate_claims']}; experimental candidate claims: {data['experimental_candidate_claims']}")
        print(f"Unplaced novel-lineage rows: {data['legacy_unplaced_novel_lineage']}; changed copies: {data['changed_copy_count']}")
        if data["changed_accessions"]:
            print("Changed accessions: " + ", ".join(data["changed_accessions"]))
            print("Affected cultures: " + "; ".join(data["changed_cultures"]))
            print("Affected named species: " + "; ".join(data["changed_organisms"]))
    print(f"\nWrote: {a.out / 'novelty_rule_v0_1.tsv'}")
    print(f"Wrote: {a.out / 'changed_decisions.tsv'}")
    print(f"Wrote: {a.out / 'summary.json'}")
    print("Frozen baseline inputs and code were read only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
