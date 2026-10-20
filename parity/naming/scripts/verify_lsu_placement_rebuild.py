#!/usr/bin/env python3
"""Verify a fresh frozen-M2b LSU placement build against archived hashes.

Uses only the Python standard library. It never downloads data and never changes
project files. On mismatches it writes placement_row_differences.tsv in the run
folder for targeted review; differences are not automatically waved through.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

EXPECTED_REFERENCE_COMMIT = "f4a0014336f49be6aea48e8374dae052360a88bd"
EXPECTED_SCRIPT_SHA256 = "359f297a4c7038327fc94edea19ec4b8dab9f91a723b64b7c5a09ca4810610cc"
EXPECTED_AMF_ONLY_SHA256 = "e6794cbbe166735a94e44b146aae2ecffb35ebfee61cf6ae8ed2e48ba9227198"
EXPECTED_QUERY_COUNT = 903
EXPECTED_REFERENCE_TIPS = 495
EXPECTED_TREE_EDGES = 987
EXPECTED_MASS = 0.95

# Hashes copied from the locked FROZEN_M2B.json manifest.
EXPECTED_ARTIFACTS = {
    "reference/reference.unaligned.fasta": "874308b806b60ddb4ddf917091893a2138c3538bda304d95703ae3b139c48d66",
    "reference/reference.aligned.fasta": "1cd6426557f842d9e0085b95f02f3156406d43fb71f840a2c6891bece769ce33",
    "reference/reference_metadata.tsv": "5c38fd7fddd3896492644503359917d5fbd1ac3d4c818c0b4c7fc8a49ff5591d",
    "reference/backbone.topology.newick": "a6ae1cb91d1ca04021028d066a4b30f2a77c1235d13cfcfdafa9598eb5f0fa66",
    "reference/fit.raxml.bestTree": "d3da113b83e3d185ebfe10fd1958aa121e401c0dd035c74350ef607a74abf11b",
    "reference/fit.raxml.bestModel": "2e7a4dd93fd75c60586bad52a2fd4e4c5a73043289e89ad87f686df21e2a6194",
    "queries.fasta": "18679043a974d1bff4fdf798c0f74655d7444e1b0916a5c0c8356c796a9cca2f",
    "queries.aligned.fasta": "c8f5d1e41e049aad46e9197396cfe24ec6accaa61deafa38a56c930fe739755f",
    "summary/edge_annotations.json": "d54eca9e489977de8c6dead886c48cf83028a4e61cdb67a766e7a7b67bd25dd6",
    "summary/lsu_placements.copies.tsv": "1de103b6ae45f0372dda9ab0799230c4d1243cb98add943e5bab67785e5e910d",
}

FLOAT_COLUMNS = {
    "retained_fraction", "best_edge_lwr", "best_pendant_length", "best_distal_length",
    "lwr_sum", "lwr_entropy_nats", "credible_mass_target", "credible_mass_retained",
    "amf_edge_mass", "external_species_clade_mass", "external_genus_clade_mass",
}
INT_COLUMNS = {
    "aligned_bases", "raw_bases", "deleted_insertions", "ambiguous_aligned_bases",
    "bridge_alignment_eligible", "best_edge", "credible_edge_count",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def tsv_diagnostics(golden: Path, actual: Path, out_dir: Path) -> tuple[bool, list[str]]:
    messages: list[str] = []
    expected_header, expected_rows = read_tsv(golden)
    actual_header, actual_rows = read_tsv(actual)
    if expected_header != actual_header:
        messages.append(f"TSV header differs: expected {expected_header!r}, got {actual_header!r}")
        return False, messages
    if len(expected_rows) != EXPECTED_QUERY_COUNT or len(actual_rows) != EXPECTED_QUERY_COUNT:
        messages.append(f"Row count mismatch: golden={len(expected_rows)}, actual={len(actual_rows)}, expected={EXPECTED_QUERY_COUNT}")

    exp_by_id = {r.get("accession", ""): r for r in expected_rows}
    act_by_id = {r.get("accession", ""): r for r in actual_rows}
    if len(exp_by_id) != len(expected_rows) or len(act_by_id) != len(actual_rows):
        messages.append("Duplicate accession(s) found in one or both TSVs")
    missing = sorted(set(exp_by_id) - set(act_by_id))
    extra = sorted(set(act_by_id) - set(exp_by_id))
    if missing:
        messages.append(f"Missing accessions ({len(missing)}): {', '.join(missing[:12])}")
    if extra:
        messages.append(f"Unexpected accessions ({len(extra)}): {', '.join(extra[:12])}")

    diffs: list[dict[str, str]] = []
    col_count: dict[str, int] = {}
    max_abs: dict[str, tuple[float, str]] = {}
    for acc in sorted(set(exp_by_id) & set(act_by_id)):
        e, a = exp_by_id[acc], act_by_id[acc]
        for col in expected_header:
            ev, av = e.get(col, ""), a.get(col, "")
            if ev == av:
                continue
            col_count[col] = col_count.get(col, 0) + 1
            delta = ""
            if col in FLOAT_COLUMNS and ev and av:
                try:
                    d = abs(float(av) - float(ev))
                    delta = format(d, ".17g")
                    prev = max_abs.get(col)
                    if prev is None or d > prev[0]:
                        max_abs[col] = (d, acc)
                except ValueError:
                    pass
            diffs.append({"accession": acc, "column": col, "golden": ev, "actual": av, "absolute_numeric_delta": delta})

    diff_path = out_dir / "placement_row_differences.tsv"
    if diffs:
        with diff_path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=["accession", "column", "golden", "actual", "absolute_numeric_delta"], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(diffs)
        messages.append(f"TSV cell differences: {len(diffs)} across {len(col_count)} columns; details: {diff_path}")
        for col, n in sorted(col_count.items(), key=lambda item: (-item[1], item[0])):
            messages.append(f"  {col}: {n} differing cells")
        for col, (d, acc) in sorted(max_abs.items()):
            messages.append(f"  {col}: max absolute delta {d:.12g} at {acc}")
        cat_cols = [c for c in col_count if c not in FLOAT_COLUMNS]
        if cat_cols:
            messages.append("Categorical/integer differences need explicit review; no cutoff-tie tolerance is auto-accepted.")
        else:
            messages.append("Only numeric cells differ; this is still NOT an automatic pass. Review the differences against the locked protocol.")
        return False, messages
    messages.append("Per-cell comparison: exact (all 903 rows and every column match)")
    return True, messages


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project", type=Path, default=Path.cwd(), help="AMFSplit repository root (default: current directory)")
    ap.add_argument("--run", type=Path, required=True, help="Fresh output directory passed to lsu_placement.py build")
    a = ap.parse_args()
    project = a.project.resolve()
    run = a.run.resolve()
    failures = 0
    print(f"Project: {project}\nRebuild: {run}\n")

    run_json_path = run / "run.json"
    try:
        info = json.loads(run_json_path.read_text(encoding="utf-8"))
    except Exception as exc:
        print(f"FAIL run.json: {exc}")
        return 2
    if info.get("reference_commit") != EXPECTED_REFERENCE_COMMIT:
        print(f"FAIL reference commit: expected {EXPECTED_REFERENCE_COMMIT}, got {info.get('reference_commit')!r}")
        failures += 1
    else:
        print(f"PASS reference commit: {EXPECTED_REFERENCE_COMMIT}")
    if info.get("script_sha256") != EXPECTED_SCRIPT_SHA256:
        print(f"FAIL frozen script SHA-256: expected {EXPECTED_SCRIPT_SHA256}, got {info.get('script_sha256')!r}")
        failures += 1
    else:
        print("PASS frozen placement script SHA-256")
    input_hashes = info.get("reference_sha256", {})
    actual_amf_sha = input_hashes.get("V18_LSUDB_052025_AMFONLY.fasta")
    if actual_amf_sha != EXPECTED_AMF_ONLY_SHA256:
        print(f"FAIL AMF-only V18 SHA-256: expected {EXPECTED_AMF_ONLY_SHA256}, got {actual_amf_sha!r}")
        failures += 1
    else:
        print("PASS AMF-only V18 SHA-256")

    for rel, expected in EXPECTED_ARTIFACTS.items():
        path = run / rel
        if not path.is_file():
            print(f"FAIL missing artifact: {rel}")
            failures += 1
            continue
        got = sha256(path)
        if got != expected:
            print(f"FAIL {rel}\n  expected {expected}\n  actual   {got}")
            failures += 1
        else:
            print(f"PASS {rel} ({got})")

    summary = info.get("summary", {})
    expected_summary = {
        "queries": EXPECTED_QUERY_COUNT,
        "reference_tips": EXPECTED_REFERENCE_TIPS,
        "tree_edges": EXPECTED_TREE_EDGES,
        "credible_mass": EXPECTED_MASS,
    }
    for key, expected in expected_summary.items():
        got = summary.get(key)
        ok = got == expected
        if isinstance(expected, float) and isinstance(got, (float, int)):
            ok = math.isclose(float(got), expected, rel_tol=0, abs_tol=1e-15)
        if ok:
            print(f"PASS summary.{key}: {got}")
        else:
            print(f"FAIL summary.{key}: expected {expected!r}, got {got!r}")
            failures += 1
    if summary.get("alignment_eligible") != EXPECTED_QUERY_COUNT:
        print(f"FAIL summary.alignment_eligible: expected 903, got {summary.get('alignment_eligible')!r}")
        failures += 1
    else:
        print("PASS summary.alignment_eligible: 903/903")

    golden = project / "parity/naming/placement/lsu_placements.copies.tsv"
    actual = run / "summary/lsu_placements.copies.tsv"
    if golden.is_file() and actual.is_file():
        exact_rows, messages = tsv_diagnostics(golden, actual, run)
        for msg in messages:
            print(("PASS " if msg.startswith("Per-cell comparison") else "INFO ") + msg)
        if not exact_rows:
            failures += 1
    else:
        print(f"FAIL golden or actual summary TSV missing: golden={golden.exists()} actual={actual.exists()}")
        failures += 1

    if failures:
        print(f"\nPLACEMENT REBUILD: FAIL ({failures} failed checks). Do not waive mismatches automatically.")
        return 1
    print("\nPLACEMENT REBUILD: PASS (all pinned artifacts and all 903 per-query summary rows match exactly)")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
