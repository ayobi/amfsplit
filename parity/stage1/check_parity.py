#!/usr/bin/env python3
"""Run the Rust `amfsplit regions` on every case in cases.tsv and compare with the golden outputs.

    python3 -I parity/stage1/check_parity.py [--bin target/release/amfsplit] [--keep DIR]

Every output file listed in golden/<case>/MANIFEST.tsv must match its SHA-256; for a mismatch
the records are compared field by field and the first differences printed. The Rust binary
also writes {out}.run.json (provenance), which the prototype does not; it is checked to exist
and otherwise ignored. Exit status 0 only if every case matches. Needs only the standard library.
"""
import argparse
import csv
import gzip
import hashlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def first_record_diffs(golden, mine, limit=8):
    """Field-level differences between two records.tsv files."""
    def load(p):
        with open(p, newline="") as fh:
            return list(csv.reader(fh, delimiter="\t"))
    g, m = load(golden), load(mine)
    out = []
    if g[0] != m[0]:
        out.append(f"header differs: {set(g[0]) ^ set(m[0]) or 'order'}")
    for i, (a, b) in enumerate(zip(g[1:], m[1:]), start=1):
        for col, x, y in zip(g[0], a, b):
            if x != y:
                out.append(f"row {i} ({a[0]}) {col}: golden {x!r} rust {y!r}")
                if len(out) >= limit:
                    return out
    if len(g) != len(m):
        out.append(f"row count: golden {len(g) - 1}, rust {len(m) - 1}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=str(ROOT / "target" / "release" / "amfsplit"))
    ap.add_argument("--keep", help="keep the Rust outputs in this directory")
    ap.add_argument("--case", action="append", help="run only these cases")
    a = ap.parse_args()
    work = Path(a.keep) if a.keep else Path(tempfile.mkdtemp(prefix="amfsplit-parity-"))
    work.mkdir(parents=True, exist_ok=True)
    failures = 0
    with open(HERE / "cases.tsv") as fh:
        cases = list(csv.DictReader(fh, delimiter="\t"))
    for c in cases:
        if a.case and c["case"] not in a.case:
            continue
        out = work / c["case"] / "out"
        out.parent.mkdir(parents=True, exist_ok=True)
        cmd = [a.bin, "regions", "-i", str(HERE / c["reads"]), "-p", str((HERE / c["panel"]).resolve()), "-o", str(out)]
        if c["anchors"]:
            cmd += ["-a", str(HERE / c["anchors"])]
        cmd += c["options"].split()
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode:
            print(f"FAIL {c['case']}: exit {r.returncode}\n{r.stderr[-1500:]}")
            failures += 1
            continue
        gold = HERE / "golden" / c["case"]
        with open(gold / "MANIFEST.tsv") as fh:
            manifest = list(csv.DictReader(fh, delimiter="\t"))
        bad = []
        for row in manifest:
            mine = Path(f"{out}.{row['file']}")
            if not mine.exists():
                bad.append(f"{row['file']}: missing")
            elif sha256(mine) != row["sha256"]:
                bad.append(f"{row['file']}: differs ({mine.stat().st_size} vs {row['bytes']} bytes)")
        produced = {p.name[len('out.'):] for p in out.parent.iterdir()}
        extra = produced - {row["file"] for row in manifest} - {"run.json"}
        if extra:
            bad.append(f"unexpected outputs: {sorted(extra)}")
        if "run.json" not in produced:
            bad.append("run.json missing")
        if bad:
            failures += 1
            print(f"FAIL {c['case']}: " + "; ".join(bad))
            if any(b.startswith("records.tsv") for b in bad):
                for d in first_record_diffs(gold / "records.tsv", f"{out}.records.tsv"):
                    print("     " + d)
        else:
            print(f"ok   {c['case']}  ({len(manifest)} files identical)")
    if not a.keep:
        shutil.rmtree(work)
    print(f"\n{failures} case(s) failed" if failures else "\nall cases identical")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
