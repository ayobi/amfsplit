#!/usr/bin/env python3
"""Golden outputs of the Stage 1 validation and triage prototypes (validate_windows.py and triage.py,
amfsplit v0.3.1) for the Rust port's parity tests.

    python3 -I parity/triage/make_golden.py --v18 V18_LSUDB_052025.fasta --vsearch /path/to/vsearch

Needs Python 3 with edlib, vsearch 2.27.0 and the Delavaux V18 LSU database *with* its outgroups
(V18_LSUDB_052025.fasta, SHA-256 28a806a6...; github.com/c383d893/AMF-LSU-Database-and-Pipeline2,
commit f4a0014). Neither vsearch nor V18 is stored here.

Inputs (written to inputs/, gzipped):
  windows.*.fasta   the prototype's Stage 1 outputs on the 903 culture sequences (stage1 golden case
                    stefani903_hifi), with test windows appended to the LR0R-FLR2 FASTA:
                      MUT07_* / MUT12_*  culture windows with 7% / 12% substitutions (seed 11)
                      RND_*              random sequences (no hit at 70% identity)
                      OUT_<V18 id>       nine non-AMF records of V18, copied exactly; stored with
                                         their sequence masked to N (V18 is not redistributed) and
                                         restored from V18 when vsearch runs
  records.tsv       Stage 1 records plus a row per test window and five rows with no window (NOWIN_*)
  hits.tsv          validate_windows.py hits of the LR0R-FLR2 windows against V18 (the triage input)
Golden outputs (golden/<case>/): every file each prototype writes, and its stdout; MANIFEST.tsv
holds the SHA-256 of each. The AMF-only FASTAs are kept as hashes only. The validate case runs on
the unmasked windows (vsearch needs the real sequences); the triage cases run on the stored, masked
inputs, which are what the parity test uses (triage only copies sequences through).

Note on the default outgroups: triage.py's default (--outgroups Mortierella,Rhodotorula,Exophiala,
Citrus) matches an older Delavaux release. V18 holds 35 non-AMF records, and only Mortierella is in
that list, so under the default the OUT_* windows other than Mortierella are called AMF. The golden
files record the prototype as it is; the port reproduces it.
"""
import argparse
import csv
import gzip
import hashlib
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
VALIDATE = HERE / "python" / "validate_windows_v0_3_1.py"
TRIAGE = HERE / "python" / "triage_v0_3_1.py"
STAGE1 = ROOT / "parity" / "stage1" / "python" / "amfsplit_v0_3_1.py"
SHA = {
    VALIDATE: "ced7f223eb0a33263e4ebcb0e6ba050a9c73a1628dafb22c99cc6b8eba171448",
    TRIAGE: "8a233d8eeba8d1dbe445f0294fa41c15358266d28ff89a1761ea1c769faa6b83",
    STAGE1: "a7d79a637095894aea6c8c68b13944a0204d187889185590ad3df44ca0aa4a05",
}
V18_SHA = "28a806a6b1930e74e413b83d5c6cc5c5d71268a740ed38488d0bb25de9a494a7"
OUTGROUP_IDS = [  # non-AMF records of V18 used as test windows
    "DQ273794.1_Mortierella_verticillata_Zygo",
    "AY198398.1_Candida_glabrata_Asco",
    "DQ273811.1_Endogone_pisiformis_Zygo",
    "AY213626.1_Rhizopus_oryzae_Zygo",
    "U13369.1_Homo_sapien",
    "X52322.1_Arabidopsis_thaliana",
    "FM876840_Geosiphon_pyriformis",
    "AM183920_Geosiphon_pyriformis",
    "Z19136.1_Schizosaccharomyces_pombe_Asco",
]
LSU = "LSU_LROR-FLR2"
VALIDATE_OPTS = ["--id", "0.70", "--threads", "1"]
CASES = [
    # name, triage options ({hits}, {out}, {records}, {prefix} are filled in)
    ("default", "--hits {hits} --out {out} --records {records} --fasta-prefix {prefix}"),
    ("amf_pattern", "--hits {hits} --out {out} --records {records} --fasta-prefix {prefix} --amf-pattern _[A-Z][a-z]+aceae_"),
    ("ont_keep_unresolved", "--hits {hits} --out {out} --records {records} --fasta-prefix {prefix} --min-id 85 --keep-unresolved"),
    ("hits_only", "--hits {hits} --out {out}"),
    ("custom_outgroups", "--hits {hits} --out {out} --records {records} --fasta-prefix {prefix} "
                         "--outgroups Mortierella,Candida,Homo,Arabidopsis,Rhizopus,Endogone,Schizosaccharomyces"),
    ("records_no_fasta", "--hits {hits} --out {out} --records {records}"),
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_fasta(path):
    out, cur = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            cur = line[1:].split()[0]
            out[cur] = []
        elif cur:
            out[cur].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def mutate(rng, s, frac):
    s = list(s)
    for i in rng.sample(range(len(s)), int(len(s) * frac)):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


def build_inputs(v18, work):
    """Stage 1 on the 903 culture sequences, then the test windows. Returns the window prefix."""
    s1 = ROOT / "parity" / "stage1"
    anchors = work / "anchors.tsv"
    with gzip.open(s1 / "inputs/stefani903/anchors.tsv.gz", "rb") as a, open(anchors, "wb") as b:
        shutil.copyfileobj(a, b)
    prefix = work / "windows"
    subprocess.run([sys.executable, "-I", str(STAGE1), "-i", str(s1 / "inputs/stefani903/reads.fasta.gz"),
                    "-a", str(anchors), "-p", str(ROOT / "panels/amf_legacy_windows.tsv"), "-o", str(prefix),
                    "--max-err-frac", "0.15", "--only-anchored", "--open-ends"], check=True, capture_output=True)
    lsu = read_fasta(f"{prefix}.{LSU}.fasta")
    rng = random.Random(11)
    ids = sorted(lsu)
    extra = []
    for acc in rng.sample(ids, 30):
        extra.append((f"MUT07_{acc}", mutate(rng, lsu[acc], 0.07)))
    for acc in rng.sample(ids, 60):
        extra.append((f"MUT12_{acc}", mutate(rng, lsu[acc], 0.12)))
    for i in range(10):
        extra.append((f"RND_{i:02d}", "".join(rng.choice("ACGT") for _ in range(700))))
    ref = read_fasta(v18)
    for oid in OUTGROUP_IDS:
        extra.append((f"OUT_{oid}", ref[oid].upper()))
    with open(f"{prefix}.{LSU}.fasta", "a") as fh:
        for rid, seq in extra:
            fh.write(f">{rid} window={LSU} test=yes\n{seq}\n")
    # records: Stage 1 rows, then one row per test window and five reads without a window
    with open(f"{prefix}.records.tsv") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    with open(f"{prefix}.records.tsv", "a") as fh:
        pad = "\t" * (len(header) - 1)
        for rid, _ in extra:
            fh.write(f"{rid}{pad}\n")
        for i in range(5):
            fh.write(f"NOWIN_{i}{pad}\n")
    for f in work.iterdir():
        if f.name.startswith("windows.") and f.name.endswith(".json"):
            f.unlink()
    return prefix


def store_inputs(prefix, hits):
    """Gzip the inputs into inputs/, with the V18 outgroup sequences masked."""
    dest = HERE / "inputs"
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir()
    for f in sorted(prefix.parent.iterdir()):
        if not f.name.startswith(prefix.name + "."):
            continue
        suffix = f.name[len(prefix.name) + 1:]
        if not (suffix.endswith(".fasta") or suffix == "records.tsv"):
            continue
        with open(f) as src, gzip.open(dest / f"windows.{suffix}.gz", "wt", compresslevel=9) as dst:
            masked = False
            for line in src:
                if line.startswith(">"):
                    masked = line[1:].startswith("OUT_")
                    dst.write(line)
                else:
                    dst.write("N" * len(line.rstrip("\n")) + "\n" if masked else line)
    shutil.copyfile(hits, dest / "hits.tsv")


def run(cmd, env, stdout_path=None):
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode:
        sys.exit(f"failed: {' '.join(cmd)}\n{r.stderr[-2000:]}")
    if stdout_path:
        Path(stdout_path).write_text(r.stdout)
    return r


def record_case(name, files, keep_full):
    dest = HERE / "golden" / name
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    with open(dest / "MANIFEST.tsv", "w") as fh:
        fh.write("file\tsha256\tbytes\n")
        for label, path in sorted(files.items()):
            fh.write(f"{label}\t{sha256(path)}\t{os.path.getsize(path)}\n")
            if keep_full(label):
                shutil.copyfile(path, dest / label)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--v18", required=True)
    ap.add_argument("--vsearch", required=True)
    a = ap.parse_args()
    for p, digest in SHA.items():
        if sha256(p) != digest:
            sys.exit(f"{p} is not the frozen prototype")
    if sha256(a.v18) != V18_SHA:
        sys.exit("--v18 is not V18_LSUDB_052025.fasta (SHA-256 28a806a6...)")
    env = dict(os.environ, PATH=str(Path(a.vsearch).resolve().parent) + os.pathsep + os.environ["PATH"])
    version = subprocess.run([a.vsearch, "--version"], capture_output=True, text=True).stderr.splitlines()[0]

    work = HERE / ".golden_tmp"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir()
    prefix = build_inputs(a.v18, work)

    # validate_windows.py: the vsearch run that feeds triage
    vout = work / "validate"
    run([sys.executable, "-I", str(VALIDATE), "-w", f"{prefix}.{LSU}.fasta", "-r", a.v18, "-o", str(vout)]
        + VALIDATE_OPTS, env, work / "validate.stdout")
    record_case("validate", {"hits.tsv": f"{vout}.hits.tsv", "summary.json": f"{vout}.summary.json",
                             "stdout.txt": work / "validate.stdout"}, lambda f: True)
    store_inputs(prefix, f"{vout}.hits.tsv")

    rows = []
    for name, opts in CASES:
        case_dir = work / name
        case_dir.mkdir()
        # a fresh copy of the stored (masked) windows: triage writes next to them, and the parity
        # test runs on exactly these files
        for f in (HERE / "inputs").glob("windows.*.gz"):
            with gzip.open(f, "rb") as src, open(case_dir / f.name[:-3], "wb") as dst:
                shutil.copyfileobj(src, dst)
        p = case_dir / "windows"
        out = case_dir / "out"
        args = opts.format(hits=HERE / "inputs" / "hits.tsv", out=out, records=f"{p}.records.tsv", prefix=p).split()
        run([sys.executable, "-I", str(TRIAGE)] + args, env, case_dir / "stdout.txt")
        files = {"stdout.txt": case_dir / "stdout.txt"}
        for f in case_dir.iterdir():
            if f.name.startswith("out."):
                files[f.name[len("out."):]] = f
            elif f.name.endswith(".amf.fasta"):
                files[f.name] = f
        record_case(f"triage_{name}", files, lambda f: not f.endswith(".amf.fasta"))
        rows.append((f"triage_{name}", opts))
    with open(HERE / "cases.tsv", "w") as fh:
        fh.write("case\toptions\n")
        for name, opts in rows:
            fh.write(f"{name}\t{opts}\n")
    with open(HERE / "golden" / "PROVENANCE.txt", "w") as fh:
        fh.write(f"validate_windows.py {SHA[VALIDATE]}\ntriage.py {SHA[TRIAGE]}\nvsearch {version}\n"
                 f"V18_LSUDB_052025.fasta {V18_SHA}\nvalidate options {' '.join(VALIDATE_OPTS)}\n"
                 f"python {sys.version.split()[0]}\n")
    shutil.rmtree(work)
    print(f"validate + {len(CASES)} triage cases written; {version}")


if __name__ == "__main__":
    main()
