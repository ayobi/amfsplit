#!/usr/bin/env python3
"""Golden outputs of the Stage 1 prototype (amfsplit.py v0.3.1) for the Rust port's parity tests.

Run once, with Python 3 and edlib 1.3.9 (the package whose C++ core the Rust port vendors):

    python3 -I parity/stage1/make_golden.py            # from the repository root

It writes
  inputs/synthetic/   deterministic synthetic reads, anchors and panels (built here, seed 7)
  golden/<case>/      the prototype's records.tsv and summary.json for each case, plus
                      MANIFEST.tsv with the SHA-256 and size of every output file (FASTAs included)
  cases.tsv           one row per case: reads, anchors, panel and the exact command-line options

and refuses to run if the prototype script is not the frozen v0.3.1 (checked by SHA-256).
The real inputs in inputs/stefani903/ are copied, not generated: the 903 retained sequences of
culture reference v4 (amf-novelty-benchmark, reference/culture_ref_v4/fasta/full.fasta) and the
ITSxRust anchors of the 904 deposited sequences (benchmarks/stefani/out/itsx/stefani.anchors.tsv).
PX215117.1, dropped from the reference, is not in this workspace; with --only-anchored its
anchors row is unused, and the 904-sequence check runs where GenBank is reachable
(parity/stage1/README.md).

The synthetic reads extend tests/test_v03.py of the prototype (same construction, same seed for
the shared reads) with minus-strand, strandless, lowercase/N, partial-chain, short and duplicate
reads, a FASTQ.gz copy and a panel with region overrides and comment lines. Primers are random,
so these cases test logic, not primer biology.
"""
import csv
import gzip
import hashlib
import json
import random
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PROTO = HERE / "python" / "amfsplit_v0_3_1.py"
PROTO_SHA256 = "a7d79a637095894aea6c8c68b13944a0204d187889185590ad3df44ca0aa4a05"
REAL_PANEL = ROOT / "panels" / "amf_legacy_windows.tsv"
SYN = HERE / "inputs" / "synthetic"
STEF = HERE / "inputs" / "stefani903"
GOLD = HERE / "golden"

# ----------------------------------------------------------------------------- synthetic reads
# Construction copied from amfsplit v0.3.1 tests/test_v03.py (same rng, same order of draws), so
# reads A to N are the reads of that test.
rng = random.Random(7)
COMP = str.maketrans("ACGT", "TGCA")


def rnd(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def rc(s):
    return s.translate(COMP)[::-1]


def mutate(s, n_edits):
    s = list(s)
    for i in rng.sample(range(len(s)), n_edits):
        s[i] = rng.choice([b for b in "ACGT" if b != s[i]])
    return "".join(s)


P = {name: rnd(20) for name in ("AML1", "NS31", "AMV45", "AM1", "AMDGR", "LROR", "FLR2", "FLR3", "FLR4", "LSUMBR")}
PANEL = [
    ("AMPLICON_AML1-LSUmBr", "AML1", P["AML1"], "LSUmBr", P["LSUMBR"], 1800, 3500),
    ("SSU_VT_NS31-AM1", "NS31", P["NS31"], "AM1", P["AM1"], 400, 650),
    ("SSU_V4_AMV4.5NF-AMDGR", "AMV4.5NF", P["AMV45"], "AMDGR", P["AMDGR"], 150, 300),
    ("LSU_LROR-FLR2", "LROR", P["LROR"], "FLR2", P["FLR2"], 500, 900),
    ("LSU_D2_FLR3-FLR4", "FLR3", P["FLR3"], "FLR4", P["FLR4"], 250, 450),
]


def put(seq, pos, site):
    return seq[:pos] + site + seq[pos + len(site):]


def make_copy():
    ssu, its, lsu = rnd(1600), rnd(465), rnd(1100)
    ssu = put(ssu, 20, P["AML1"])
    ssu = put(ssu, 1000, P["NS31"])
    ssu = put(ssu, 1100, P["AMV45"])
    ssu = put(ssu, 1300, rc(P["AMDGR"]))
    ssu = put(ssu, 1500, rc(P["AM1"]))
    lsu = put(lsu, 40, P["LROR"])
    lsu = put(lsu, 500, P["FLR3"])
    lsu = put(lsu, 780, rc(P["FLR2"]))
    lsu = put(lsu, 900, rc(P["FLR4"]))
    lsu = put(lsu, 1050, rc(P["LSUMBR"]))
    return ssu, its, lsu


def anchors_row(rid, n_ssu=1600, n_its=465, strand="+"):
    a, b = n_ssu + 1, n_ssu + n_its
    return {"read_id": rid, "its1_start": a, "its1_end": a + 108, "its2_start": a + 265, "its2_end": b,
            "full_start": a, "full_end": b, "confidence": "confident", "ssu_strand": strand, "s58s_strand": strand,
            "s58e_strand": strand, "lsu_strand": strand}


def build():
    ssu, its, lsu = make_copy()
    full = ssu + its + lsu
    reads, anchors = [], {}

    def add(name, seq, anchored=True, anc=None):
        reads.append((name, seq))
        if anchored:
            anchors[name] = anc or anchors_row(name)

    # --- reads of tests/test_v03.py, in its order of construction (the rng draws must match)
    add("A_clean", full)
    add("B_am1_absent", put(full, 1500, mutate(rc(P["AM1"]), 9)))
    add("C_am1_absent_decoy", put(put(full, 1500, mutate(rc(P["AM1"]), 9)), 1600 + 465 + 600, rc(P["AM1"])))
    cut = 1010
    add("D_nested_5p", full[cut:], anc=anchors_row("D_nested_5p", 1600 - cut))
    add("E_trunc_3p", full[:1600 + 465 + 700])
    add("F_lror_absent", put(full, 1600 + 465 + 40, mutate(P["LROR"], 9)))
    add("G_concat_tandem", full + mutate(full, 30), anc=anchors_row("G_concat_tandem"))
    add("H_concat_inverted", full + rc(full), anc=anchors_row("H_concat_inverted"))
    add("I_unanchored", full, anchored=False)
    add("J_offtarget", rnd(3200), anchored=False)
    noisy = put(full, 150, mutate(P["NS31"], 3))
    noisy = put(noisy, 1900, mutate(P["FLR3"], 3))
    noisy = put(put(noisy, 1300, mutate(P["AML1"], 3)), 2300, mutate(P["AML1"], 3))
    noisy = put(noisy, 400, mutate(rc(P["FLR2"]), 3))
    add("K_spurious_repeats", noisy)
    ssu_, its_, lsu_ = full[:1600], full[1600:2065], full[2065:]
    other = mutate(ssu_, 80) + rnd(len(its_)) + mutate(lsu_, 88)
    add("N_concat_different_molecule", full + other, anc=anchors_row("N_concat_different_molecule"))
    tail = full[-330:]
    add("L_partial_tandem", full + tail, anc=anchors_row("L_partial_tandem"))
    add("M_partial_foldback", full + rc(tail), anc=anchors_row("M_partial_foldback"))

    # --- additions for the port (drawn after the reads above, so those stay identical)
    # minus-strand read, anchors in rRNA-sense coordinates as ITSxRust reports them
    add("O_minus_strand", rc(mutate(full, 12)), anc=anchors_row("O_minus_strand", strand="-"))
    # anchored, but no strand in the anchors row: orientation from the primers, regions from the row
    row = anchors_row("P_no_strand")
    for key in ("ssu_strand", "s58s_strand", "s58e_strand", "lsu_strand"):
        row[key] = ""
    add("P_no_strand", rc(full), anc=row)
    # lowercase input with N in two primer sites (N equals any base, the primer's IUPAC codes do not)
    q = put(full, 1003, "NN")
    q = put(q, 1600 + 465 + 45, "N")
    add("Q_lowercase_N", q.lower())
    # unanchored minus-strand read: orientation from the primers
    add("S_minus_unanchored", rc(full), anchored=False)
    # anchored, empty confidence
    row = anchors_row("T_empty_confidence")
    row["confidence"] = ""
    add("T_empty_confidence", full, anc=row)
    # partial chain: full ITS known, ITS2 end missing (no 5.8S, ITS2 or ITS2-derived regions)
    row = anchors_row("U_partial_chain")
    row["its2_start"] = ""
    row["its2_end"] = ""
    add("U_partial_chain", full, anc=row)
    # anchors cut the read at base 1: no SSU flank
    add("V_ssu_at_start", full[1600:], anc=anchors_row("V_ssu_at_start", n_ssu=0))
    # read far shorter than any primer region, and one shorter than a primer
    add("W_short_read", full[990:1060], anchored=False)
    add("X_tiny_read", "ACGTACGTAC", anchored=False)
    # duplicate read ID: both records are processed, the anchors row is shared
    add("A_clean", full)
    # a read whose two LSU primers are swapped in order (reverse primer before forward)
    swapped = put(put(full, 1600 + 465 + 40, rc(P["FLR2"])), 1600 + 465 + 780, P["LROR"])
    add("Y_out_of_order", swapped)
    # SSU V4 window too short (forward primer moved 50 nt downstream)
    tshort = put(put(full, 1100, rnd(20)), 1150, P["AMV45"])
    add("Z_too_short", tshort)
    return reads, anchors


ANCHOR_COLS = ["read_id", "its1_start", "its1_end", "its2_start", "its2_end", "full_start", "full_end",
               "confidence", "ssu_strand", "s58s_strand", "s58e_strand", "lsu_strand"]


def write_synthetic():
    reads, anchors = build()
    SYN.mkdir(parents=True, exist_ok=True)
    with open(SYN / "reads.fa", "w") as fh:
        for k, v in reads:
            # wrapped at 70 columns, as many FASTA writers do
            fh.write(f">{k} synthetic read\n" + "\n".join(v[i:i + 70] for i in range(0, len(v), 70)) + "\n")
    with gzip.open(SYN / "reads.fq.gz", "wt", compresslevel=6) as fh:
        for k, v in reads:
            fh.write(f"@{k} synthetic\n{v}\n+\n{'I' * len(v)}\n")
    with open(SYN / "anchors.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ANCHOR_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(anchors.values())
    with open(SYN / "panel.tsv", "w") as fh:
        fh.write("window\ttarget\tfwd_name\tfwd_seq\trev_name\trev_seq\tmin_len\tmax_len\n")
        for row in PANEL:
            fh.write("\t".join(map(str, row[:1] + ("test",) + row[1:])) + "\n")
    # panel with comments, blank lines and explicit region columns (one override, one left empty)
    with open(SYN / "panel_regions.tsv", "w") as fh:
        fh.write("# synthetic panel with region columns\n\n")
        fh.write("window\ttarget\tfwd_name\tfwd_seq\trev_name\trev_seq\tmin_len\tmax_len\tfwd_region\trev_region\n")
        for i, row in enumerate(PANEL):
            fr, rr = ("any", "any") if i == 0 else (("", "") if i == 1 else ("", "lsu" if row[0].startswith("LSU") else ""))
            fh.write("\t".join(map(str, row[:1] + ("test",) + row[1:])) + f"\t{fr}\t{rr}\n")
            if i == 2:
                fh.write("# comment between windows\n")
    # an empty read, alone (the prototype crashes on it in the concatemer screen; tested without it)
    with open(SYN / "empty.fa", "w") as fh:
        fh.write(">E_empty\n\n>A_clean\n" + reads[0][1] + "\n")


# ----------------------------------------------------------------------------- cases
CASES = [
    # name, reads, anchors, panel, options
    ("syn_default", "inputs/synthetic/reads.fa", "inputs/synthetic/anchors.tsv", "inputs/synthetic/panel.tsv",
     "--open-ends --max-err-frac 0.15"),
    ("syn_keep_concatemers", "inputs/synthetic/reads.fa", "inputs/synthetic/anchors.tsv", "inputs/synthetic/panel.tsv",
     "--open-ends --keep-concatemers"),
    ("syn_v02_like", "inputs/synthetic/reads.fa", "inputs/synthetic/anchors.tsv", "inputs/synthetic/panel.tsv",
     "--open-ends --no-region-limit --no-concat-screen"),
    ("syn_keep_primers_closed", "inputs/synthetic/reads.fa", "inputs/synthetic/anchors.tsv", "inputs/synthetic/panel.tsv",
     "--keep-primers"),
    ("syn_only_anchored_margin0", "inputs/synthetic/reads.fa", "inputs/synthetic/anchors.tsv", "inputs/synthetic/panel.tsv",
     "--open-ends --only-anchored --region-margin 0"),
    ("syn_fastq_concat_params", "inputs/synthetic/reads.fq.gz", "inputs/synthetic/anchors.tsv", "inputs/synthetic/panel.tsv",
     "--open-ends --max-err-frac 0.2 --concat-min-primers 3 --concat-strong 2 --concat-err-frac 0.1 "
     "--concat-sim 0.3 --concat-tol 0.1 --concat-min-sep 400"),
    ("syn_panel_regions", "inputs/synthetic/reads.fa", "inputs/synthetic/anchors.tsv", "inputs/synthetic/panel_regions.tsv",
     "--open-ends --region-margin 120"),
    ("syn_no_anchors", "inputs/synthetic/reads.fa", "", "inputs/synthetic/panel.tsv",
     "--open-ends"),
    ("syn_empty_read", "inputs/synthetic/empty.fa", "", "inputs/synthetic/panel.tsv",
     "--open-ends --no-concat-screen"),
    ("stefani903_hifi", "inputs/stefani903/reads.fasta.gz", "inputs/stefani903/anchors.tsv.gz", "../../panels/amf_legacy_windows.tsv",
     "--max-err-frac 0.15 --only-anchored --open-ends"),
    ("stefani903_ont", "inputs/stefani903/reads.fasta.gz", "inputs/stefani903/anchors.tsv.gz", "../../panels/amf_legacy_windows.tsv",
     "--max-err-frac 0.20 --only-anchored --open-ends"),
    ("stefani903_keep_primers", "inputs/stefani903/reads.fasta.gz", "inputs/stefani903/anchors.tsv.gz", "../../panels/amf_legacy_windows.tsv",
     "--keep-primers --keep-concatemers"),
]
STORE_FASTA = lambda case: case.startswith("syn_")  # real-data FASTAs are checked by hash only


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def run_case(name, reads, anchors, panel, options, tmp):
    """The prototype reads the anchors with open(), so gzipped anchors are decompressed first."""
    a = ""
    if anchors:
        a = str(HERE / anchors)
        if a.endswith(".gz"):
            plain = tmp / f"{name}.anchors.tsv"
            with gzip.open(a, "rb") as src, open(plain, "wb") as dst:
                shutil.copyfileobj(src, dst)
            a = str(plain)
    out = tmp / name / "out"
    out.parent.mkdir(parents=True)
    cmd = [sys.executable, "-I", str(PROTO), "-i", str(HERE / reads), "-p", str((HERE / panel).resolve()),
           "-o", str(out)] + (["-a", a] if a else []) + options.split()
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"{name}: prototype failed\n{r.stderr[-2000:]}")
    dest = GOLD / name
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    rows = []
    for f in sorted(out.parent.iterdir()):
        suffix = f.name[len("out."):]
        rows.append((suffix, sha256(f), f.stat().st_size))
        if suffix in ("records.tsv", "summary.json") or STORE_FASTA(name):
            shutil.copyfile(f, dest / suffix)
    with open(dest / "MANIFEST.tsv", "w") as fh:
        fh.write("file\tsha256\tbytes\n")
        for suffix, digest, size in rows:
            fh.write(f"{suffix}\t{digest}\t{size}\n")
    summ = json.load(open(dest / "summary.json"))
    return summ["reads"], summ["reads_processed"], summ["concatemers"]["flagged"]


def main():
    import edlib
    if sha256(PROTO) != PROTO_SHA256:
        sys.exit(f"{PROTO} is not the frozen v0.3.1 prototype")
    write_synthetic()
    tmp = HERE / ".golden_tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir()
    with open(HERE / "cases.tsv", "w") as fh:
        fh.write("case\treads\tanchors\tpanel\toptions\n")
        for name, reads, anchors, panel, options in CASES:
            n, done, conc = run_case(name, reads, anchors, panel, options, tmp)
            fh.write(f"{name}\t{reads}\t{anchors}\t{panel}\t{options}\n")
            print(f"{name:<28} reads {n:>4}  processed {done:>4}  concatemers {conc}")
    shutil.rmtree(tmp)
    prov = {"prototype": PROTO.name, "prototype_sha256": PROTO_SHA256,
            "python": sys.version.split()[0], "edlib": edlib.__version__ if hasattr(edlib, "__version__") else "",
            "inputs": {str(p.relative_to(HERE)): sha256(p) for p in sorted((HERE / "inputs").rglob("*")) if p.is_file()},
            "panel": {"panels/amf_legacy_windows.tsv": sha256(REAL_PANEL)}}
    try:
        import importlib.metadata as md
        prov["edlib"] = md.version("edlib")
    except Exception:
        pass
    with open(GOLD / "PROVENANCE.json", "w") as fh:
        json.dump(prov, fh, indent=1)
        fh.write("\n")


if __name__ == "__main__":
    main()
