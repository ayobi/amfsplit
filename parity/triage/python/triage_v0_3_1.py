#!/usr/bin/env python3
"""triage.py: decide which reads are Glomeromycota before any AMF classification.

ITSxRust anchors tell you where the regions are, not whose DNA it is: Mortierella-like
and other non-AMF fungi anchor perfectly well. This step labels every read from the best
hit of one of its windows (normally LSU_LROR-FLR2) against a reference that contains
non-AMF lineages, then writes AMF-only copies of the window FASTAs.

Two ways to say what counts as AMF:
  --outgroups  genera in the target header that are NOT AMF (Delavaux-style headers
               ACC_Genus_species; default: the four outgroups of the Delavaux LSU DB)
  --amf-pattern  regex that marks an AMF target (e.g. 'Glomeromycota' for a reference
               whose headers carry the lineage, such as a EUKARYOME LSU release)

Calls per read:
  AMF          best hit is AMF at >= --min-id
  unresolved   best hit is AMF but below --min-id (divergent AMF or a lineage the
               reference lacks: candidate dark matter, kept out of the AMF set by default)
  non_AMF      best hit is a non-AMF target
  no_hit       window present but no hit at the validation identity floor
  no_window    read has no usable window (only with --records)
"""
import argparse
import csv
import glob
import json
import re
from collections import Counter

DELAVAUX_OUTGROUPS = "Mortierella,Rhodotorula,Exophiala,Citrus"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--hits", required=True, help="hits TSV from validate_windows.py")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--records", help="amfsplit records TSV (adds no_window calls)")
    ap.add_argument("--fasta-prefix", help="amfsplit output prefix; writes <window>.amf.fasta subsets")
    ap.add_argument("--outgroups", default=DELAVAUX_OUTGROUPS,
                    help=f"comma-separated non-AMF genera (default: {DELAVAUX_OUTGROUPS})")
    ap.add_argument("--amf-pattern", help="regex marking AMF targets; overrides --outgroups")
    ap.add_argument("--min-id", type=float, default=90.0,
                    help="identity (%%) for a confident AMF call (default 90; ~85 for ONT R9.4.1)")
    ap.add_argument("--keep-unresolved", action="store_true",
                    help="include unresolved reads in the AMF subsets")
    args = ap.parse_args()

    outgroups = {g.strip() for g in args.outgroups.split(",") if g.strip()}
    amf_re = re.compile(args.amf_pattern) if args.amf_pattern else None

    best = {}
    for line in open(args.hits):
        q, t, ident = line.rstrip("\n").split("\t")[:3]
        ident = float(ident)
        if q not in best or ident > best[q][1]:
            best[q] = (t, ident)

    def is_amf(target):
        if amf_re:
            return bool(amf_re.search(target))
        parts = target.split("_")
        return len(parts) < 2 or parts[1] not in outgroups

    calls = {}
    for q, (t, ident) in best.items():
        if not is_amf(t):
            calls[q] = ("non_AMF", t, ident)
        elif ident >= args.min_id:
            calls[q] = ("AMF", t, ident)
        else:
            calls[q] = ("unresolved", t, ident)

    universe = list(best)
    if args.records:
        with open(args.records) as fh:
            universe = [row["read_id"] for row in csv.DictReader(fh, delimiter="\t")]
    # reads that had a window but no hit are in the window FASTA, not in the hits file
    windowed = set(best)
    if args.fasta_prefix:
        for path in glob.glob(f"{args.fasta_prefix}.LSU_LROR-FLR2.fasta"):
            windowed |= {l[1:].split()[0] for l in open(path) if l.startswith(">")}
    for rid in universe:
        if rid not in calls:
            calls[rid] = ("no_hit" if rid in windowed else "no_window", "", "")

    keep = {"AMF", "unresolved"} if args.keep_unresolved else {"AMF"}
    with open(f"{args.out}.triage.tsv", "w") as fh:
        fh.write("read_id\tcall\tbest_target\tidentity\n")
        for rid in universe if args.records else sorted(calls):
            c, t, i = calls[rid]
            fh.write(f"{rid}\t{c}\t{t}\t{i}\n")
    amf_ids = {rid for rid, (c, _, _) in calls.items() if c in keep}
    with open(f"{args.out}.amf.ids", "w") as fh:
        fh.write("".join(f"{rid}\n" for rid in sorted(amf_ids)))

    counts = Counter(c for c, _, _ in calls.values())
    n = sum(counts.values())
    non_amf_targets = Counter(t.split("_")[1] if "_" in t else t
                              for c, t, _ in calls.values() if c == "non_AMF")
    summary = {"reads": n,
               "calls": {k: {"n": v, "pct": round(100 * v / max(1, n), 1)} for k, v in counts.most_common()},
               "non_AMF_best_hits": non_amf_targets.most_common(6),
               "params": {"min_id": args.min_id, "amf_pattern": args.amf_pattern,
                          "outgroups": sorted(outgroups) if not amf_re else None,
                          "keep_unresolved": args.keep_unresolved}}

    if args.fasta_prefix:
        written = {}
        for path in sorted(glob.glob(f"{args.fasta_prefix}.*.fasta")):
            if path.endswith(".amf.fasta"):
                continue
            out_path = path[:-len(".fasta")] + ".amf.fasta"
            n_out, keep_rec = 0, False
            with open(path) as src, open(out_path, "w") as dst:
                for line in src:
                    if line.startswith(">"):
                        keep_rec = line[1:].split()[0] in amf_ids
                        n_out += keep_rec
                    if keep_rec:
                        dst.write(line)
            written[out_path.split("/")[-1]] = n_out
        summary["amf_fastas"] = written

    with open(f"{args.out}.triage.summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps({k: summary[k] for k in ("reads", "calls", "non_AMF_best_hits")}, indent=1))


if __name__ == "__main__":
    main()
