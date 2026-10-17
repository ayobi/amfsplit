#!/usr/bin/env python3
"""amfsplit v0.3.1 (Stage 1 prototype).

Turns each long AMF rDNA read (SSU-ITS-LSU amplicon, PacBio HiFi or ONT) into a
multi-marker record:

  * ITS regions (ITS1, 5.8S, ITS2, full ITS) and the SSU / LSU flanks, taken
    from ITSxRust's per-read anchor table (--anchors-tsv), and
  * "legacy barcode windows": the exact short-read regions that the MaarjAM,
    SSU-V4, LSU-D2 and Delavaux LSU literature was built on, cut from the long
    read at the same primer sites (IUPAC-aware, mismatch-tolerant search).

Cutting windows at the original primer sites keeps every window directly
comparable with its legacy reference database, which is what the Stage 3
VT <-> SH <-> LSU-clade bridge needs.

Coordinates in the output are 1-based, inclusive, in the rRNA-sense
orientation of the read (the same convention ITSxRust uses for minus-strand
reads).

v0.3 changes (the fixes queued after the Stefani and ONT checks):
  * Region-limited primer search. Each primer is searched only inside the part of
    the read where it can lie: SSU primers before the ITSxRust full-ITS start, LSU
    primers after its end, plus --region-margin (default 50 nt). Reads without
    anchors fall back to the whole read (region_limited = no). Panel columns
    fwd_region / rev_region (ssu, lsu, any) override the defaults, which follow the
    window name (SSU_*: ssu; LSU_*: lsu; AMPLICON_*: ssu then lsu).
  * missing_5p / missing_3p. A primer that is absent although the read covers the
    place where it should lie is no longer relabelled open or too_long. A site
    counts as beyond the read (open_5p / open_3p) only where the search region is
    cut by the read end and the window would still be shorter than the panel's
    max_len. Previously AM1 failures on full-length copies came out as too_long.
  * Explicit concatemer screen, over the whole read and independent of the region
    limit: a read is flagged when at least --concat-min-primers different primers
    each occur twice, at least --concat-min-sep apart, with the same spacing (a
    tandem copy), or each occur on the opposite strand with the same mirror sum
    (inverted repeat). Agreement of the spacing across primers is what separates a
    concatemer from chance primer hits in a noisy read (v0.3.0 flagged any repeated
    primers and called single-length ONT reads concatemers). Flagged reads stay in
    records.tsv (concatemer = yes, every window status = concatemer) and are left
    out of the FASTA outputs unless --keep-concatemers is set. concat_kind separates
    a full repeat (copy spacing at least 0.4 of the read length, or the mirror centre in
    the middle third: two comparable copies) from a partial one (a short segment
    repeated, typically an end duplicated or folded back, which makes the read only
    10-15% longer).
  * records.tsv uses LF line endings (the csv default is CRLF).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import statistics
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass

import edlib

# ----------------------------------------------------------------------------- sequence helpers
IUPAC = {"R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
         "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
EQUALITIES = [(code, base) for code, bases in IUPAC.items() for base in bases]
_COMP = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def open_text(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_sequences(path: str):
    """Yield (read_id, sequence) from FASTA or FASTQ, gzipped or not."""
    with open_text(path) as fh:
        first = fh.readline()
        if not first:
            return
        if first.startswith("@"):
            header = first
            while header:
                seq = fh.readline().strip()
                fh.readline()
                fh.readline()
                yield header[1:].split()[0], seq.upper()
                header = fh.readline()
        elif first.startswith(">"):
            rid, chunks = first[1:].split()[0], []
            for line in fh:
                if line.startswith(">"):
                    yield rid, "".join(chunks).upper()
                    rid, chunks = line[1:].split()[0], []
                else:
                    chunks.append(line.strip())
            yield rid, "".join(chunks).upper()
        else:
            sys.exit(f"error: {path} is neither FASTA nor FASTQ")


# ----------------------------------------------------------------------------- panel
@dataclass
class Window:
    name: str
    target: str
    fwd_name: str
    fwd: str
    rev_name: str
    rev: str
    min_len: int
    max_len: int
    fwd_region: str = "any"
    rev_region: str = "any"

    @property
    def rev_site(self) -> str:
        """Reverse primer as it appears in the rRNA-sense read."""
        return revcomp(self.rev)


REGIONS = ("ssu", "lsu", "any")


def default_regions(name: str) -> tuple[str, str]:
    """Where each primer of a window lies relative to the ITS, from the window name."""
    n = name.upper()
    if n.startswith("SSU"):
        return "ssu", "ssu"
    if n.startswith("LSU"):
        return "lsu", "lsu"
    if n.startswith("AMPLICON"):
        return "ssu", "lsu"
    return "any", "any"


def load_panel(path: str) -> list[Window]:
    windows = []
    with open(path) as fh:
        rows = (line for line in fh if line.strip() and not line.startswith("#"))
        for row in csv.DictReader(rows, delimiter="\t"):
            fr, rr = default_regions(row["window"])
            fr = (row.get("fwd_region") or fr).lower()
            rr = (row.get("rev_region") or rr).lower()
            if fr not in REGIONS or rr not in REGIONS:
                sys.exit(f"error: panel window {row['window']}: region must be one of {REGIONS}")
            windows.append(Window(row["window"], row["target"], row["fwd_name"],
                                  row["fwd_seq"].upper(), row["rev_name"],
                                  row["rev_seq"].upper(), int(row["min_len"]),
                                  int(row["max_len"]), fr, rr))
    return windows


def find_site(site: str, seq: str, max_err: int):
    """Best infix match of a primer site in seq -> (start0, end0, edits, n_best) or None."""
    r = edlib.align(site, seq, mode="HW", task="locations", k=max_err,
                    additionalEqualities=EQUALITIES)
    if r["editDistance"] < 0:
        return None
    start, end = r["locations"][0]
    return start, end, r["editDistance"], len(r["locations"])


def max_edits(primer: str, frac: float) -> int:
    return max(1, int(len(primer) * frac))


def find_site_in(site: str, seq: str, max_err: int, lo: int, hi: int):
    """find_site restricted to seq[lo:hi]; coordinates stay 0-based on the whole read."""
    if hi - lo < len(site) - max_err:
        return None
    hit = find_site(site, seq[lo:hi], max_err)
    if hit is None:
        return None
    start, end, edits, n_best = hit
    return start + lo, end + lo, edits, n_best


def region_bounds(region: str, full, read_len: int, margin: int):
    """0-based half-open search bounds of a primer's region.

    ssu: before the ITSxRust full-ITS start; lsu: after its end; each widened by `margin`.
    Without a full-ITS anchor (or region 'any') the whole read is searched.
    """
    if full is None or region == "any":
        return 0, read_len
    if region == "ssu":
        return 0, min(read_len, full[0] - 1 + margin)
    return max(0, full[1] - margin), read_len


def find_all_sites(site: str, seq: str, max_err: int, min_sep: int, max_sites: int = 3):
    """Sites of one primer at least min_sep apart: take the best match, mask it and its
    surroundings with 'X' (which no primer base or IUPAC code equals), search again."""
    masked, sites = seq, []
    for _ in range(max_sites):
        hit = find_site(site, masked, max_err)
        if hit is None:
            break
        start, end, edits, _ = hit
        sites.append((start, end, edits))
        lo, hi = max(0, start - min_sep), min(len(masked), end + min_sep + 1)
        masked = masked[:lo] + "X" * (hi - lo) + masked[hi:]
    return sites


def _best_cluster(named_values: dict, tol_frac: float, tol_min: int):
    """Largest set of primers whose values agree within max(tol_min, tol_frac * value)."""
    best = []
    for v0 in named_values.values():
        tol = max(tol_min, tol_frac * abs(v0))
        members = [n for n, v in named_values.items() if abs(v - v0) <= tol]
        if len(members) > len(best):
            best = members
    return sorted(best)


def _norm_ed(a: str, b: str) -> float:
    """Global edit distance per base of the longer sequence (0 = identical, ~0.5+ for unrelated DNA)."""
    if not a or not b:
        return 1.0
    return edlib.align(a, b, mode="NW", task="distance")["editDistance"] / max(len(a), len(b))


def concatemer_screen(seq: str, primers: dict, frac: float, min_sep: int, min_primers: int,
                      tol_frac: float = 0.05, max_dist: float = 0.40, min_unit: int = 100,
                      strong_primers: int = 3):
    """(is_concatemer, evidence, kind). primers: {name: rRNA-sense site sequence}.

    Two copies in one read repeat every primer site with the same spacing:
      tandem    each primer occurs at two sites >= min_sep apart, and at least min_primers of
                them share one offset (within max(40 nt, tol_frac of it): the repeat unit);
      inverted  each primer also occurs on the opposite strand, and at least min_primers of them
                share one mirror sum (sense position + antisense position + site length).
    Chance primer hits in a noisy read can agree on a spacing by accident, so a cluster of only
    min_primers primers must also pass a primer-free test: the two supposed copies (up to 800 nt)
    must resemble each other, at most max_dist edits per base (unrelated DNA is about 0.5; two
    ONT reads of the same sequence are about 0.3). A cluster of strong_primers or more primers
    agreeing on one spacing is accepted without it: agreement of three primers by chance is rare,
    and the two units of a real concatemer are often different molecules (other taxa or rDNA
    variants) whose ITS does not resemble each other. kind is 'full' for two comparable units
    and 'partial' when only a segment is repeated."""
    n = len(seq)
    rc = revcomp(seq)
    offsets, mirrors = {}, {}
    for name, site in primers.items():
        k = max_edits(site, frac)
        sites = sorted(x[0] for x in find_all_sites(site, seq, k, min_sep, max_sites=3))
        if len(sites) >= 2:
            offsets[name] = (sites[1] - sites[0], sites[0])
        anti = find_site(site, rc, k)
        if anti and sites:
            mirrors[name] = (sites[0] + (n - 1 - anti[1]) + len(site), sites[0])
    evidence, kinds = [], []
    if offsets:
        cluster = _best_cluster({k: v[0] for k, v in offsets.items()}, tol_frac, 40)
        if len(cluster) >= min_primers:
            period = int(statistics.median(offsets[t][0] for t in cluster))
            pos = int(statistics.median(offsets[t][1] for t in cluster))
            w = min(period, n - (pos + period), 800)
            strong = len(cluster) >= strong_primers
            d = _norm_ed(seq[pos:pos + w], seq[pos + period:pos + period + w]) if (w >= min_unit and not strong) else None
            if strong or (d is not None and d <= max_dist):
                evidence.append(f"tandem:{','.join(cluster)}@{period}" + ("" if d is None else f"~{d:.2f}"))
                kinds.append("full" if period >= 0.4 * n else "partial")
    if mirrors:
        cluster = _best_cluster({k: v[0] for k, v in mirrors.items()}, tol_frac, 40)
        if len(cluster) >= min_primers:
            centre = int(statistics.median(mirrors[t][0] for t in cluster) / 2)
            w = min(centre, n - centre, 800)
            strong = len(cluster) >= strong_primers
            d = _norm_ed(seq[centre - w:centre], revcomp(seq[centre:centre + w])) if (w >= min_unit and not strong) else None
            if strong or (d is not None and d <= max_dist):
                evidence.append(f"inverted:{','.join(cluster)}@{centre}" + ("" if d is None else f"~{d:.2f}"))
                kinds.append("full" if 0.35 * n <= centre <= 0.65 * n else "partial")
    kind = "" if not kinds else ("full" if "full" in kinds else "partial")
    return bool(evidence), ";".join(evidence), kind


# ----------------------------------------------------------------------------- ITSxRust anchors
def load_anchors(path: str | None) -> dict[str, dict]:
    if not path:
        return {}
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["read_id"]] = row
    return out


def _int(v: str | None):
    return int(v) if v not in (None, "") else None


def itsx_regions(row: dict, read_len: int) -> dict[str, tuple[int, int] | None]:
    """ITS sub-regions and flanks (1-based inclusive, oriented) from an anchors row."""
    reg: dict[str, tuple[int, int] | None] = {}
    i1 = (_int(row.get("its1_start")), _int(row.get("its1_end")))
    i2 = (_int(row.get("its2_start")), _int(row.get("its2_end")))
    fl = (_int(row.get("full_start")), _int(row.get("full_end")))
    reg["its1"] = i1 if None not in i1 else None
    reg["its2"] = i2 if None not in i2 else None
    reg["full_its"] = fl if None not in fl else None
    reg["s58"] = (i1[1] + 1, i2[0] - 1) if reg["its1"] and reg["its2"] else None
    reg["ssu_flank"] = (1, fl[0] - 1) if reg["full_its"] and fl[0] > 1 else None
    reg["lsu_flank"] = (fl[1] + 1, read_len) if reg["full_its"] and fl[1] < read_len else None
    return reg


def orientation_from_anchors(row: dict) -> str | None:
    for key in ("ssu_strand", "s58s_strand", "s58e_strand", "lsu_strand"):
        if row.get(key) in ("+", "-"):
            return row[key]
    return None


# ----------------------------------------------------------------------------- core
ITSX_REGIONS = ["ssu_flank", "its1", "s58", "its2", "full_its", "lsu_flank"]
USABLE = ("ok", "open_5p", "open_3p")


def orient_by_primers(seq: str, windows: list[Window], frac: float) -> str:
    """Fallback orientation when ITSxRust gave no anchors: majority of primer sites."""
    rc = revcomp(seq)
    score = {"+": 0, "-": 0}
    for w in windows:
        for site in (w.fwd, w.rev_site):
            k = max_edits(site, frac)
            if find_site(site, seq, k):
                score["+"] += 1
            if find_site(site, rc, k):
                score["-"] += 1
    if score["+"] == score["-"] == 0:
        return "?"
    return "+" if score["+"] >= score["-"] else "-"


def cut_window(w: Window, seq: str, frac: float, keep_primers: bool, open_ends: bool = False,
               fwd_bounds=None, rev_bounds=None):
    n = len(seq)
    fb = fwd_bounds or (0, n)
    rb = rev_bounds or (0, n)
    f = find_site_in(w.fwd, seq, max_edits(w.fwd, frac), *fb)
    r = find_site_in(w.rev_site, seq, max_edits(w.rev_site, frac), *rb)
    rec = {"fwd_ed": f[2] if f else "", "rev_ed": r[2] if r else "",
           "multi": int(bool(f and f[3] > 1) or bool(r and r[3] > 1))}
    if not f and not r:
        return {**rec, "status": "both_missing"}
    open_tag = ""
    if not f or not r:
        # One primer site not found. It is "open" (beyond the read) only if the read cuts its
        # search region and the window would still be shorter than the panel's max_len;
        # otherwise the read covers the place where the primer should be: missing.
        if not f:
            s0, e0 = 0, (r[1] if keep_primers else r[0] - 1)
            beyond_read = fb[0] == 0 and (e0 - s0 + 1) <= w.max_len
            missing, tag = "missing_5p", "open_5p"
        else:
            s0, e0 = (f[0] if keep_primers else f[1] + 1), n - 1
            beyond_read = rb[1] == n and (e0 - s0 + 1) <= w.max_len
            missing, tag = "missing_3p", "open_3p"
        if not (open_ends and beyond_read):
            return {**rec, "status": missing}
        open_tag = tag
    elif keep_primers:
        s0, e0 = f[0], r[1]
    else:
        s0, e0 = f[1] + 1, r[0] - 1
    if e0 < s0:
        return {**rec, "status": "out_of_order"}
    length = e0 - s0 + 1
    min_len = w.min_len // 2 if open_tag else w.min_len
    status = open_tag or "ok"
    if length < min_len:
        status = "too_short"
    elif length > w.max_len:
        status = "too_long"
    return {**rec, "status": status, "start": s0 + 1, "end": e0 + 1, "len": length}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--reads", required=True, help="reads (FASTA/FASTQ, .gz ok)")
    ap.add_argument("-a", "--anchors", help="ITSxRust --anchors-tsv for the same reads")
    ap.add_argument("-p", "--panel", required=True, help="primer-window panel TSV")
    ap.add_argument("-o", "--out", required=True, help="output prefix")
    ap.add_argument("--max-err-frac", type=float, default=0.15,
                    help="max edits per primer as a fraction of its length (default 0.15)")
    ap.add_argument("--keep-primers", action="store_true",
                    help="include primer sites in windows (default: primer-trimmed)")
    ap.add_argument("--open-ends", action="store_true",
                    help="when one primer site lies outside the amplicon, cut the window to the read "
                         "end and flag it open_5p/open_3p (e.g. NS31 windows on nested amplicons "
                         "that start downstream of NS31)")
    ap.add_argument("--only-anchored", action="store_true",
                    help="skip reads ITSxRust could not anchor (e.g. off-target host rDNA)")
    ap.add_argument("--region-margin", type=int, default=50,
                    help="search each primer inside its ITSxRust region widened by this many nt (default 50)")
    ap.add_argument("--no-region-limit", action="store_true",
                    help="search every primer over the whole read (the v0.2 behaviour)")
    ap.add_argument("--no-concat-screen", action="store_true", help="skip the concatemer screen")
    ap.add_argument("--concat-min-primers", type=int, default=2,
                    help="different primers that must show a repeat (or an inverted site) to flag a read (default 2)")
    ap.add_argument("--concat-min-sep", type=int, default=150,
                    help="min distance between two sites of one primer to count as two copies (default 150; "
                         "partial repeats of the read end can be short, and agreement of the spacing "
                         "across primers, not this distance, guards against chance hits)")
    ap.add_argument("--concat-sim", type=float, default=0.40,
                    help="max edits per base between the two supposed copies (default 0.40; unrelated DNA is ~0.5)")
    ap.add_argument("--concat-strong", type=int, default=3,
                    help="primers agreeing on one spacing that need no copy-similarity check (default 3)")
    ap.add_argument("--concat-tol", type=float, default=0.05,
                    help="primers must agree on the copy spacing within this fraction (default 0.05, min 40 nt)")
    ap.add_argument("--concat-err-frac", type=float,
                    help="max edits per primer for the concatemer screen (default: --max-err-frac)")
    ap.add_argument("--keep-concatemers", action="store_true",
                    help="write windows and regions of flagged concatemers too (they stay flagged in records.tsv)")
    args = ap.parse_args(argv)

    windows = load_panel(args.panel)
    anchors = load_anchors(args.anchors)

    region_names = ITSX_REGIONS + [w.name for w in windows]
    fasta = {name: open(f"{args.out}.{name}.fasta", "w") for name in region_names}
    rec_fh = open(f"{args.out}.records.tsv", "w", newline="")
    cols = ["read_id", "read_len", "orient", "orient_source", "itsx_confidence"]
    cols += [f"{r}_coords" for r in ITSX_REGIONS]
    for w in windows:
        cols += [f"{w.name}_{k}" for k in ("status", "start", "end", "len", "fwd_ed", "rev_ed", "multi")]
    cols += ["windows_ok", "flank_consistent", "region_limited", "concatemer", "concat_evidence", "concat_kind"]
    writer = csv.DictWriter(rec_fh, fieldnames=cols, delimiter="\t", extrasaction="ignore",
                            lineterminator="\n")
    writer.writeheader()

    by_site: dict = {}
    for w in windows:
        by_site.setdefault(w.fwd, w.fwd_name)
        by_site.setdefault(w.rev_site, w.rev_name)
    screen_primers = {name: site for site, name in by_site.items()}
    concat_frac = args.concat_err_frac if args.concat_err_frac is not None else args.max_err_frac

    n_reads = 0
    n_processed = 0
    n_concat = 0
    n_limited = 0
    concat_kinds = Counter()
    status_counts = {w.name: Counter() for w in windows}
    lengths = defaultdict(list)
    orient_src = Counter()
    region_counts = Counter()
    inconsistent = Counter()

    for rid, raw in read_sequences(args.reads):
        n_reads += 1
        row = anchors.get(rid)
        if row is None and args.only_anchored:
            orient_src["skipped_unanchored"] += 1
            continue
        n_processed += 1
        orient = orientation_from_anchors(row) if row else None
        source = "itsxrust" if orient else "primers"
        if orient is None:
            orient = orient_by_primers(raw, windows, args.max_err_frac)
        orient_src[source if orient != "?" else "unresolved"] += 1
        seq = revcomp(raw) if orient == "-" else raw

        rec = {"read_id": rid, "read_len": len(seq), "orient": orient,
               "orient_source": source if orient != "?" else "unresolved",
               "itsx_confidence": (row or {}).get("confidence", "unanchored")}

        regions = itsx_regions(row, len(seq)) if row else {r: None for r in ITSX_REGIONS}
        full = regions.get("full_its")
        limited = full is not None and not args.no_region_limit
        rec["region_limited"] = "yes" if limited else "no"
        n_limited += limited

        # explicit concatemer screen: whole read, independent of the region limit
        is_concat, evidence, kind = False, "", ""
        if not args.no_concat_screen and orient != "?":
            is_concat, evidence, kind = concatemer_screen(seq, screen_primers, concat_frac,
                                                          args.concat_min_sep, args.concat_min_primers,
                                                          args.concat_tol, args.concat_sim,
                                                          strong_primers=args.concat_strong)
        rec["concatemer"] = "yes" if is_concat else "no"
        rec["concat_evidence"] = evidence
        rec["concat_kind"] = kind
        n_concat += is_concat
        if is_concat:
            concat_kinds[kind] += 1
        skip_output = is_concat and not args.keep_concatemers

        for name in ITSX_REGIONS:
            coords = regions.get(name)
            rec[f"{name}_coords"] = f"{coords[0]}-{coords[1]}" if coords else ""
            if coords and coords[1] >= coords[0] and not skip_output:
                region_counts[name] += 1
                lengths[name].append(coords[1] - coords[0] + 1)
                fasta[name].write(f">{rid} region={name} coords={coords[0]}-{coords[1]} "
                                  f"orient={orient}\n{seq[coords[0] - 1:coords[1]]}\n")

        if skip_output:
            # flagged, not dropped: the read stays in records.tsv with every window marked
            for w in windows:
                status_counts[w.name]["concatemer"] += 1
                rec[f"{w.name}_status"] = "concatemer"
            rec["windows_ok"] = 0
            rec["flank_consistent"] = "n/a"
            writer.writerow(rec)
            continue

        n_ok = 0
        flank_ok = True
        for w in windows:
            if orient == "?":
                res = {"status": "unoriented"}
            else:
                fb = region_bounds(w.fwd_region, full, len(seq), args.region_margin) if limited else None
                rb = region_bounds(w.rev_region, full, len(seq), args.region_margin) if limited else None
                res = cut_window(w, seq, args.max_err_frac, args.keep_primers, args.open_ends, fb, rb)
            status_counts[w.name][res["status"]] += 1
            for k in ("status", "start", "end", "len", "fwd_ed", "rev_ed", "multi"):
                rec[f"{w.name}_{k}"] = res.get(k, "")
            if res["status"] in USABLE:
                n_ok += 1
                lengths[w.name].append(res["len"])
                fasta[w.name].write(f">{rid} window={w.name} coords={res['start']}-{res['end']} "
                                    f"status={res['status']} orient={orient} "
                                    f"fwd_ed={res['fwd_ed']} rev_ed={res['rev_ed']}\n"
                                    f"{seq[res['start'] - 1:res['end']]}\n")
                # cross-check primer-derived windows against HMM anchors
                if full:
                    if w.name.startswith("SSU") and res["end"] >= full[0]:
                        flank_ok = False
                        inconsistent[w.name] += 1
                    if w.name.startswith("LSU") and res["start"] <= full[1]:
                        flank_ok = False
                        inconsistent[w.name] += 1
        rec["windows_ok"] = n_ok
        rec["flank_consistent"] = "yes" if flank_ok else "no"
        writer.writerow(rec)

    for fh in fasta.values():
        fh.close()
    rec_fh.close()

    def dist(v):
        if not v:
            return None
        q = statistics.quantiles(v, n=4) if len(v) > 1 else [v[0]] * 3
        return {"n": len(v), "median": statistics.median(v), "q1": q[0], "q3": q[2],
                "min": min(v), "max": max(v)}

    summary = {
        "reads": n_reads,
        "reads_processed": n_processed,
        "orientation_source": dict(orient_src),
        "itsxrust_regions": {r: {"reads": region_counts[r], "length": dist(lengths[r])}
                             for r in ITSX_REGIONS},
        "windows": {w.name: {"target": w.target,
                             "primers": f"{w.fwd_name}/{w.rev_name}",
                             "ok": status_counts[w.name]["ok"],
                             "open": status_counts[w.name]["open_5p"] + status_counts[w.name]["open_3p"],
                             "usable": sum(status_counts[w.name][k] for k in USABLE),
                             "ok_pct": round(100 * status_counts[w.name]["ok"] / max(1, n_processed), 1),
                             "usable_pct": round(100 * sum(status_counts[w.name][k] for k in USABLE)
                                                 / max(1, n_processed), 1),
                             "status": dict(status_counts[w.name]),
                             "length": dist(lengths[w.name]),
                             "crosses_its_anchor": inconsistent[w.name]}
                    for w in windows},
        "concatemers": {"flagged": n_concat, "pct": round(100 * n_concat / max(1, n_processed), 2),
                        "full": concat_kinds["full"], "partial": concat_kinds["partial"],
                        "kept_in_fasta": args.keep_concatemers},
        "region_limited_reads": n_limited,
        "params": {"version": "0.3.1", "max_err_frac": args.max_err_frac, "keep_primers": args.keep_primers,
                   "open_ends": args.open_ends, "only_anchored": args.only_anchored,
                   "region_margin": args.region_margin, "region_limit": not args.no_region_limit,
                   "concat_screen": not args.no_concat_screen, "concat_min_primers": args.concat_min_primers,
                   "concat_min_sep": args.concat_min_sep, "concat_tol": args.concat_tol, "concat_sim": args.concat_sim, "concat_strong": args.concat_strong, "concat_err_frac": concat_frac},
    }
    with open(f"{args.out}.summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)

    # short console report
    print(f"{n_reads} reads ({n_processed} processed)  orientation: {dict(orient_src)}",
          file=sys.stderr)
    print(f"  concatemers flagged: {n_concat} ({summary['concatemers']['pct']}%; "
          f"{concat_kinds['full']} full, {concat_kinds['partial']} partial)  "
          f"region-limited reads: {n_limited}", file=sys.stderr)
    for r in ITSX_REGIONS:
        d = summary["itsxrust_regions"][r]
        med = d["length"]["median"] if d["length"] else "-"
        print(f"  {r:<24} {d['reads']:>6} reads  median len {med}", file=sys.stderr)
    for w in windows:
        d = summary["windows"][w.name]
        med = d["length"]["median"] if d["length"] else "-"
        print(f"  {w.name:<24} {d['usable']:>6} usable ({d['usable_pct']}%; {d['open']} open)  "
              f"median len {med}  "
              f"{dict(status_counts[w.name])}", file=sys.stderr)


if __name__ == "__main__":
    main()
