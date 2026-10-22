#!/usr/bin/env python3
"""Development estimates for protocol v3 (draft), from stored outputs of amf-novelty-benchmark only.

Every number is a development estimate: these data were seen before the rules were written.
Reads per-sequence score tables, locked predictions and follow-up tables; no sequences, so it runs
anywhere the repository is checked out (LongGloDB tables hold identifiers and derived numbers only).
Before each estimate the frozen count it starts from is reproduced; any mismatch exits 1.

Rules: R1a unnamed cultures (inside only an unnamed culture -> U1); R2 inside gate (a known call
stands only if the sequence lies inside a called species on all three markers); R3b placement flag
(pendant length of the best LSU placement above the 99th percentile of main-type culture copies in
the LSU-removal rounds); R5 motif step (YTYGGRGARTATG carrier whose matched VT links only to
Rhizoglomus, or, as a fallback, whose VT is not matched and whose nearest SSU culture is Rhizoglomus
within the within-genus SSU limit); R6 triage (reads below --min-id against V18, by edit-distance
identity of LR0R-FLR2 to the nearest V18 record, which approximates triage's vsearch identity);
R7 genus-clade fallback (a new-lineage call takes the genus of its LSU genus clade when two or more
training cultures of that one genus support it: at unit level where per-round placements of the
training cultures are stored, otherwise from the clade's label); R8 single-marker novelty (a new
lineage resting on one unmatched arm, the other two no call; the demoting variant calls it
insufficient evidence). Genus names are normalised as the frozen method does (normalise_genus).

Usage: python3 v3_dev_estimates.py --repo /path/to/amf-novelty-benchmark [--ssu-genus-limit 2.108]
       [--out v3_dev_estimates.tsv]
"""
import argparse, csv, gzip, os, statistics, sys
from collections import Counter

KNOWN = ("known species", "known species group")
ABSENT_SP, ABSENT_LIN = "novel species, genus in v4", "novel lineage, genus not in v4"
GIVES_GENUS = ("novel species, genus known", "known species", "known species group", "divergent class")
RHIZO = ("Rhizophagus", "Rhizoglomus", "[Rhizoglomus]")


def lab_of(roles, k):
    o = roles[k]["organism"]
    return o.split()[0] if o else ""


def norm(name):
    """Species key across the repository's spellings: Rhizophagus = Rhizoglomus, epithet stem."""
    p = name.replace("[", "").replace("]", "").split()
    if len(p) < 2:
        return name.strip()
    g = "Rhizoglomus" if p[0] in ("Rhizophagus", "Rhizoglomus") else p[0]
    return g + " " + p[1][:6].lower()


def names(field):
    return {norm(x) for x in field.replace("; ", ";").split(";") if x.strip()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--repo", required=True)
    ap.add_argument("--ssu-genus-limit", type=float, default=2.108, help="R5 fallback cut (v3_limits.py, SSU p99)")
    ap.add_argument("--out", default="v3_dev_estimates.tsv")
    a = ap.parse_args()
    J = lambda *p: os.path.join(a.repo, *p)

    def rd(*p):
        path = J(*p)
        op = gzip.open if path.endswith(".gz") else open
        return list(csv.DictReader(op(path, "rt"), delimiter="\t"))

    out, bad = [], []

    def row(rule, data, outcome, frozen, with_rule, note=""):
        out.append(dict(rule=rule, data=data, outcome=outcome, frozen=frozen, with_rule=with_rule, note=note))
        print(f"{rule:4s} | {data:44s} | {outcome:34s} | {frozen!s:>9s} -> {with_rule!s:<12s} {note}")

    def check(label, got, want):
        if got != want:
            bad.append(label)
        print(f"     check {label}: {got} (frozen {want}) {'ok' if got == want else 'MISMATCH'}")

    # ---------- EUKARYOME v2.1 ----------
    H = ("benchmarks", "m4_external_holdout")
    roles = {r["asv"]: r for r in rd("benchmarks", "eukaryome", "overlap.asvs.tsv")}
    pred = {r["asv"]: r for r in rd("benchmarks", "eukaryome", "predictions.predictions.tsv")}
    eA = {r["accession"]: r for r in rd(*H, "runs_euk", "06_scores", "A_primary_v2.copies.tsv.gz")}
    eB = {r["accession"]: r for r in rd(*H, "runs_euk", "06_scores", "B_primary_v2.copies.tsv.gz")}
    em = {k: eB.get(k) or eA[k] for k in eA}  # main: species/genus removed from V18 where V18 holds it
    lab = lambda k: roles[k]["organism"].split()[0] if roles[k]["organism"] else ""
    inside = lambda tab, k, p: bool(names(tab[k]["call_species"]) & names(p[k]["inside"]))
    fk = [k for k in em if roles[k]["role"] in (ABSENT_SP, ABSENT_LIN) and em[k]["status"] in KNOWN]
    check("EUKARYOME false known", len(fk), 13)
    g = [k for k in fk if inside(em, k, pred)]
    row("R2", "EUKARYOME, absent taxa", "false known", f"{len(fk)}/2706", f"{len(g)}/2706", ";".join(g))
    kn = [k for k in eA if roles[k]["role"] == "known species"]
    kept = [k for k in kn if eA[k]["status"] in KNOWN and eA[k]["correct"] in ("yes", "contains")]
    check("EUKARYOME kept", len(kept), 152)
    row("R2", "EUKARYOME, known species", "kept", f"{len(kept)}/{len(kn)}",
        f"{sum(inside(eA, k, pred) for k in kept)}/{len(kn)}")
    go = [k for k in eA if roles[k]["role"] == "genus only" and eA[k]["status"] in KNOWN]
    row("R2", "EUKARYOME, genus-only reads", "known calls outside called species", f"{len(go)} calls",
        f"{sum(not inside(eA, k, pred) for k in go)} demoted", "species unknown; descriptive")
    fkd = {r["accession"]: r for r in rd("paper", "analysis", "euk_followup", "false_known_distances.tsv")}
    u1 = [k for k in fk if not fkd[k]["inside_named"] and fkd[k]["inside_genus_only"]]
    row("R1a", "EUKARYOME, absent taxa", "false known", len(fk), len(fk) - len(u1),
        f"{len(u1)} inside only an unnamed culture -> U1; with R2: {len([k for k in g if k not in u1])}")
    gi = [r["asv"] for r in rd("benchmarks", "eukaryome", "genus_inside.reads.tsv") if r["inside_genus_only_culture"]]
    st = Counter(em[k]["status"] for k in gi if not pred[k]["inside"])
    row("R1a", "EUKARYOME, absent lineages inside unnamed cultures", "calls", len(gi),
        f"{sum(v for s, v in st.items() if s not in ('divergent class', 'artifact', 'discordant'))} -> U1",
        "; ".join(f"{s} {v}" for s, v in st.most_common()))
    # R5
    motif = [r["accession"] for r in rd(*H, "runs_euk", "report", "motif_reads.tsv")]
    def vt_rhizo(r):
        gs = {x for x in r["VT_genera"].split(";") if x}
        return r["VT_state"] == "matched" and bool(gs) and gs <= set(RHIZO)
    def fallback(r):
        nc = r["nearest_culture"].split(" | ")[0].split()[0] if r["nearest_culture"] else ""
        return (r["VT_state"] != "matched" and nc in RHIZO and r["nearest_ssu_pct"] != ""
                and float(r["nearest_ssu_pct"]) <= a.ssu_genus_limit)
    rz = [k for k in motif if lab(k) in RHIZO]
    other = [k for k in motif if k not in rz]
    now = sum(em[k]["status"] == "divergent class" for k in rz)
    check("EUKARYOME Rhizoglomus motif carriers", len(rz), 95)
    for name, f in (("VT only", vt_rhizo), (f"VT or SSU fallback {a.ssu_genus_limit}%", lambda r: vt_rhizo(r) or fallback(r))):
        rec = sum(em[k]["status"] == "divergent class" or f(em[k]) for k in rz)
        wrong = Counter(lab(k) for k in other if f(em[k]))
        row("R5", f"EUKARYOME motif carriers ({name})", "Rhizoglomus called divergent; others reclassified",
            f"{now}/{len(rz)}", f"{rec}/{len(rz)}; {sum(wrong.values())}/{len(other)}",
            "; ".join(f"{s} {v}" for s, v in wrong.items()))
    # R6
    nl = [k for k in em if roles[k]["role"] == ABSENT_LIN and lab(k) != "Nanoglomus" and em[k]["status"] == "novel lineage"]
    for t in (85, 88, 90, 92):
        n = sum(1 for k in nl if roles[k]["v18_best_pct"] and float(roles[k]["v18_best_pct"]) < t)
        ctrl = sum(1 for k in roles if roles[k]["role"] in ("known species", ABSENT_SP)
                   and roles[k]["v18_best_pct"] and float(roles[k]["v18_best_pct"]) < t)
        row("R6", f"EUKARYOME new-lineage calls, --min-id {t}", "dropped by filtering triage", f"0/{len(nl)}",
            f"{n}/{len(nl)}", f"known and absent-species reads below the cut: {ctrl}")
    # R3b
    sv2 = rd(*H, "aggregate_v2", "species.v2.copies.tsv")
    pend = [float(r["lsu_best_pendant_length"]) for r in sv2 if r["copy_class"] == "main" and r["lsu_best_pendant_length"]]
    thr = statistics.quantiles(pend, n=100)[98]
    pl = lambda r: float(r["lsu_best_pendant_length"]) if r["lsu_best_pendant_length"] else 0.0
    wrong_g = [k for k in em if roles[k]["role"] == ABSENT_LIN and lab(k) != "Nanoglomus" and em[k]["status"] in GIVES_GENUS]
    check("EUKARYOME absent genera given a genus", len(wrong_g), 268)
    row("R3b", "EUKARYOME, absent genera given a genus", "flagged", "no flag",
        f"{sum(pl(em[k]) > thr for k in wrong_g)}/{len(wrong_g)}", f"cut {thr:.5f} from {len(pend)} copies")
    # absent species with a genus outcome, Dentiscutata erythropus not counted (genus uncertain), as in the paper
    abs_sp = [k for k in em if roles[k]["role"] == ABSENT_SP and em[k]["status"] in GIVES_GENUS
              and not roles[k]["organism"].startswith(("Dentiscutata erythropus", "Fuscutata erythropus"))]
    right = [k for k in abs_sp if em[k]["final_genus_compatible"] in ("yes", "True", "1")]
    check("EUKARYOME absent species given the right genus", len(right), 49)
    row("R3b", "EUKARYOME, absent species, right genus", "false flags", "no flag",
        f"{sum(pl(em[k]) > thr for k in right)}/{len(right)}")
    row("R3b", "EUKARYOME, known species kept", "false flags", "no flag", f"{sum(pl(eA[k]) > thr for k in kept)}/{len(kept)}")

    # ---------- LongGloDB v1 (identifiers and derived numbers only) ----------
    lp = {r["asv"]: r for r in rd("benchmarks", "longglodb", "predictions.predictions.tsv")}
    lA = {r["accession"]: r for r in rd(*H, "runs_lgdb", "06_scores", "A_primary_v2.copies.tsv")}
    lB = {r["accession"]: r for r in rd(*H, "runs_lgdb", "06_scores", "B_primary_v2.copies.tsv")}
    lm = {k: lB.get(k) or lA[k] for k in lA}
    lfk = [k for k in lm if lm[k]["test_role"] in (ABSENT_SP, ABSENT_LIN) and lm[k]["status"] in KNOWN]
    check("LongGloDB false known", len(lfk), 4)
    row("R2", "LongGloDB, absent taxa", "false known", f"{len(lfk)}/102", f"{sum(inside(lm, k, lp) for k in lfk)}/102")
    lkn = [k for k in lA if lA[k]["test_role"] == "known species"]
    lkept = [k for k in lkn if lA[k]["status"] in KNOWN and lA[k]["correct"] in ("yes", "contains")]
    check("LongGloDB kept", len(lkept), 136)
    lost = [f"{lA[k]['organism']} {lA[k]['test_code']}" for k in lkept if not inside(lA, k, lp)]
    row("R2", "LongGloDB, known species", "kept", f"{len(lkept)}/{len(lkn)}", f"{len(lkept) - len(lost)}/{len(lkn)}", "; ".join(lost))

    # ---------- culture reference v4 ----------
    fw = {r["accession"]: r for r in rd("paper", "analysis", "foldwise.copies.tsv")}
    gate = lambda r: bool(names(r["call_species"]) & names(fw[r["accession"]]["inside_prospective"]))
    for tag, fn, want in (("V18 removed", "species.v2.copies.tsv", 53), ("V18 intact", "species.control_v2.copies.tsv", 9)):
        rows = [r for r in rd(*H, "aggregate_v2", fn) if r["copy_class"] == "main"
                and r["development_fold"] in ("", "0", "False", "no") and not r["organism"].startswith("Ambispora callosa")]
        f = [r for r in rows if r["status"] in KNOWN]
        check(f"LSU removal false known, {tag}", len(f), want)
        row("R2", f"v4, LSU removal, {tag}", "false known", f"{len(f)}/{len(rows)}", f"{sum(map(gate, f))}/{len(rows)}")
    sp = [r for r in rd(*H, "frozen_m2b", "comparison", "epa_species_clades_physical_species.copies.tsv")
          if r["copy_class"] == "main" and r["accession"] in fw]
    f = [r for r in sp if r["status"] in KNOWN]
    check("species withheld, placement arm, false known", len(f), 71)
    row("R2", "v4, species withheld, V18 intact (no N6)", "false known", f"{len(f)}/{len(sp)}",
        f"{sum(map(gate, f))}/{len(sp)}", "N6 alone leaves 18/803")

    # ---------- R7 genus-clade fallback, with R3b flags on what it adds ----------
    def gnorm(g):
        g = (g or "").strip().strip("[]")
        return {"Rhizoglomus": "Rhizophagus", "Claroideoglomus": "Entrophospora"}.get(g, g)
    cul = rd("reference", "culture_ref_v4", "cultures.tsv")
    gcount = Counter(gnorm(c["genus"]) for c in cul)
    cgenus = {c["culture"]: gnorm(c["genus"]) for c in cul}
    corg = {c["culture"]: c["organism"] for c in cul}
    cop = {r["accession"]: r for r in rd("reference", "culture_ref_v4", "copies.tsv") if r["status"] == "retained"}
    plc = {r["accession"]: r for r in rd(*H, "frozen_m2b", "placement", "summary", "lsu_placements.copies.tsv")}
    holders = {}
    for acc, r in plc.items():
        if r["external_genus_clade_unit"]:
            holders.setdefault(r["external_genus_clade_unit"], set()).add(cop[acc]["culture"])
    flagged = lambda r, key="lsu_best_pendant_length": float(r[key] or 0) > thr
    # (a) culture reference, species withheld, V18 intact: unit level
    sp = rd(*H, "frozen_m2b", "comparison", "epa_species_clades_physical_species.copies.tsv")
    res, fl = Counter(), Counter()
    for r in sp:
        if r["status"] != "novel lineage" or r["copy_class"] != "main":
            continue
        k = r["accession"]; u = plc[k]["external_genus_clade_unit"]; lab = gnorm(plc[k]["external_genus_clade"])
        org = cop[k]["organism"]; named = " sp." not in org and not org.endswith(" sp")
        held = {c for c in holders.get(u, ()) if (corg[c] != org if named else c != cop[k]["culture"])}
        gens = {cgenus[c] for c in held}
        v = ("right" if lab == gnorm(cop[k]["genus"]) else "wrong") if len(held) >= 2 and gens == {lab} else "abstain"
        res[v] += 1
        if v != "abstain" and flagged(plc[k], "best_pendant_length"):
            fl[v] += 1
    row("R7", "v4, species withheld, V18 intact (unit level)", "genus for F's new-lineage calls",
        f"0/{sum(res.values())}", f"right {res['right']}, wrong {res['wrong']}, abstain {res['abstain']}",
        f"R3b flags {fl['right']} right, {fl['wrong']} wrong")
    # (b) LSU removal, V18 removed: label level (round placements of training cultures are not stored)
    t = [r for r in sv2 if r["copy_class"] == "main" and r["development_fold"] in ("", "0", "False", "no")
         and not r["organism"].startswith("Ambispora callosa")]
    spc = Counter(c["organism"] for c in cul)
    res, fl = Counter(), Counter()
    for r in t:
        if r["status"] != "novel lineage":
            continue
        lab = gnorm(r["lsu_external_genus_clade"]); true = gnorm(r["genus"])
        support = gcount[lab] - (spc[r["organism"]] if lab == true else 0)
        v = ("right" if lab == true else "wrong") if lab and support >= 2 else "abstain"
        res[v] += 1
        if v != "abstain" and flagged(r):
            fl[v] += 1
    kept_genus = sum(1 for r in t if r["status"] != "novel lineage" and r["call_genera"])
    check("LSU removal, genus abstentions (genera without another culture)", res["abstain"], 186)
    check("LSU removal, genus another culture holds", kept_genus + res["right"] + res["wrong"], 487)
    row("R7", "v4, LSU removal, V18 removed (label level)", "genus for F's new-lineage calls",
        f"{kept_genus}/487 genus kept", f"right {res['right']}, wrong {res['wrong']}, abstain {res['abstain']}",
        f"R3b flags {fl['right']} right (in-sample cut), {fl['wrong']} wrong")
    # (c) EUKARYOME main and (d) LongGloDB main: label level
    for name, tab, role in (("EUKARYOME", em, lambda k: roles[k]["role"]), ("LongGloDB", lm, lambda k: lm[k]["test_role"])):
        res, fl = Counter(), Counter()
        for k, r in tab.items():
            if r["status"] != "novel lineage":
                continue
            ro = role(k)
            if name == "EUKARYOME" and lab_of(roles, k) == "Nanoglomus":
                continue
            labg = gnorm(r["lsu_external_genus_clade"])
            if not (labg and gcount[labg] >= 2):
                res[(ro, "abstain")] += 1
                continue
            org = roles[k]["organism"] if name == "EUKARYOME" else r["organism"]
            true = gnorm(org.split()[0]) if org else ""
            v = "wrong" if ro == ABSENT_LIN else ("right" if labg == true else "wrong") if ro in (ABSENT_SP, "known species") else "label only"
            res[(ro, v)] += 1
            if flagged(r):
                fl[(ro, v)] += 1
        for ro in (ABSENT_SP, ABSENT_LIN, "known species"):
            got = {v: res[(ro, v)] for v in ("right", "wrong", "abstain") if res[(ro, v)]}
            if got:
                row("R7", f"{name}, {ro}", "genus for F's new-lineage calls", "no genus",
                    ", ".join(f"{v} {n}" for v, n in got.items()),
                    f"R3b flags {fl[(ro, 'right')]} right, {fl[(ro, 'wrong')]} wrong")

    # ---------- R8 single-marker novelty, on the main arm ----------
    ARMS = ("VT_state", "SH_state", "LSU_state")
    def r8(rows, data, stratum):
        nl = [r for r in rows if r["status"] == "novel lineage"]
        one = [r for r in nl if sum(r[x] == "unmatched" for x in ARMS) < 2]
        row("R8", data, "new lineage on one unmatched arm", f"{len(nl)} new lineage", f"{len(one)} flagged",
            "; ".join(f"{s} {n}" for s, n in Counter(map(stratum, one)).most_common()))
    for mode in ("physical_culture", "physical_species"):
        r8(rd(*H, "frozen_m2b", "comparison", f"epa_species_clades_{mode}.copies.tsv"),
           f"v4, species clades, {mode}", lambda r: r["copy_class"])
    r8([r for r in sv2 if r["development_fold"] in ("", "0", "False", "no") and not r["organism"].startswith("Ambispora callosa")],
       "v4, LSU removal, V18 removed", lambda r: r["copy_class"])
    r8(list(em.values()), "EUKARYOME main", lambda r: roles[r["accession"]]["role"])
    r8(list(lm.values()), "LongGloDB main", lambda r: r["test_role"])

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(out)
    if bad:
        print("MISMATCH:", ", ".join(bad)); return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
