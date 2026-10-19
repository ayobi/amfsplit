#!/usr/bin/env python3
"""Phase 3: reconcile the VT, SH and LSU sets of every copy into one call, scored out of sample.

Each arm's compatible sets come from arm_sets.py (or vt_sets.py for VT). A unit (VT, SH group, LSU
unit) is linked to the species and genera of the cultures whose culture set holds it. To call a copy,
its own culture is left out of those links (--holdout culture, the default), so no copy is scored
against itself; --holdout species leaves out every culture of the copy's species instead, which tests
what happens to a species the reference does not know (the Phase 5 leave-one-lineage-out in small).

Per arm a copy is:
  matched    its set holds units linked to other cultures: candidate species and genera
  unmatched  its set holds only units no other culture has (unknown to the reference)
  no call    empty set (below the arm's identity threshold)
The call combines the matched arms: genera are intersected over all matched arms, species over the
matched arms that name species (an arm whose units belong only to 'sp.' cultures constrains genus only).
  species        one species left
  species group  several species left (e.g. species sharing VTX00199)
  genus          no species survives, one or more genera do
  discordant     matched arms share no genus
  unplaced       no arm matched
Known causes attached to a copy: divergent class (short 5.8S or same-length), and held out (no other
culture of its species in the reference, so a species call is impossible and a genus call is the best
outcome). In environmental data the divergent class has to be recognised from the read itself (Phase 4);
here the reference's copy classes stand in for that test.

Report:
  1. arm states under the holdout
  2. main-type copies, by group:
       named, species elsewhere  the M3 negative control: the true species should be called or kept
       named, species held out   no other culture of the species: a species call here is a false species
       unnamed ('sp.')           genus is the best possible call
  3. each arm alone vs the combined call, for 'named, species elsewhere'
  4. copies that fail the negative control, with their arms
  5. divergent-class copies
  6. novelty rules: an 'unmatched' arm says the copy is unlike anything else in the reference. The
     combined call ignores that; the rules below demote a species call to genus when arms are
     unmatched, and the table shows the trade-off: true species kept for 'species elsewhere' against
     false species for 'held out'. Run with --holdout species as well for the fuller novelty test.
Writes {out}.copies.tsv.

Usage (from amfsplit_stage3):
  python3 reconcile.py --ref ../ref/culture_ref/v4 --vt ../bench/stefani/vt_sets_v4 \\
      --sh ../bench/stefani/sets_SH_v4 --lsu ../bench/stefani/sets_LSU_v4 --out ../bench/stefani/reconcile_v4
"""
import argparse
import collections
import os
import sys

from arms import GENUS_SYNONYMS, is_named, read_tsv, write_tsv

ARMS = ("VT", "SH", "LSU")


def split_set(s):
    return frozenset(u for u in s.split(",") if u) if s else frozenset()


def load_arm(prefix):
    cop = read_tsv(f"{prefix}.copies.tsv")
    col = "unit_set" if "unit_set" in cop[0] else "vt_set"
    per_copy = {r["accession"]: split_set(r[col]) for r in cop}
    cul = {r["culture"]: split_set(r["culture_set"]) for r in read_tsv(f"{prefix}.cultures.tsv")}
    holders = collections.defaultdict(set)
    for c, units in cul.items():
        for u in units:
            holders[u].add(c)
    return per_copy, holders


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True, help="frozen culture reference directory")
    ap.add_argument("--vt", help="prefix of the VT sets (vt_sets.py or arm_sets.py --arm VT)")
    ap.add_argument("--sh", help="prefix of the SH sets (arm_sets.py --arm SH)")
    ap.add_argument("--lsu", help="prefix of the LSU sets (arm_sets.py --arm LSU)")
    ap.add_argument("--holdout", choices=["culture", "species"], default="culture")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = [r for r in read_tsv(os.path.join(a.ref, "copies.tsv")) if r["status"] == "retained"]
    cul_meta = {r["culture"]: r for r in read_tsv(os.path.join(a.ref, "cultures.tsv"))}
    arms = {name: load_arm(p) for name, p in (("VT", a.vt), ("SH", a.sh), ("LSU", a.lsu)) if p}
    if not arms:
        raise SystemExit("give at least one of --vt, --sh, --lsu")

    genus = {c: GENUS_SYNONYMS.get(m["genus"], m["genus"]) for c, m in cul_meta.items()}
    species = {c: m["organism"] for c, m in cul_meta.items()
               if is_named(m["organism"]) and m.get("name_status", "ok") == "ok"}
    by_species = collections.defaultdict(set)
    for c, sp in species.items():
        by_species[sp].add(c)

    def held_out(c):
        if a.holdout == "species" and c in species:
            return by_species[species[c]]
        return {c}

    out = []
    for r in rows:
        x, c = r["accession"], r["culture"]
        drop = held_out(c)
        per_arm = {}
        for name, (per_copy, holders) in arms.items():
            units = per_copy.get(x, frozenset())
            if not units:
                per_arm[name] = ("no call", set(), set())
                continue
            others = set().union(*(holders.get(u, set()) for u in units)) - drop
            sp = {species[o] for o in others if o in species}
            gn = {genus[o] for o in others if genus[o]}
            per_arm[name] = ("matched" if gn else "unmatched", sp, gn)
        matched = [n for n in per_arm if per_arm[n][0] == "matched"]
        if not matched:
            rank, call_sp, call_gn = "unplaced", set(), set()
        else:
            call_gn = set.intersection(*(per_arm[n][2] for n in matched))
            sp_arms = [n for n in matched if per_arm[n][1]]
            call_sp = set.intersection(*(per_arm[n][1] for n in sp_arms)) if sp_arms else set()
            call_sp = {s for s in call_sp if any(genus[o] in call_gn for o in by_species[s])} if call_gn else set()
            if not call_gn:
                rank = "discordant"
            elif len(call_sp) == 1:
                rank = "species"
            elif call_sp:
                rank = "species group"
            else:
                rank = "genus"
        true_sp = species.get(c, "")
        true_gn = genus.get(c, "")
        elsewhere = len(by_species[true_sp] - {c}) if true_sp else 0
        if a.holdout == "species" and true_sp:
            group = "named, species held out"
        elif true_sp and elsewhere:
            group = "named, species elsewhere"
        elif true_sp:
            group = "named, species held out"
        else:
            group = "unnamed ('sp.')"
        causes = []
        if r["copy_class"] == "divergent":
            causes.append(f"divergent ({r.get('div_kind') or 'unknown kind'})")
        if group == "named, species held out":
            causes.append("held out")
        if rank == "species":
            correct = "yes" if true_sp and call_sp == {true_sp} else ("genus only" if true_gn in call_gn else "no")
        elif rank == "species group":
            correct = "contains" if true_sp in call_sp else ("genus only" if true_gn in call_gn else "no")
        elif rank == "genus":
            correct = "genus" if true_gn in call_gn else "no"
        else:
            correct = ""
        row = {"accession": x, "culture": c, "organism": r["organism"], "genus": true_gn,
               "copy_class": r["copy_class"], "div_kind": r.get("div_kind", ""), "group": group,
               "other_cultures_of_species": elsewhere}
        for n in ARMS:
            st, sp, gn = per_arm.get(n, ("absent", set(), set()))
            row[f"{n}_state"] = st
            row[f"{n}_species"] = ";".join(sorted(sp))
            row[f"{n}_genera"] = ";".join(sorted(gn))
        row.update({"rank": rank, "call_species": ";".join(sorted(call_sp)), "call_genera": ";".join(sorted(call_gn)),
                    "correct": correct, "causes": ";".join(causes),
                    "unmatched_arms": ";".join(n for n in ARMS if per_arm.get(n, ("",))[0] == "unmatched")})
        out.append(row)

    # ---- report
    main_rows = [r for r in out if r["copy_class"] == "main"]
    print(f"{len(out)} copies ({len(main_rows)} main-type), {len(cul_meta)} cultures, {len(by_species)} named species; "
          f"arms: {', '.join(arms)}; holdout: {a.holdout}")

    print("\n1. arm states (main-type copies | divergent copies)")
    print(f"  {'arm':6s}{'matched':>18s}{'unmatched':>18s}{'no call':>18s}")
    for n in arms:
        cells = []
        for st in ("matched", "unmatched", "no call"):
            m = sum(r[f"{n}_state"] == st for r in main_rows)
            d = sum(r[f"{n}_state"] == st for r in out if r["copy_class"] == "divergent")
            cells.append(f"{m} | {d}")
        print(f"  {n:6s}" + "".join(f"{c:>18s}" for c in cells))

    def bucket(r):
        if r["rank"] == "species":
            return "sp. right" if r["correct"] == "yes" else "sp. other"
        if r["rank"] == "species group":
            return "group+" if r["correct"] == "contains" else "group-"
        if r["rank"] == "genus":
            return "genus right" if r["correct"] == "genus" else "genus wrong"
        return r["rank"]
    heads = ["sp. right", "sp. other", "group+", "group-", "genus right", "genus wrong", "discordant", "unplaced"]
    print("\n2. main-type copies: rank of the call and whether it fits the culture's name")
    print(f"  {'group':28s}{'copies':>7s}{'cult.':>6s}" + "".join(f"{h:>12s}" for h in heads))
    groups = ["named, species elsewhere", "named, species held out", "unnamed ('sp.')"]
    for g in groups:
        rs = [r for r in main_rows if r["group"] == g]
        if not rs:
            continue
        cnt = collections.Counter(bucket(r) for r in rs)
        print(f"  {g:28s}{len(rs):>7d}{len({r['culture'] for r in rs}):>6d}" + "".join(f"{cnt[h]:>12d}" for h in heads))
    print("  sp. other = one species that is not the culture's; group+/- = species group with/without it.")
    print("  expected: 'species elsewhere' at sp. right or group+; 'held out' at genus or unplaced (a species call there"
          " is a false species); 'sp.' cultures at genus, or a species call worth checking as a name")

    ref_rows = [r for r in main_rows if r["group"] == "named, species elsewhere"]
    if ref_rows:
        print("\n3. named species present elsewhere: each arm alone vs combined (main-type copies)")
        print(f"  {'':10s}{'exact species':>15s}{'group holding it':>18s}{'wrong species':>15s}{'genus/none':>12s}")
        for n in list(arms) + ["combined"]:
            ex = grp = wrong = rest = 0
            for r in ref_rows:
                if n == "combined":
                    sp = set(r["call_species"].split(";")) - {""}
                else:
                    sp = set(r[f"{n}_species"].split(";")) - {""} if r[f"{n}_state"] == "matched" else set()
                if sp == {r["organism"]}:
                    ex += 1
                elif r["organism"] in sp:
                    grp += 1
                elif sp:
                    wrong += 1
                else:
                    rest += 1
            k = len(ref_rows)
            print(f"  {n:10s}{ex:>9d} ({100 * ex / k:3.0f}%){grp:>11d} ({100 * grp / k:3.0f}%){wrong:>8d} ({100 * wrong / k:3.0f}%)"
                  f"{rest:>6d} ({100 * rest / k:3.0f}%)")

        fails = [r for r in ref_rows if r["rank"] in ("discordant", "unplaced")
                 or (r["rank"] in ("species", "species group") and r["correct"] in ("no", "genus only"))]
        by_c = collections.Counter(r["culture"] for r in fails)
        print(f"\n4. negative control (M3): main-type copies of species present elsewhere that are placed wrongly, "
              f"discordant or unplaced: {len(fails)} of {len(ref_rows)} in {len(by_c)} cultures")
        for cu, n in by_c.most_common(15):
            r0 = next(r for r in fails if r["culture"] == cu)
            arms_txt = "; ".join(f"{n2} {r0[f'{n2}_state']}" + (f" {r0[f'{n2}_species'] or r0[f'{n2}_genera']}"
                                  if r0[f'{n2}_state'] == 'matched' else '') for n2 in arms)
            print(f"      {cu} ({n}): {r0['rank']} {r0['call_species'] or r0['call_genera']} | {arms_txt}"[:230])

    held = [r for r in main_rows if r["group"] == "named, species held out" and r["rank"] in ("species", "species group")]
    if held:
        by_c = collections.Counter(r["culture"] for r in held)
        print(f"\n4b. false species: main-type copies of a species with no other culture in the reference, called as a "
              f"known species or species group: {len(held)} of "
              f"{sum(r['group'] == 'named, species held out' for r in main_rows)} in {len(by_c)} cultures")
        for cu, n in by_c.most_common(15):
            r0 = next(r for r in held if r["culture"] == cu)
            arms_txt = "; ".join(f"{n2} {r0[f'{n2}_state']}" + (f" {r0[f'{n2}_species'] or r0[f'{n2}_genera']}"
                                  if r0[f'{n2}_state'] == 'matched' else '') for n2 in arms)
            print(f"      {cu} ({n}): {r0['rank']} {r0['call_species']} | {arms_txt}"[:230])

    print("\n5. divergent-class copies (combined call)")
    for kind in ("short_5.8S", "same_length", "longer", ""):
        rs = [r for r in out if r["copy_class"] == "divergent" and r["div_kind"] == kind]
        if not rs:
            continue
        cnt = collections.Counter(r["rank"] for r in rs)
        ok = sum(r["correct"] in ("yes", "contains", "genus") for r in rs)
        wrong = sum(r["correct"] in ("no", "genus only") and r["rank"] in ("species", "species group") for r in rs)
        unm = collections.Counter(n for r in rs for n in arms if r[f"{n}_state"] == "unmatched")
        print(f"  {kind or 'no 5.8S length'}: {len(rs)} copies; " + ", ".join(f"{k} {v}" for k, v in cnt.most_common())
              + f"; consistent with the culture {ok}, wrong species {wrong}; arms unmatched: "
              + (", ".join(f"{k} {v}" for k, v in unm.most_common()) or "none"))

    # ---- 6. novelty rules
    rules = [("matched arms only", lambda um: False),
             ("demote if VT or LSU unmatched", lambda um: bool(um & {"VT", "LSU"})),
             ("demote if >= 2 arms unmatched", lambda um: len(um) >= 2),
             ("demote if any arm unmatched", lambda um: bool(um))]
    print("\n6. novelty rules on main-type copies: a species call is demoted to genus when the rule fires")
    print(f"  {'rule':32s}{'species elsewhere: kept':>26s}{'demoted':>10s}{'held out: false species':>26s}")
    ref_n = sum(r["group"] == "named, species elsewhere" for r in main_rows)
    held_n = sum(r["group"] == "named, species held out" for r in main_rows)
    for label, fires in rules:
        kept = demoted = false_sp = 0
        for r in main_rows:
            um = set(r["unmatched_arms"].split(";")) - {""}
            is_sp = r["rank"] in ("species", "species group")
            dem = is_sp and fires(um)
            if r["group"] == "named, species elsewhere":
                if is_sp and not dem and r["correct"] in ("yes", "contains"):
                    kept += 1
                demoted += dem
            elif r["group"] == "named, species held out" and is_sp and not dem:
                false_sp += 1
        k = f"{kept}/{ref_n} ({100 * kept / ref_n:.0f}%)" if ref_n else "-"
        f = f"{false_sp}/{held_n} ({100 * false_sp / held_n:.0f}%)" if held_n else "-"
        print(f"  {label:32s}{k:>26s}{demoted:>10d}{f:>26s}")

    write_tsv(f"{a.out}.copies.tsv", out, list(out[0]))
    print(f"\ntable: {a.out}.copies.tsv")


if __name__ == "__main__":
    sys.exit(main())
