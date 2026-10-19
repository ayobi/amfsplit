"""Fixed Phase 4 decision rules from the prior baseline-consistency experiment.

Decision block copied from novelty_baseline_consistent.py, with training-only
arm inputs and an SSU nearest-neighbor cache. No thresholds fitted here.
"""
from collections import defaultdict
from functools import lru_cache
from types import SimpleNamespace
import statistics
from novelty_baseline_consistent import dist as original_dist, split

dist = lru_cache(maxsize=100000)(original_dist)


def evaluate(data, query_ids, fits, drop, rec):
    cop, fa = data.cop, data.fa
    main_58 = data.main58
    vt_units, vt_best, vt_anchor, vt_holders = fits["VT"].sets, fits["VT"].best, fits["VT"].anchors, fits["VT"].holders
    arm_anchor = {n:f.anchors for n,f in fits.items()}
    para_arms = ["VT", "SH"]
    reach = {n: (defaultdict(set), defaultdict(set)) for n in para_arms}
    for n in para_arms:
        for y,u in arm_anchor[n].items():
            if not u or cop[y]["culture"] in drop: continue
            if cop[y]["copy_class"] == "main": reach[n][0][u].add(cop[y]["culture"])
            elif cop[y].get("div_kind") == "short_5.8S": reach[n][1][u].add(cop[y]["culture"])
    a = SimpleNamespace(artifact_min_ssu=1., artifact_ratio=.5, short_58s=3, subst_58s=3., vt_min_id=97., boundary=1.)
    out=[]
    for x in query_ids:
        r=cop[x]; c=r["culture"]
        reasons,flags=[],[]
        near,near_d,nearest=data.nearest(x, drop)
        its_d = dist(fa["full_its"][x], fa["full_its"][near]) if near and x in fa["full_its"] and near in fa["full_its"] else None

        # 1. artifact
        # SSU often has ties. Require complete, concordant evidence across ALL
        # tied nearest copies before the A1 screen fires. The representative
        # pair remains in the legacy report columns; it is not the sole test.
        artifact = False
        artifact_its = []
        if near_d is not None and near_d >= a.artifact_min_ssu:
            artifact_its = [dist(fa["full_its"][x], fa["full_its"][y])
                            for y in nearest if x in fa["full_its"] and y in fa["full_its"]]
            votes = [near_d >= a.artifact_ratio * d for d in artifact_its]
            artifact = len(votes) == len(nearest) and bool(votes) and all(votes)
            if votes and any(votes) and (not all(votes) or len(votes) != len(nearest)):
                flags.append("A2 artifact screen unresolved across SSU-tied nearest copies")
        if artifact:
            reasons.append(f"A1 screening flag: SSU {near_d:.1f}% vs ITS "
                           f"{min(artifact_its):.1f}-{max(artifact_its):.1f}% across "
                           f"{len(nearest)} SSU-tied nearest copies; not proof of an artifact")

        # 2. divergent class from the copy's own sequence
        short = paralog = False
        ncul = cop[near]["culture"] if near else ""
        comp = sorted(cu for cu in vt_holders.get(vt_anchor.get(x, ""), set()) - drop if main_58.get(cu))
        comp_src = "cultures sharing its VT"
        if not comp and near_d is not None and near_d <= 1.0:
            comp = sorted({cop[y]["culture"] for y in nearest if main_58.get(cop[y]["culture"])})
            comp_src = f"SSU-tied nearest cultures ({near_d:.2f}% SSU)"
        s58_info = ""
        if r.get("s58_len") and comp and x in fa["s58"]:
            ref_len = statistics.median(statistics.median(int(cop[y]["s58_len"]) for y in main_58[cu]) for cu in comp)
            delta = int(r["s58_len"]) - ref_len
            # Development ablation: when testing a short query against ref_len,
            # the sequence comparator must represent that length baseline too.
            # Do not let a shorter main-labelled sequence supply a different
            # baseline for the sequence-distance half of the SAME D1 test.
            compare_ids = [y for cu in comp for y in main_58[cu]]
            if delta <= -a.short_58s:
                compare_ids = [y for y in compare_ids if len(fa["s58"][y]) >= ref_len]
            d58 = min(dist(fa["s58"][x], fa["s58"][y]) for y in compare_ids)
            s58_info = f"5.8S {r['s58_len']} bp vs {ref_len:g} in {comp_src}, {d58:.1f}% different"
            if delta <= -a.short_58s and d58 >= a.subst_58s:
                short = True
                reasons.append(f"D1 {s58_info}")
            elif delta <= -a.short_58s:
                flags.append(f"D4 {s58_info} (shorter, below divergence threshold: inspect boundaries)")
            elif abs(delta) < a.short_58s and d58 >= a.subst_58s:
                flags.append(f"D3 {s58_info}")
        for n in para_arms:
            u = arm_anchor[n].get(x, "")
            if not u:
                continue
            m_c = reach[n][0].get(u, set()) - drop
            d_c = reach[n][1].get(u, set()) - drop
            if len(d_c) >= 2 and not m_c:
                paralog = True
                reasons.append(f"D2 paralog unit {u} ({n}; short-5.8S copies of {len(d_c)} other cultures)")
        divergent = short or paralog

        # 3. boundary
        rc = rec[x]
        vt_state = rc["VT_state"]
        best_vt, best_id = vt_best.get(x, ("", 0.0))
        boundary = bool(best_vt) and vt_state == "no call" and a.vt_min_id - a.boundary <= best_id < a.vt_min_id
        off_vt = bool(best_vt) and vt_state == "no call" and best_id < a.vt_min_id - a.boundary
        if vt_state == "no call" and best_vt and best_id >= a.vt_min_id:
            flags.append("B2 VT identity passes threshold but no call: inspect coverage, filtering and provenance")
        if vt_state == "no call" and not best_vt:
            flags.append("Q1 no recorded best VT hit; distance from all VTs is not established")

        # 4. placement and novelty
        rank, um = rc["rank"], set(split(rc["unmatched_arms"], ";"))
        if artifact:
            status = "artifact"
        elif divergent:
            status = "divergent class"
            place = rc["VT_species"] or rc["VT_genera"]
            reasons.append(f"placed by VT: {place or 'no VT match'}")
        elif boundary:
            status = "boundary"
            reasons.append(f"B1 best VT {best_vt} {best_id:.2f}%")
        elif rank in ("species", "species group"):
            if um & {"VT", "LSU"}:
                status = "novel species, genus known"
                reasons += [f"{'N1' if n == 'VT' else 'N2'} {n} unmatched" for n in ("VT", "LSU") if n in um]
            else:
                status = "known " + rank
                reasons.append("K1")
        elif rank == "genus":
            status = "novel species, genus known"
            reasons.append("N3 no species shared by the arms")
        elif rank == "discordant":
            status = "discordant"
            reasons.append("X1 arms share no genus")
        else:
            status = "novel lineage"
            reasons.append("N5 recorded best VT identity below boundary" if off_vt else "N4 no arm matched")
        placement = (rc["call_species"] if status.startswith("known") else
                     (rc["VT_species"] or rc["VT_genera"]) if status == "divergent class" else
                     rc["call_genera"] if status == "novel species, genus known" else "")
        row = {"accession": x, "culture": c, "organism": r["organism"], "copy_class": r["copy_class"],
                    "div_kind": r.get("div_kind", ""), "group": rc["group"], "status": status,
                    "novelty_claim": "candidate_only" if status.startswith("novel") else "not_asserted",
                    "artifact_claim": "screen_flag_only" if status == "artifact" else "not_asserted",
                    "placement": placement, "reasons": "; ".join(reasons), "flags": "; ".join(flags),
                    "s58": s58_info, "nearest_culture": ncul, "nearest_ssu_pct": round(near_d, 2) if near_d is not None else "",
                    "nearest_its_pct": round(its_d, 2) if its_d is not None else "",
                    "nearest_ssu_copy_count": len(nearest),
                    "nearest_ssu_cultures": ";".join(sorted({cop[y]["culture"] for y in nearest})),
                    "artifact_comparator_its_pct": ";".join(f"{d:.8g}" for d in sorted(artifact_its)),
                    "best_vt": best_vt, "best_vt_pident": best_id, "reconcile_rank": rank,
                    "unmatched_arms": rc["unmatched_arms"]}
        row.update(sequence_type_signal="short_58S" if short else "paralog_unit" if paralog else "unresolved",
                   d1_comparator_cultures=";".join(comp), d1_comparator_copies=";".join(compare_ids) if r.get("s58_len") and comp and x in fa["s58"] else "",
                   d1_short=int(short), d2_paralog=int(paralog), artifact_screen=int(artifact),
                   placement_rank_before_type_filter=rc["rank"], placement_species_before_type_filter=rc["call_species"])
        out.append(row)

    return out
