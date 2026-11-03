//! Protocol v3 rules under test (opt-in; the frozen method stays the default).
//!
//! R2, the inside gate: after N6, a known species or species group stands only if the sequence lies
//! inside the called species, or one member of the called group, on the markers it holds. Inside means
//! within the marker's within-species limit of the nearest main-type sequence of that species among
//! the training cultures, in full ITS, LR0R-FLR2 and the SSU flank. The limits are the complex rule's
//! (`species_groups.py`): per marker, the 99th percentile (exclusive method) of each unique
//! sequence-culture pair's distance to the nearest pair of another training culture of the same named
//! species, over species with two or more training cultures, refitted per round. Markers the sequence
//! lacks are reported and left unchecked; a call checked on one marker or none fails (decided
//! 2 November 2026). A failing call becomes "novel species, genus known", code N7.
//!
//! R1a, unnamed cultures as taxa (inside): after the decision steps (and before R2 when both run), a
//! sequence inside at least one unnamed training culture on all three markers, and inside no named
//! species, is called "unnamed culture" (code U1): the cultures' genus and codes, no species.
//! Divergent, artifact and discordant calls keep their call and gain flag U. A known call inside both
//! a named species and an unnamed culture keeps its species, written as a group that includes the
//! culture. Unnamed means a culture without a usable species name, as in the frozen data.
//!
//! R1b, unnamed cultures as taxa (block), the paper's alternative: a matched VT, SH or LSU unit held
//! only by unnamed training cultures blocks a species. A known call becomes genus only, as N3 does
//! ("novel species, genus known"), with code U2 and flag U naming the cultures. The holders are the
//! naming code's own fits per round: the VT and SH trainers, and the LSU placement fit on the
//! species-clade units of the round's placements.
//!
//! R5, the Rhizoglomus motif step: a new decision step D4, after D1 and D2 and before B1. A sequence
//! whose 5.8S (the ITSxRust window) contains YTYGGRGARTATG exactly becomes a divergent sequence,
//! named by its VT alone, when every matched VT unit is linked only to Rhizoglomus training cultures
//! (v4 labels Rhizophagus and Rhizoglomus). The SSU-fallback variant (r5s) also accepts a sequence
//! whose VT is not matched when its nearest training culture by SSU flank is Rhizoglomus and within
//! the within-genus SSU limit (R3a's definition, refitted per round). A motif without that evidence
//! is flag M1 and changes nothing. Artifact and divergent calls (A1, D1, D2) come first.
//!
//! R3, the provisional-genus flag: G1, which never changes a call, on any call that gives a genus.
//! R3a (distance) flags when the sequence lies beyond the within-genus limit from every training
//! culture of the given genus on any marker (r3a) or on all markers (r3a_all, the sensitivity).
//! R3b (placement) flags when the best LSU placement's pendant length exceeds R3B_CUT, the 99th
//! percentile of that length among the 803 culture main-type sequences of the 37 LSU-removal rounds.
//!
//! R7, the genus-clade fallback: a new-lineage call takes the genus of the V18 genus clade its LSU
//! places in, when that clade's unit is held by two or more training cultures, all of one genus that
//! agrees with the clade's label (code G2; never a species). Holders are training cultures whose
//! bridge-eligible main-type copies place in the unit, as the frozen LSU arm counts them. It runs only
//! with an R3 arm, after R2 and before the flag.
//!
//! R8, novelty on one marker: a new lineage resting on one unmatched arm, with the other two arms no
//! call, keeps its call and carries flag N8 (r8). The variant R8d (r8d) demotes it to "insufficient
//! evidence", neither known nor new. It runs after R7 and before the R3 flag.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

use anyhow::{anyhow, bail, Result};

use super::*;

pub const INSIDE_MARKERS: [(&str, &str); 3] = [("full_its", "ITS"), ("LSU_LROR-FLR2", "LSU"), ("ssu_flank", "SSU")];
pub const R2_COLUMNS: [&str; 4] = ["r2_inside", "r2_failing_markers", "r2_unchecked_markers", "r2_limits_pct"];
pub const R1A_COLUMNS: [&str; 3] = ["r1a_unnamed_inside", "r1a_named_inside", "r1a_action"];
pub const R1B_COLUMNS: [&str; 2] = ["r1b_blocking_units", "r1b_unnamed_holders"];
pub const R5_COLUMNS: [&str; 3] = ["r5_motif", "r5_rhizoglomus_evidence", "r5_action"];
pub const R3_COLUMNS: [&str; 3] = ["r3_genus", "r3_evidence", "r3_flag"];
pub const R7_COLUMNS: [&str; 2] = ["r7_genus_clade_unit", "r7_action"];
pub const R8_COLUMNS: [&str; 1] = ["r8_action"];
/// In protocol order: the decision step R5 (or its SSU-fallback variant r5s), then R1a or R1b, then
/// R2; the R3 flags (r3a, r3a_all, r3b) read the final calls.
pub const RULES: [&str; 11] = ["r5", "r5s", "r1a", "r1b", "r2", "r7", "r8", "r8d", "r3a", "r3a_all", "r3b"];

/// R8 (flag) or R8d (demote) on rows; returns the number of calls flagged or demoted.
fn apply_r8(demote: bool, rows: &mut [Row]) -> Result<usize> {
    let mut n = 0;
    for row in rows {
        row.insert("r8_action".into(), String::new());
        if row_value(row, "status")? != "novel lineage" { continue; }
        let states: Vec<(&str, String)> = ["VT", "SH", "LSU"].iter().map(|a| Ok((*a, row_value(row, &format!("{a}_state"))?.to_owned()))).collect::<Result<_>>()?;
        let unmatched: Vec<&str> = states.iter().filter(|(_, s)| s == "unmatched").map(|(a, _)| *a).collect();
        let no_call = states.iter().filter(|(_, s)| s == "no call").count();
        if unmatched.len() != 1 || no_call != 2 { continue; }
        let text = format!("N8 new lineage resting on one unmatched arm ({}); the other two arms make no call", unmatched[0]);
        let col = if demote { "reasons" } else { "flags" };
        let existing = row_value(row, col)?.to_owned();
        row.insert(col.into(), if existing.is_empty() { text } else { format!("{existing}; {text}") });
        if demote {
            row.insert("status".into(), "insufficient evidence".into());
            row.insert("placement".into(), String::new());
            row.insert("novelty_claim".into(), "not_asserted".into());
            row.insert("r8_action".into(), "R8d".into());
        } else {
            row.insert("r8_action".into(), "N8".into());
        }
        n += 1;
    }
    Ok(n)
}

/// R7's genus-clade unit holders for one round: unit -> training cultures.
fn genus_clade_holders(data: &FrozenData, training: &StringSet, placements: &BTreeMap<String, Row>) -> BTreeMap<String, StringSet> {
    let mut h: BTreeMap<String, StringSet> = BTreeMap::new();
    for (acc, c) in &data.copies {
        if c.copy_class != "main" || !training.contains(&c.culture) { continue; }
        let Some(p) = placements.get(acc) else { continue };
        if p.get("bridge_alignment_eligible").map(String::as_str) != Some("1") { continue; }
        if let Some(u) = p.get("external_genus_clade_unit").filter(|u| !u.is_empty()) { h.entry(u.clone()).or_default().insert(c.culture.clone()); }
    }
    h
}

/// R7 on rows; returns the number of genera given.
fn apply_r7(data: &FrozenData, holders: &BTreeMap<String, StringSet>, rows: &mut [Row], placement: &dyn Fn(&str) -> Option<Row>) -> Result<usize> {
    let mut given = 0;
    for row in rows {
        for c in R7_COLUMNS { row.insert(c.into(), String::new()); }
        if row_value(row, "status")? != "novel lineage" { continue; }
        let acc = row_value(row, "accession")?.to_owned();
        let Some(p) = placement(&acc) else { continue };
        let unit = p.get("external_genus_clade_unit").cloned().unwrap_or_default();
        if unit.is_empty() { continue; }
        row.insert("r7_genus_clade_unit".into(), unit.clone());
        let label = normalise_genus(p.get("external_genus_clade").map(String::as_str).unwrap_or(""));
        let hs = holders.get(&unit).cloned().unwrap_or_default();
        let genera: BTreeSet<String> = hs.iter().map(|c| data.cultures.get(c).map(|m| m.genus.clone()).unwrap_or_default()).collect();
        if hs.len() < 2 || genera.len() != 1 || !genera.contains(&label) { row.insert("r7_action".into(), "abstain".into()); continue; }
        row.insert("status".into(), "novel species, genus known".into());
        row.insert("call_genera".into(), label.clone());
        row.insert("placement".into(), label.clone());
        let reason = format!("G2 genus {label} from the V18 genus clade {unit}, held by {} training cultures of {label}", hs.len());
        let existing = row_value(row, "reasons")?.to_owned();
        row.insert("reasons".into(), if existing.is_empty() { reason } else { format!("{existing}; {reason}") });
        row.insert("r7_action".into(), "G2".into());
        given += 1;
    }
    Ok(given)
}
/// R3b's cut, by its definition (Protocol v3): the 99th percentile (exclusive) of the best LSU
/// placement's pendant length among the 803 main-type query rows of the 37 species LSU-removal rounds
/// (release 0.3.2, external_holdout_v2). A test re-derives it from the published rounds.
pub const R3B_CUT: f64 = 0.074_880_850_445;
const GENUS_MARKERS: [(&str, &str); 3] = [("ssu_flank", "SSU"), ("LSU_LROR-FLR2", "LSU"), ("full_its", "ITS")];

/// R3a for one round: within-genus limits per marker, and each genus's training main-type sequences.
pub(crate) struct GenusRule {
    limits: BTreeMap<&'static str, Option<f64>>,
    seqs: BTreeMap<(&'static str, String), Vec<String>>,
}

impl GenusRule {
    pub(crate) fn fit(data: &FrozenData, training: &StringSet, memo: &mut DistanceMemo) -> Result<Self> {
        let mut limits = BTreeMap::new();
        for (m, _) in GENUS_MARKERS { limits.insert(m, within_genus_limit(data, training, m, memo)?); }
        let mut seqs: BTreeMap<(&'static str, String), Vec<String>> = BTreeMap::new();
        for (acc, c) in &data.copies {
            if c.copy_class != "main" || !training.contains(&c.culture) { continue; }
            let genus = data.cultures.get(&c.culture).map(|m| m.genus.clone()).unwrap_or_default();
            if genus.is_empty() { continue; }
            for (m, _) in GENUS_MARKERS {
                if let Some(s) = data.fa.get(m).and_then(|f| f.get(acc)) { seqs.entry((m, genus.clone())).or_default().push(s.clone()); }
            }
        }
        Ok(Self { limits, seqs })
    }

    /// The markers on which the sequence lies beyond the limit from every training culture of
    /// `genera`, and how many could be checked (None when none could).
    fn beyond(&self, genera: &[&str], query: &dyn Fn(&str) -> Option<String>, memo: &mut DistanceMemo) -> Result<Option<(Vec<&'static str>, usize)>> {
        let (mut beyond, mut checked) = (Vec::new(), 0);
        for (m, lab) in GENUS_MARKERS {
            let (Some(w), Some(t)) = (query(m).filter(|w| !w.is_empty()), self.limits[m]) else { continue };
            let mut best = f64::INFINITY;
            for g in genera { for s in self.seqs.get(&(m, (*g).to_owned())).into_iter().flatten() { best = best.min(memo.d(&w, s)?); } }
            if !best.is_finite() { continue; }
            checked += 1;
            if best > t { beyond.push(lab); }
        }
        Ok((checked > 0).then_some((beyond, checked)))
    }
}

/// R3 on rows (a flag only): `variant` is r3a, r3a_all or r3b. Returns the number of calls flagged.
pub(crate) fn apply_r3(variant: &str, genus_rule: Option<&GenusRule>, rows: &mut [Row], query: &dyn Fn(&str, &str) -> Option<String>, memo: &mut DistanceMemo) -> Result<usize> {
    let mut flagged = 0;
    for row in rows {
        for c in R3_COLUMNS { row.insert(c.into(), String::new()); }
        let genera_text = row_value(row, "call_genera")?.to_owned();
        let genera: Vec<&str> = genera_text.split(';').filter(|g| !g.is_empty()).collect();
        if genera.is_empty() { continue; }
        row.insert("r3_genus".into(), genera_text.clone());
        let acc = row_value(row, "accession")?.to_owned();
        let (flag, evidence) = match variant {
            "r3b" => {
                let p = row_value(row, "lsu_best_pendant_length")?.to_owned();
                match p.parse::<f64>() {
                    Ok(v) => (v > R3B_CUT, format!("LSU pendant length {p} against {R3B_CUT}")),
                    Err(_) => (false, "no LSU placement".into()),
                }
            }
            _ => {
                let rule = genus_rule.ok_or_else(|| anyhow!("R3a needs the round's genus rule"))?;
                match rule.beyond(&genera, &|m| query(m, &acc), memo)? {
                    None => (false, "no marker to check".into()),
                    Some((beyond, checked)) => {
                        let fire = if variant == "r3a_all" { beyond.len() == checked } else { !beyond.is_empty() };
                        (fire, if beyond.is_empty() { format!("inside the genus limit on {checked} marker(s)") } else { format!("beyond the genus limit on {}", beyond.join(", ")) })
                    }
                }
            }
        };
        row.insert("r3_evidence".into(), evidence.clone());
        if !flag { continue; }
        row.insert("r3_flag".into(), "G1".into());
        let existing = row_value(row, "flags")?.to_owned();
        let text = format!("G1 provisional genus {genera_text}: {evidence}");
        row.insert("flags".into(), if existing.is_empty() { text } else { format!("{existing}; {text}") });
        flagged += 1;
    }
    Ok(flagged)
}
const RHIZOGLOMUS: [&str; 2] = ["Rhizophagus", "Rhizoglomus"];

/// YTYGGRGARTATG, exactly (IUPAC: Y = C or T, R = A or G).
pub(crate) fn has_rhizoglomus_motif(s58: &str) -> bool {
    const M: &[u8] = b"YTYGGRGARTATG";
    let ok = |p: u8, b: u8| match p { b'Y' => b == b'C' || b == b'T', b'R' => b == b'A' || b == b'G', _ => p == b };
    let b: Vec<u8> = s58.bytes().map(|c| c.to_ascii_uppercase()).collect();
    b.windows(M.len()).any(|w| w.iter().zip(M).all(|(x, p)| ok(*p, *x)))
}

/// R3a's within-genus limit per marker: the 99th percentile (exclusive) of the distance from each
/// main-type sequence of a named training species to the nearest main-type sequence of another named
/// training species of its genus, over genera with two or more named training species.
pub(crate) fn within_genus_limit(data: &FrozenData, training: &StringSet, marker: &str, memo: &mut DistanceMemo) -> Result<Option<f64>> {
    let mut by_genus: BTreeMap<String, Vec<(String, String)>> = BTreeMap::new();
    for (acc, c) in &data.copies {
        if c.copy_class != "main" || !training.contains(&c.culture) { continue; }
        let (Some(sp), Some(seq)) = (data.species.get(&c.culture), data.fa.get(marker).and_then(|f| f.get(acc))) else { continue };
        let genus = data.cultures.get(&c.culture).map(|m| m.genus.clone()).unwrap_or_default();
        by_genus.entry(genus).or_default().push((seq.clone(), sp.clone()));
    }
    let mut vals = Vec::new();
    for items in by_genus.values() {
        if items.iter().map(|(_, sp)| sp).collect::<BTreeSet<_>>().len() < 2 { continue; }
        for (s1, sp1) in items {
            let mut best = f64::INFINITY;
            for (s2, sp2) in items { if sp2 != sp1 { best = best.min(memo.d(s1, s2)?); } }
            vals.push(best);
        }
    }
    Ok((vals.len() >= 2).then(|| quantiles_exclusive(&vals, 100)[98]))
}

/// R5 for one round: the SSU-fallback limit and the training main-type SSU sequences.
pub(crate) struct MotifRule {
    fallback: bool,
    ssu_limit: Option<f64>,
    ssu_refs: Vec<(String, String)>,
}

impl MotifRule {
    pub(crate) fn fit(data: &FrozenData, training: &StringSet, fallback: bool, memo: &mut DistanceMemo) -> Result<Self> {
        let ssu_limit = if fallback { within_genus_limit(data, training, "ssu_flank", memo)? } else { None };
        let ssu_refs = data.copies.iter().filter(|(_, c)| c.copy_class == "main" && training.contains(&c.culture))
            .filter_map(|(acc, c)| data.fa.get("ssu_flank").and_then(|f| f.get(acc)).map(|s| (s.clone(), c.culture.clone()))).collect();
        Ok(Self { fallback, ssu_limit, ssu_refs })
    }

    /// Returns the number of calls turned divergent (D4).
    pub(crate) fn apply(&self, data: &FrozenData, rows: &mut [Row], query: &dyn Fn(&str, &str) -> Option<String>, memo: &mut DistanceMemo) -> Result<usize> {
        let mut turned = 0;
        let genus = |c: &str| data.cultures.get(c).map(|m| m.genus.clone()).unwrap_or_default();
        for row in rows {
            for c in R5_COLUMNS { row.insert(c.into(), String::new()); }
            let acc = row_value(row, "accession")?.to_owned();
            let Some(s58) = query("s58", &acc) else { continue };
            if !has_rhizoglomus_motif(&s58) { row.insert("r5_motif".into(), "0".into()); continue; }
            row.insert("r5_motif".into(), "1".into());
            let status = row_value(row, "status")?.to_owned();
            let vt_genera: BTreeSet<String> = row_value(row, "VT_genera")?.split(';').filter(|g| !g.is_empty()).map(str::to_owned).collect();
            let vt_matched = row_value(row, "VT_state")? == "matched";
            let mut evidence: Option<(String, String)> = None;   // (kind, placement)
            if vt_matched && !vt_genera.is_empty() && vt_genera.iter().all(|g| RHIZOGLOMUS.contains(&g.as_str())) {
                evidence = Some(("VT".into(), row_value(row, "VT_species")?.to_owned()));
            } else if self.fallback && !vt_matched {
                if let (Some(w), Some(limit)) = (query("ssu_flank", &acc).filter(|w| !w.is_empty()), self.ssu_limit) {
                    let mut best = f64::INFINITY;
                    let mut nearest: BTreeSet<String> = BTreeSet::new();
                    for (s, c) in &self.ssu_refs {
                        let d = memo.d(&w, s)?;
                        if d < best { best = d; nearest.clear(); }
                        if d == best { nearest.insert(c.clone()); }
                    }
                    if !nearest.is_empty() && best <= limit && nearest.iter().all(|c| RHIZOGLOMUS.contains(&genus(c).as_str())) {
                        let genera: BTreeSet<String> = nearest.iter().map(|c| genus(c)).collect();
                        evidence = Some((format!("SSU fallback: nearest culture at {best:.2}%, within the within-genus SSU limit {limit:.3}%"),
                                         genera.into_iter().collect::<Vec<_>>().join(";")));
                    }
                }
            }
            let add = |row: &mut Row, col: &str, text: String| -> Result<()> {
                let existing = row_value(row, col)?.to_owned();
                row.insert(col.into(), if existing.is_empty() { text } else { format!("{existing}; {text}") });
                Ok(())
            };
            match evidence {
                Some((kind, placement)) if status != "artifact" && status != "divergent class" => {
                    row.insert("status".into(), "divergent class".into());
                    row.insert("placement".into(), placement.clone());
                    row.insert("novelty_claim".into(), "not_asserted".into());
                    add(row, "reasons", format!("D4 5.8S carries the Rhizoglomus motif YTYGGRGARTATG; Rhizoglomus by {kind}; placed by {}: {placement}",
                                                if kind == "VT" { "VT" } else { "genus" }))?;
                    row.insert("r5_rhizoglomus_evidence".into(), kind);
                    row.insert("r5_action".into(), "D4".into());
                    turned += 1;
                }
                Some((kind, _)) => { row.insert("r5_rhizoglomus_evidence".into(), kind); }
                None => {
                    add(row, "flags", "M1 5.8S carries the Rhizoglomus motif without SSU evidence for Rhizoglomus".into())?;
                    row.insert("r5_action".into(), "M1".into());
                }
            }
        }
        Ok(turned)
    }
}

/// R1b's unit holders for one round: arm -> unit -> training cultures holding it.
pub(crate) struct BlockRule {
    holders: BTreeMap<&'static str, BTreeMap<String, StringSet>>,
}

impl BlockRule {
    pub(crate) fn fit(data: &FrozenData, trainers: &BTreeMap<String, Trainer>, psets: &BTreeMap<String, StringSet>, training: &StringSet) -> Result<Self> {
        let mut holders = BTreeMap::new();
        for arm in ["VT", "SH"] {
            let t = trainers.get(arm).ok_or_else(|| anyhow!("no {arm} trainer"))?;
            holders.insert(arm, t.fit(data, training)?.holders);
        }
        holders.insert("LSU", placement_fit(data, training, psets).holders);
        Ok(Self { holders })
    }

    /// Returns the number of known calls blocked.
    pub(crate) fn apply(&self, data: &FrozenData, rows: &mut [Row]) -> Result<usize> {
        let code = |c: &str| c.split(" | ").nth(1).unwrap_or(c).to_owned();
        let mut blocked = 0;
        for row in rows {
            for c in R1B_COLUMNS { row.insert(c.into(), String::new()); }
            let status = row_value(row, "status")?.to_owned();
            if status != "known species" && status != "known species group" { continue; }
            let (mut units, mut cultures) = (Vec::new(), BTreeSet::new());
            for arm in ["VT", "SH", "LSU"] {
                if row_value(row, &format!("{arm}_state"))? != "matched" { continue; }
                for u in row_value(row, &format!("{arm}_units"))?.split(',').filter(|u| !u.is_empty()) {
                    let Some(h) = self.holders[arm].get(u) else { continue };
                    if !h.is_empty() && h.iter().all(|c| !data.species.contains_key(c)) {
                        units.push(format!("{arm} {u}"));
                        cultures.extend(h.iter().map(|c| code(c)));
                    }
                }
            }
            if units.is_empty() { continue; }
            let codes = cultures.into_iter().collect::<Vec<_>>().join(", ");
            row.insert("r1b_blocking_units".into(), units.join("; "));
            row.insert("r1b_unnamed_holders".into(), codes.clone());
            let genus = row_value(row, "call_genera")?.to_owned();
            row.insert("status".into(), "novel species, genus known".into());
            row.insert("placement".into(), genus);
            row.insert("novelty_claim".into(), "candidate_only".into());
            for (col, text) in [("reasons", format!("U2 species blocked: {} held only by unnamed training culture(s) {codes} (genus only, as N3)", units.join("; "))),
                                ("flags", format!("U unit held only by unnamed culture(s) {codes}"))] {
                let existing = row_value(row, col)?.to_owned();
                row.insert(col.into(), if existing.is_empty() { text } else { format!("{existing}; {text}") });
            }
            blocked += 1;
        }
        Ok(blocked)
    }
}

/// Unique main-type sequences per (marker, species), training cultures only, each with its culture.
fn species_items(data: &FrozenData, training: &StringSet) -> BTreeMap<(&'static str, String), BTreeSet<(String, String)>> {
    let mut items: BTreeMap<(&'static str, String), BTreeSet<(String, String)>> = BTreeMap::new();
    for (acc, c) in &data.copies {
        if c.copy_class != "main" || !training.contains(&c.culture) { continue; }
        let Some(species) = data.species.get(&c.culture) else { continue };
        for (m, _) in INSIDE_MARKERS {
            if let Some(seq) = data.fa.get(m).and_then(|f| f.get(acc)) {
                items.entry((m, species.clone())).or_default().insert((seq.clone(), c.culture.clone()));
            }
        }
    }
    items
}

/// Every distance R2 computes, once: sequences are interned, and a pair keeps its argument order
/// (the distance's direction matters for equal lengths, as in the frozen code).
#[derive(Default)]
pub struct DistanceMemo {
    ids: HashMap<String, usize>,
    memo: HashMap<(usize, usize), f64>,
}

impl DistanceMemo {
    fn id(&mut self, s: &str) -> usize {
        if let Some(i) = self.ids.get(s) { return *i; }
        let i = self.ids.len();
        self.ids.insert(s.to_owned(), i);
        i
    }
    fn d(&mut self, a: &str, b: &str) -> Result<f64> {
        let key = (self.id(a), self.id(b));
        if let Some(v) = self.memo.get(&key) { return Ok(*v); }
        let v = percent_distance(a, b, -1)?.ok_or_else(|| anyhow!("no distance"))?;
        self.memo.insert(key, v);
        Ok(v)
    }
}

pub struct InsideRule {
    items: BTreeMap<(&'static str, String), BTreeSet<(String, String)>>,
    /// unnamed training cultures: (marker, culture) -> main-type sequences
    unnamed: BTreeMap<(&'static str, String), BTreeSet<String>>,
    /// per marker: the limit, or None when fewer than two values exist
    pub limits: BTreeMap<&'static str, Option<f64>>,
}

impl InsideRule {
    /// species_groups.py: the 99th percentile of the nearest other-culture distance within a species.
    pub(crate) fn fit(data: &FrozenData, training: &StringSet, memo: &mut DistanceMemo) -> Result<Self> {
        let items = species_items(data, training);
        let mut unnamed: BTreeMap<(&'static str, String), BTreeSet<String>> = BTreeMap::new();
        for (acc, c) in &data.copies {
            if c.copy_class != "main" || !training.contains(&c.culture) || data.species.contains_key(&c.culture) { continue; }
            for (m, _) in INSIDE_MARKERS {
                if let Some(seq) = data.fa.get(m).and_then(|f| f.get(acc)) { unnamed.entry((m, c.culture.clone())).or_default().insert(seq.clone()); }
            }
        }
        let mut limits = BTreeMap::new();
        for (m, _) in INSIDE_MARKERS {
            let mut vals = Vec::new();
            for ((mm, _), set) in &items {
                if *mm != m || set.iter().map(|(_, c)| c).collect::<BTreeSet<_>>().len() < 2 { continue; }
                for (s1, c1) in set {
                    let mut best: Option<f64> = None;
                    for (s2, c2) in set {
                        if c2 == c1 { continue; }
                        let d = memo.d(s1, s2)?;
                        best = Some(best.map_or(d, |b: f64| b.min(d)));
                    }
                    if let Some(b) = best { vals.push(b); }
                }
            }
            limits.insert(m, (vals.len() >= 2).then(|| quantiles_exclusive(&vals, 100)[98]));
        }
        Ok(Self { items, unnamed, limits })
    }

    /// Distance from a query window to the nearest training main-type sequence of a species.
    fn nearest(&self, marker: &'static str, species: &str, query: &str, memo: &mut DistanceMemo) -> Result<Option<f64>> {
        let Some(set) = self.items.get(&(marker, species.to_owned())) else { return Ok(None) };
        let mut best: Option<f64> = None;
        for (s, _) in set {
            let d = memo.d(query, s)?;
            best = Some(best.map_or(d, |b: f64| b.min(d)));
        }
        Ok(best)
    }

    fn limits_text(&self) -> String {
        INSIDE_MARKERS.iter().map(|(m, lab)| format!("{lab} {}", self.limits[m].map(|t| py_float(py_round_float(t, 3))).unwrap_or_default()))
            .collect::<Vec<_>>().join("; ")
    }

    /// Inside on all three markers: each window present and within the marker's limit of the nearest
    /// of `seqs(marker)`.
    fn inside_all(&self, windows: &[(&'static str, Option<String>)], seqs: &dyn Fn(&'static str) -> Vec<String>, memo: &mut DistanceMemo) -> Result<bool> {
        for (m, w) in windows {
            let (Some(w), Some(t)) = (w, self.limits[m]) else { return Ok(false) };
            let set = seqs(m);
            if set.is_empty() { return Ok(false); }
            let mut best = f64::INFINITY;
            for s in &set { best = best.min(memo.d(w, s)?); }
            if best > t { return Ok(false); }
        }
        Ok(true)
    }

    /// R1a on rows named with these training cultures. Returns the number of U1 calls.
    pub(crate) fn apply_r1a(&self, data: &FrozenData, rows: &mut [Row], query: &dyn Fn(&str, &str) -> Option<String>, memo: &mut DistanceMemo) -> Result<usize> {
        let species: BTreeSet<&String> = self.items.keys().map(|(_, s)| s).collect();
        let cultures: BTreeSet<&String> = self.unnamed.keys().map(|(_, c)| c).collect();
        let code = |c: &str| c.split(" | ").nth(1).unwrap_or(c).to_owned();
        let label = |c: &str| c.replace(" | ", " ");
        let mut u1 = 0;
        for row in rows {
            for c in R1A_COLUMNS { row.insert(c.into(), String::new()); }
            let acc = row_value(row, "accession")?.to_owned();
            let windows: Vec<(&'static str, Option<String>)> = INSIDE_MARKERS.iter().map(|(m, _)| (*m, query(m, &acc).filter(|w| !w.is_empty()))).collect();
            if windows.iter().any(|(_, w)| w.is_none()) { continue; }
            let mut inside_u = Vec::new();
            for c in &cultures {
                let seqs = |m: &'static str| self.unnamed.get(&(m, (*c).clone())).map(|s| s.iter().cloned().collect()).unwrap_or_default();
                if self.inside_all(&windows, &seqs, memo)? { inside_u.push((*c).clone()); }
            }
            if inside_u.is_empty() { continue; }
            let mut inside_s = Vec::new();
            for sp in &species {
                let seqs = |m: &'static str| self.items.get(&(m, (*sp).clone())).map(|s| s.iter().map(|(q, _)| q.clone()).collect()).unwrap_or_default();
                if self.inside_all(&windows, &seqs, memo)? { inside_s.push((*sp).clone()); }
            }
            let codes: Vec<String> = inside_u.iter().map(|c| code(c)).collect();
            row.insert("r1a_unnamed_inside".into(), codes.join(", "));
            row.insert("r1a_named_inside".into(), inside_s.join(";"));
            let status = row_value(row, "status")?.to_owned();
            let add = |row: &mut Row, col: &str, text: String| -> Result<()> {
                let existing = row_value(row, col)?.to_owned();
                row.insert(col.into(), if existing.is_empty() { text } else { format!("{existing}; {text}") });
                Ok(())
            };
            let flag = format!("U inside unnamed culture(s) {}", codes.join(", "));
            if inside_s.is_empty() && !matches!(status.as_str(), "divergent class" | "artifact" | "discordant") {
                let genera: BTreeSet<String> = inside_u.iter().map(|c| data.cultures.get(c).map(|m| m.genus.clone()).unwrap_or_default()).collect();
                row.insert("status".into(), "unnamed culture".into());
                row.insert("call_species".into(), String::new());
                row.insert("call_genera".into(), genera.into_iter().collect::<Vec<_>>().join(";"));
                row.insert("placement".into(), inside_u.iter().map(|c| label(c)).collect::<Vec<_>>().join(";"));
                row.insert("novelty_claim".into(), "not_asserted".into());
                add(row, "reasons", format!("U1 inside unnamed culture(s) {} on ITS, LSU and SSU, and inside no named species", codes.join(", ")))?;
                row.insert("r1a_action".into(), "U1".into());
                u1 += 1;
            } else if !inside_s.is_empty() && (status == "known species" || status == "known species group") {
                let mut members: BTreeSet<String> = row_value(row, "call_species")?.split(';').filter(|s| !s.is_empty()).map(str::to_owned).collect();
                members.extend(inside_u.iter().map(|c| label(c)));
                let joined = members.into_iter().collect::<Vec<_>>().join(";");
                row.insert("status".into(), "known species group".into());
                row.insert("call_species".into(), joined.clone());
                row.insert("placement".into(), joined);
                add(row, "flags", flag)?;
                row.insert("r1a_action".into(), "group".into());
            } else {
                add(row, "flags", flag)?;
                row.insert("r1a_action".into(), "flag".into());
            }
        }
        Ok(u1)
    }

    /// Apply the gate to rows named with these training cultures. `query(marker, accession)` gives a
    /// sequence's window. Returns the number of calls demoted.
    pub(crate) fn apply(&self, rows: &mut [Row], query: &dyn Fn(&str, &str) -> Option<String>, memo: &mut DistanceMemo) -> Result<usize> {
        let limits = self.limits_text();
        let mut demoted = 0;
        for row in rows {
            for c in R2_COLUMNS { row.insert(c.into(), String::new()); }
            let status = row_value(row, "status")?.to_owned();
            if status != "known species" && status != "known species group" { continue; }
            let acc = row_value(row, "accession")?.to_owned();
            let members: Vec<String> = row_value(row, "call_species")?.split(';').filter(|s| !s.is_empty()).map(str::to_owned).collect();
            if members.is_empty() { bail!("{acc}: a known call without species"); }
            let windows: Vec<(&'static str, &'static str, Option<String>)> = INSIDE_MARKERS.iter()
                .map(|(m, lab)| (*m, *lab, query(m, &acc).filter(|w| !w.is_empty()))).collect();
            let unchecked: Vec<&str> = windows.iter().filter(|(_, _, w)| w.is_none()).map(|(_, lab, _)| *lab).collect();
            let mut stands = false;
            let mut failing: Vec<String> = Vec::new();
            for sp in &members {
                let mut fails = Vec::new();
                let mut checked = 0;
                for (m, lab, w) in &windows {
                    let (Some(w), Some(t)) = (w, self.limits[m]) else { continue };
                    let Some(d) = self.nearest(m, sp, w, memo)? else { continue };
                    checked += 1;
                    if d > t { fails.push(*lab); }
                }
                if checked >= 2 && fails.is_empty() { stands = true; break; }
                if checked < 2 { fails.push("too few markers"); }
                failing.push(if members.len() > 1 { format!("{sp}: {}", fails.join(",")) } else { fails.join(",") });
            }
            row.insert("r2_inside".into(), if stands { "1" } else { "0" }.into());
            row.insert("r2_unchecked_markers".into(), unchecked.join(","));
            row.insert("r2_limits_pct".into(), limits.clone());
            if stands { continue; }
            row.insert("r2_failing_markers".into(), failing.join("; "));
            let genus = row_value(row, "call_genera")?.to_owned();
            row.insert("status".into(), "novel species, genus known".into());
            row.insert("placement".into(), genus);
            row.insert("novelty_claim".into(), "candidate_only".into());
            let reason = format!("N7 outside the called species ({}; within-species limits {limits})", failing.join("; "));
            let existing = row_value(row, "reasons")?.to_owned();
            row.insert("reasons".into(), if existing.is_empty() { reason } else { format!("{existing}; {reason}") });
            demoted += 1;
        }
        Ok(demoted)
    }
}

pub struct ArmReport {
    /// every table's rows before and after the rule, for scoring
    pub results: Vec<super::scoring::TableResult>,
    pub rows: usize,
    pub known_before: usize,
    pub demoted: usize,
    pub rounds: usize,
    pub limits: Vec<String>,
}

/// R2 on named-row tables (any frozen or published v1/v2 table with `held_cultures`): each row's
/// training set is every reference culture outside its held cultures, so each round gets its own
/// limits. Windows come from the culture reference, then from `extra_windows` (marker -> sequences)
/// for sequences outside it. Writes each table to its output with the R2 columns appended; one
/// distance memo serves every table and round.
pub fn run_r2_on_tables(fixture_dir: &Path, tables: &[(std::path::PathBuf, std::path::PathBuf)], extra_windows: &BTreeMap<String, SeqMap>) -> Result<ArmReport> {
    run_rules_on_tables(fixture_dir, &["r2"], tables, extra_windows)
}

/// The rules (`RULES`, applied in protocol order: R1a, then R2) on named-row tables; see
/// `run_r2_on_tables`. `demoted` counts R2 demotions plus R1a's U1 calls.
pub fn run_rules_on_tables(fixture_dir: &Path, rules: &[&str], tables: &[(std::path::PathBuf, std::path::PathBuf)], extra_windows: &BTreeMap<String, SeqMap>) -> Result<ArmReport> {
    run_rules_with_placements(fixture_dir, rules, tables, extra_windows, None)
}

/// As `run_rules_on_tables`, with the queries' own placement table for sequences outside the culture
/// reference (R7 reads each query's genus clade from it).
pub fn run_rules_with_placements(fixture_dir: &Path, rules: &[&str], tables: &[(std::path::PathBuf, std::path::PathBuf)], extra_windows: &BTreeMap<String, SeqMap>, query_placements: Option<&Path>) -> Result<ArmReport> {
    if rules.is_empty() || rules.iter().any(|r| !RULES.contains(r)) { bail!("rules must be among {}", RULES.join(", ")); }
    if rules.contains(&"r7") && !rules.iter().any(|r| r.starts_with("r3")) { bail!("R7 runs only with an R3 arm (r3a, r3a_all or r3b)"); }
    if rules.contains(&"r8") && rules.contains(&"r8d") { bail!("r8 and r8d are variants of one rule; choose one"); }
    let qplace: BTreeMap<String, Row> = match query_placements { Some(p) => read_row_map(&read_table_any(p)?), None => BTreeMap::new() };
    let data = FrozenData::load(fixture_dir.join("inputs/benchmark"))?;
    let all: StringSet = data.cultures.keys().cloned().collect();
    let query = |m: &str, acc: &str| data.fa.get(m).and_then(|f| f.get(acc)).or_else(|| extra_windows.get(m).and_then(|f| f.get(acc))).cloned();
    let mut memo = DistanceMemo::default();
    let mut fits: BTreeMap<String, InsideRule> = BTreeMap::new();
    let trainers: BTreeMap<String, Trainer> = if rules.contains(&"r1b") {
        ["VT", "SH"].iter().map(|a| Ok(((*a).to_owned(), Trainer::new(&data, a)?))).collect::<Result<_>>()?
    } else { BTreeMap::new() };
    let frozen_placements = fixture_dir.join("placement/lsu_placements.copies.tsv");
    let mut psets_by: BTreeMap<std::path::PathBuf, BTreeMap<String, StringSet>> = BTreeMap::new();
    let mut blocks: BTreeMap<(String, std::path::PathBuf), BlockRule> = BTreeMap::new();
    if rules.contains(&"r5") && rules.contains(&"r5s") { bail!("r5 and r5s are variants of one rule; choose one"); }
    let motif_variant = if rules.contains(&"r5s") { Some(true) } else if rules.contains(&"r5") { Some(false) } else { None };
    let mut motifs: BTreeMap<String, MotifRule> = BTreeMap::new();
    let r3: Vec<&str> = rules.iter().copied().filter(|r| r.starts_with("r3")).collect();
    if r3.len() > 1 { bail!("r3a, r3a_all and r3b are separate arms; choose one"); }
    let mut genus_rules: BTreeMap<String, GenusRule> = BTreeMap::new();
    let (mut n_rows, mut known_before, mut demoted, mut rounds, mut limits) = (0, 0, 0, BTreeSet::new(), BTreeSet::new());
    let mut results = Vec::new();
    for (table, out) in tables {
        let t = read_table_any(table)?;
        if !t.headers.iter().any(|h| h == "held_cultures") { bail!("{} has no held_cultures column", table.display()); }
        let mut by_round: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (i, r) in t.rows.iter().enumerate() { by_round.entry(row_value(r, "held_cultures")?.to_owned()).or_default().push(i); }
        let mut rows = t.rows.clone();
        let frozen_rows = t.rows.clone();
        // LSU units come from the round's own placements when the table sits in an M4 round folder
        let round_placements = table.parent().and_then(|p| p.parent()).map(|r| r.join("05_summary/lsu_placements.copies.tsv")).filter(|p| p.is_file());
        let placements_path = round_placements.unwrap_or_else(|| frozen_placements.clone());
        let place_rows: BTreeMap<String, Row> = if rules.contains(&"r7") { read_row_map(&read_table(&placements_path)?) } else { BTreeMap::new() };
        if rules.contains(&"r1b") && !psets_by.contains_key(&placements_path) {
            let placements = read_row_map(&read_table(&placements_path)?);
            psets_by.insert(placements_path.clone(), species_clade_sets(&data, &placements)?);
        }
        known_before += rows.iter().filter(|r| r.get("status").is_some_and(|s| s.starts_with("known"))).count();
        for (held, idx) in &by_round {
            if !fits.contains_key(held) {
                let held_set: StringSet = held.split(';').filter(|s| !s.is_empty()).map(str::to_owned).collect();
                let training: StringSet = all.difference(&held_set).cloned().collect();
                fits.insert(held.clone(), InsideRule::fit(&data, &training, &mut memo)?);
            }
            let rule = &fits[held];
            limits.insert(rule.limits_text());
            rounds.insert(held.clone());
            let mut part: Vec<Row> = idx.iter().map(|i| rows[*i].clone()).collect();
            if let Some(fallback) = motif_variant {
                if !motifs.contains_key(held) {
                    let held_set: StringSet = held.split(';').filter(|s| !s.is_empty()).map(str::to_owned).collect();
                    let training: StringSet = all.difference(&held_set).cloned().collect();
                    motifs.insert(held.clone(), MotifRule::fit(&data, &training, fallback, &mut memo)?);
                }
                demoted += motifs[held].apply(&data, &mut part, &query, &mut memo)?;
            }
            if rules.contains(&"r1a") { demoted += rule.apply_r1a(&data, &mut part, &query, &mut memo)?; }
            if rules.contains(&"r1b") {
                let key = (held.clone(), placements_path.clone());
                if !blocks.contains_key(&key) {
                    let held_set: StringSet = held.split(';').filter(|s| !s.is_empty()).map(str::to_owned).collect();
                    let training: StringSet = all.difference(&held_set).cloned().collect();
                    blocks.insert(key.clone(), BlockRule::fit(&data, &trainers, &psets_by[&placements_path], &training)?);
                }
                demoted += blocks[&key].apply(&data, &mut part)?;
            }
            if rules.contains(&"r2") { demoted += rule.apply(&mut part, &query, &mut memo)?; }
            if rules.contains(&"r7") {
                let held_set: StringSet = held.split(';').filter(|s| !s.is_empty()).map(str::to_owned).collect();
                let training: StringSet = all.difference(&held_set).cloned().collect();
                let holders = genus_clade_holders(&data, &training, &place_rows);
                demoted += apply_r7(&data, &holders, &mut part, &|acc| qplace.get(acc).or_else(|| place_rows.get(acc)).cloned())?;
            }
            if rules.contains(&"r8") || rules.contains(&"r8d") { demoted += apply_r8(rules.contains(&"r8d"), &mut part)?; }
            if let Some(v) = r3.first() {
                if *v != "r3b" && !genus_rules.contains_key(held) {
                    let held_set: StringSet = held.split(';').filter(|s| !s.is_empty()).map(str::to_owned).collect();
                    let training: StringSet = all.difference(&held_set).cloned().collect();
                    genus_rules.insert(held.clone(), GenusRule::fit(&data, &training, &mut memo)?);
                }
                demoted += apply_r3(v, genus_rules.get(held), &mut part, &query, &mut memo)?;
            }
            for (i, r) in idx.iter().zip(part) { rows[*i] = r; }
        }
        let mut headers = t.headers.clone();
        if motif_variant.is_some() { headers.extend(R5_COLUMNS.iter().map(|c| (*c).to_owned())); }
        if rules.contains(&"r1a") { headers.extend(R1A_COLUMNS.iter().map(|c| (*c).to_owned())); }
        if rules.contains(&"r1b") { headers.extend(R1B_COLUMNS.iter().map(|c| (*c).to_owned())); }
        if rules.contains(&"r2") { headers.extend(R2_COLUMNS.iter().map(|c| (*c).to_owned())); }
        if rules.contains(&"r7") { headers.extend(R7_COLUMNS.iter().map(|c| (*c).to_owned())); }
        if rules.contains(&"r8") || rules.contains(&"r8d") { headers.extend(R8_COLUMNS.iter().map(|c| (*c).to_owned())); }
        if !r3.is_empty() { headers.extend(R3_COLUMNS.iter().map(|c| (*c).to_owned())); }
        std::fs::write(out, serialize_table(&headers, &rows))?;
        n_rows += rows.len();
        results.push(super::scoring::TableResult { table: table.clone(), frozen: frozen_rows, arm: rows });
    }
    Ok(ArmReport { results, rows: n_rows, known_before, demoted, rounds: rounds.len(), limits: limits.into_iter().collect() })
}

/// The three inside-marker windows of sequences outside the culture reference, from a windows prefix
/// (`<prefix>.full_its.fasta`, `<prefix>.LSU_LROR-FLR2.fasta`, `<prefix>.ssu_flank.fasta`).
pub fn read_windows(prefix: &Path) -> Result<BTreeMap<String, SeqMap>> {
    let name = prefix.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    let mut out: BTreeMap<String, SeqMap> = INSIDE_MARKERS.iter()
        .map(|(m, _)| Ok(((*m).to_owned(), read_fasta(&prefix.with_file_name(format!("{name}.{m}.fasta")))?))).collect::<Result<_>>()?;
    let s58 = prefix.with_file_name(format!("{name}.s58.fasta"));
    if s58.is_file() { out.insert("s58".into(), read_fasta(&s58)?); }
    Ok(out)
}

/// Score an arm's tables under one design (see `scoring`), against the culture reference.
pub fn score_results(fixture_dir: &Path, design: &str, results: &[super::scoring::TableResult], euk_predictions: Option<&Path>) -> Result<Vec<super::scoring::OutcomeScore>> {
    let data = FrozenData::load(fixture_dir.join("inputs/benchmark"))?;
    super::scoring::score(&data, design, results, euk_predictions)
}
