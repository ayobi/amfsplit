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

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

use anyhow::{anyhow, bail, Result};

use super::*;

pub const INSIDE_MARKERS: [(&str, &str); 3] = [("full_its", "ITS"), ("LSU_LROR-FLR2", "LSU"), ("ssu_flank", "SSU")];
pub const R2_COLUMNS: [&str; 4] = ["r2_inside", "r2_failing_markers", "r2_unchecked_markers", "r2_limits_pct"];

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
    /// per marker: the limit, or None when fewer than two values exist
    pub limits: BTreeMap<&'static str, Option<f64>>,
}

impl InsideRule {
    /// species_groups.py: the 99th percentile of the nearest other-culture distance within a species.
    pub(crate) fn fit(data: &FrozenData, training: &StringSet, memo: &mut DistanceMemo) -> Result<Self> {
        let items = species_items(data, training);
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
        Ok(Self { items, limits })
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
    let data = FrozenData::load(fixture_dir.join("inputs/benchmark"))?;
    let all: StringSet = data.cultures.keys().cloned().collect();
    let query = |m: &str, acc: &str| data.fa.get(m).and_then(|f| f.get(acc)).or_else(|| extra_windows.get(m).and_then(|f| f.get(acc))).cloned();
    let mut memo = DistanceMemo::default();
    let mut fits: BTreeMap<String, InsideRule> = BTreeMap::new();
    let (mut n_rows, mut known_before, mut demoted, mut rounds, mut limits) = (0, 0, 0, BTreeSet::new(), BTreeSet::new());
    for (table, out) in tables {
        let t = read_table_any(table)?;
        if !t.headers.iter().any(|h| h == "held_cultures") { bail!("{} has no held_cultures column", table.display()); }
        let mut by_round: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (i, r) in t.rows.iter().enumerate() { by_round.entry(row_value(r, "held_cultures")?.to_owned()).or_default().push(i); }
        let mut rows = t.rows.clone();
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
            demoted += rule.apply(&mut part, &query, &mut memo)?;
            for (i, r) in idx.iter().zip(part) { rows[*i] = r; }
        }
        let mut headers = t.headers.clone();
        headers.extend(R2_COLUMNS.iter().map(|c| (*c).to_owned()));
        std::fs::write(out, serialize_table(&headers, &rows))?;
        n_rows += rows.len();
    }
    Ok(ArmReport { rows: n_rows, known_before, demoted, rounds: rounds.len(), limits: limits.into_iter().collect() })
}
