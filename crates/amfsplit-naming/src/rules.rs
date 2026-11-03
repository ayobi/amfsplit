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

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

use anyhow::{anyhow, bail, Result};

use super::*;

pub const INSIDE_MARKERS: [(&str, &str); 3] = [("full_its", "ITS"), ("LSU_LROR-FLR2", "LSU"), ("ssu_flank", "SSU")];
pub const R2_COLUMNS: [&str; 4] = ["r2_inside", "r2_failing_markers", "r2_unchecked_markers", "r2_limits_pct"];
pub const R1A_COLUMNS: [&str; 3] = ["r1a_unnamed_inside", "r1a_named_inside", "r1a_action"];
pub const RULES: [&str; 2] = ["r1a", "r2"];

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
    if rules.is_empty() || rules.iter().any(|r| !RULES.contains(r)) { bail!("rules must be among {}", RULES.join(", ")); }
    let data = FrozenData::load(fixture_dir.join("inputs/benchmark"))?;
    let all: StringSet = data.cultures.keys().cloned().collect();
    let query = |m: &str, acc: &str| data.fa.get(m).and_then(|f| f.get(acc)).or_else(|| extra_windows.get(m).and_then(|f| f.get(acc))).cloned();
    let mut memo = DistanceMemo::default();
    let mut fits: BTreeMap<String, InsideRule> = BTreeMap::new();
    let (mut n_rows, mut known_before, mut demoted, mut rounds, mut limits) = (0, 0, 0, BTreeSet::new(), BTreeSet::new());
    let mut results = Vec::new();
    for (table, out) in tables {
        let t = read_table_any(table)?;
        if !t.headers.iter().any(|h| h == "held_cultures") { bail!("{} has no held_cultures column", table.display()); }
        let mut by_round: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for (i, r) in t.rows.iter().enumerate() { by_round.entry(row_value(r, "held_cultures")?.to_owned()).or_default().push(i); }
        let mut rows = t.rows.clone();
        let frozen_rows = t.rows.clone();
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
            if rules.contains(&"r1a") { demoted += rule.apply_r1a(&data, &mut part, &query, &mut memo)?; }
            if rules.contains(&"r2") { demoted += rule.apply(&mut part, &query, &mut memo)?; }
            for (i, r) in idx.iter().zip(part) { rows[*i] = r; }
        }
        let mut headers = t.headers.clone();
        if rules.contains(&"r1a") { headers.extend(R1A_COLUMNS.iter().map(|c| (*c).to_owned())); }
        if rules.contains(&"r2") { headers.extend(R2_COLUMNS.iter().map(|c| (*c).to_owned())); }
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
    INSIDE_MARKERS.iter().map(|(m, _)| Ok(((*m).to_owned(), read_fasta(&prefix.with_file_name(format!("{name}.{m}.fasta")))?))).collect()
}

/// Score an arm's tables under one design (see `scoring`), against the culture reference.
pub fn score_results(fixture_dir: &Path, design: &str, results: &[super::scoring::TableResult], euk_predictions: Option<&Path>) -> Result<Vec<super::scoring::OutcomeScore>> {
    let data = FrozenData::load(fixture_dir.join("inputs/benchmark"))?;
    super::scoring::score(&data, design, results, euk_predictions)
}
