//! Scoring for the Protocol v3 arms: the paper's outcomes, frozen against arm, row by row, with the
//! protocol's paired interval (95% cluster bootstrap resampling species, 10,000 replicates, fixed seed).
//!
//! Designs and their truth:
//! - `culture`: culture holdout. Known species kept: a main-type copy of a named species that another
//!   training culture holds, called a known species or group containing its own species (755).
//!   Exact species: the same copies called exactly their own species.
//! - `species`: species holdout. False known: a main-type copy of a named species called a known
//!   species or group (803).
//! - Both holdouts also report short-5.8S copies recognised: copies of the short-5.8S divergent class
//!   called divergent (43), the paper's measure of recognising a paralog from its own sequence.
//! - `rounds`: the M4 species rounds, first-outcome folds only (no development fold, own placement):
//!   false known among main-type query rows, V18 removed and V18 intact (673 each).
//! - `longglodb`: named test cultures; false known among those whose species v4 lacks (102), kept
//!   among those whose species v4 holds (164).
//! - `eukaryome`: the locked roles of the EUKARYOME test (`predictions.predictions.tsv`): false known
//!   among absent taxa (2,706), kept among known species (199).

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Result};
use serde_json::Value;

use super::*;

pub const BOOTSTRAP_REPLICATES: usize = 10_000;
pub const BOOTSTRAP_SEED: u64 = 20_261_102;
const KNOWN: [&str; 2] = ["known species", "known species group"];

/// One table run through an arm: its rows before and after, in the same order.
pub struct TableResult {
    pub table: PathBuf,
    pub frozen: Vec<Row>,
    pub arm: Vec<Row>,
}

#[derive(Clone, Debug)]
pub struct OutcomeScore {
    pub design: String,
    pub outcome: String,
    pub n: usize,
    pub clusters: usize,
    pub frozen: usize,
    pub arm: usize,
    /// arm minus frozen, in percentage points: 95% cluster-bootstrap interval
    pub diff_low: f64,
    pub diff_high: f64,
}

struct Unit { cluster: String, frozen: bool, arm: bool }

fn known(r: &Row) -> bool { r.get("status").is_some_and(|s| KNOWN.contains(&s.as_str())) }
fn calls(r: &Row, sp: &str) -> bool { r.get("call_species").is_some_and(|c| c.split(';').any(|x| x == sp)) }
fn exact(r: &Row, sp: &str) -> bool { r.get("status").map(String::as_str) == Some("known species") && r.get("call_species").map(String::as_str) == Some(sp) }
fn col<'a>(r: &'a Row, k: &str) -> &'a str { r.get(k).map(String::as_str).unwrap_or("") }

/// SplitMix64: a small, documented generator, so the intervals reproduce anywhere.
struct SplitMix(u64);
impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn below(&mut self, n: usize) -> usize { (self.next() % n as u64) as usize }
}

/// The paired difference (arm minus frozen, in points): resample species with replacement, and take
/// the 2.5th and 97.5th percentiles of the replicate differences.
fn paired_interval(units: &[Unit]) -> (f64, f64, usize) {
    let mut by: BTreeMap<&str, (usize, i64)> = BTreeMap::new();
    for u in units {
        let e = by.entry(u.cluster.as_str()).or_default();
        e.0 += 1;
        e.1 += i64::from(u.arm) - i64::from(u.frozen);
    }
    let clusters: Vec<(usize, i64)> = by.into_values().collect();
    let k = clusters.len();
    if k == 0 { return (0.0, 0.0, 0); }
    let mut rng = SplitMix(BOOTSTRAP_SEED);
    let mut diffs: Vec<f64> = (0..BOOTSTRAP_REPLICATES).map(|_| {
        let (mut n, mut d) = (0usize, 0i64);
        for _ in 0..k { let (cn, cd) = clusters[rng.below(k)]; n += cn; d += cd; }
        if n > 0 { 100.0 * d as f64 / n as f64 } else { 0.0 }
    }).collect();
    diffs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let lo = diffs[(0.025 * BOOTSTRAP_REPLICATES as f64) as usize];
    let hi = diffs[(0.975 * BOOTSTRAP_REPLICATES as f64) as usize - 1];
    (lo, hi, k)
}

fn summarize(design: &str, outcome: &str, units: Vec<Unit>) -> OutcomeScore {
    let (diff_low, diff_high, clusters) = paired_interval(&units);
    OutcomeScore {
        design: design.into(), outcome: outcome.into(), n: units.len(), clusters,
        frozen: units.iter().filter(|u| u.frozen).count(), arm: units.iter().filter(|u| u.arm).count(), diff_low, diff_high,
    }
}

fn pairs(t: &TableResult) -> Result<impl Iterator<Item = (&Row, &Row)>> {
    if t.frozen.len() != t.arm.len() { bail!("{}: frozen and arm rows differ in number", t.table.display()); }
    for (f, a) in t.frozen.iter().zip(&t.arm) {
        if f.get("accession") != a.get("accession") { bail!("{}: frozen and arm rows are not paired", t.table.display()); }
    }
    Ok(t.frozen.iter().zip(&t.arm))
}

/// Short-5.8S divergent copies called divergent, clustered by species (or culture, when unnamed).
fn short_recognised(data: &FrozenData, results: &[TableResult]) -> Result<Vec<Unit>> {
    let mut units = Vec::new();
    for t in results {
        for (f, a) in pairs(t)? {
            if col(f, "copy_class") != "divergent" || col(f, "div_kind") != "short_5.8S" { continue; }
            let c = col(f, "culture");
            let cluster = data.species.get(c).cloned().unwrap_or_else(|| c.to_owned());
            units.push(Unit { cluster, frozen: col(f, "status") == "divergent class", arm: col(a, "status") == "divergent class" });
        }
    }
    Ok(units)
}

/// Score an arm's tables under one design. `euk_predictions` is the EUKARYOME test's locked
/// `predictions.predictions.tsv` (design `eukaryome` only).
pub fn score(data: &FrozenData, design: &str, results: &[TableResult], euk_predictions: Option<&Path>) -> Result<Vec<OutcomeScore>> {
    let v4_species: StringSet = data.species.values().cloned().collect();
    let mut out = Vec::new();
    match design {
        "culture" => {
            let (mut kept, mut ex) = (Vec::new(), Vec::new());
            for t in results {
                for (f, a) in pairs(t)? {
                    let Some(sp) = data.species.get(col(f, "culture")) else { continue };
                    if col(f, "copy_class") != "main" { continue; }
                    let held: StringSet = col(f, "held_cultures").split(';').map(str::to_owned).collect();
                    let other = data.by_species.get(sp).is_some_and(|cs| cs.iter().any(|c| !held.contains(c)));
                    if !other { continue; }
                    kept.push(Unit { cluster: sp.clone(), frozen: known(f) && calls(f, sp), arm: known(a) && calls(a, sp) });
                    ex.push(Unit { cluster: sp.clone(), frozen: exact(f, sp), arm: exact(a, sp) });
                }
            }
            out.push(summarize("culture holdout", "known species kept", kept));
            out.push(summarize("culture holdout", "exact species calls", ex));
            out.push(summarize("culture holdout", "short-5.8S copies recognised", short_recognised(data, results)?));
        }
        "species" => {
            let mut fk = Vec::new();
            for t in results {
                for (f, a) in pairs(t)? {
                    let Some(sp) = data.species.get(col(f, "culture")) else { continue };
                    if col(f, "copy_class") != "main" { continue; }
                    fk.push(Unit { cluster: sp.clone(), frozen: known(f), arm: known(a) });
                }
            }
            out.push(summarize("species holdout", "false known", fk));
            out.push(summarize("species holdout", "short-5.8S copies recognised", short_recognised(data, results)?));
        }
        "rounds" => {
            let mut by_kind: BTreeMap<String, Vec<Unit>> = BTreeMap::new();
            let mut manifests: BTreeMap<PathBuf, Value> = BTreeMap::new();
            for t in results {
                // rounds/<key>/06_evaluation_v2/<kind>.copies.tsv, with metrics.json beside it
                let dir = t.table.parent().ok_or_else(|| anyhow!("{}: not in a round folder", t.table.display()))?;
                let round = dir.parent().ok_or_else(|| anyhow!("{}: not in a round folder", t.table.display()))?;
                let key = round.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                let mpath = round.parent().map(|p| p.join("MANIFEST.json")).unwrap_or_default();
                if !manifests.contains_key(&mpath) { manifests.insert(mpath.clone(), serde_json::from_str(&fs::read_to_string(&mpath)?)?); }
                let own = manifests[&mpath].pointer(&format!("/rounds/{key}/own_placement")).and_then(Value::as_bool).unwrap_or(false);
                let metrics: Value = serde_json::from_str(&fs::read_to_string(dir.join("metrics.json"))?)?;
                let dev = metrics.get("development_fold").and_then(Value::as_bool).unwrap_or(false);
                if !key.starts_with("species_") || dev || !own { continue; }
                let file = t.table.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                let kind = if file.starts_with("external_holdout") { "false known, V18 removed" } else if file.starts_with("full_reference") { "false known, V18 intact" } else { bail!("{}: unknown round table", t.table.display()) };
                for (f, a) in pairs(t)? {
                    if col(f, "copy_class") != "main" { continue; }
                    by_kind.entry(kind.into()).or_default().push(Unit { cluster: key.clone(), frozen: known(f), arm: known(a) });
                }
            }
            for (kind, units) in by_kind { out.push(summarize("LSU-removal rounds", &kind, units)); }
        }
        "longglodb" => {
            let (mut fk, mut kept) = (Vec::new(), Vec::new());
            for t in results {
                for (f, a) in pairs(t)? {
                    let o = col(f, "organism");
                    if !is_named(o) { continue; }
                    if v4_species.contains(o) { kept.push(Unit { cluster: o.into(), frozen: known(f) && calls(f, o), arm: known(a) && calls(a, o) }); }
                    else { fk.push(Unit { cluster: o.into(), frozen: known(f), arm: known(a) }); }
                }
            }
            out.push(summarize("LongGloDB", "false known (absent taxa)", fk));
            out.push(summarize("LongGloDB", "known species kept", kept));
        }
        "eukaryome" => {
            let path = euk_predictions.ok_or_else(|| anyhow!("design eukaryome needs the test's predictions.predictions.tsv"))?;
            let roles: BTreeMap<String, Row> = read_table(path)?.rows.into_iter().map(|r| Ok((row_value(&r, "asv")?.to_owned(), r))).collect::<Result<_>>()?;
            let (mut fk, mut kept) = (Vec::new(), Vec::new());
            for t in results {
                for (f, a) in pairs(t)? {
                    let x = col(f, "accession");
                    let p = roles.get(x).ok_or_else(|| anyhow!("{x} has no EUKARYOME role"))?;
                    let own = col(p, "own_species_in_v4");
                    let cluster = col(p, "organism").to_owned();
                    if !own.is_empty() { kept.push(Unit { cluster, frozen: known(f) && calls(f, own), arm: known(a) && calls(a, own) }); }
                    else if !col(p, "role").starts_with("genus") { fk.push(Unit { cluster, frozen: known(f), arm: known(a) }); }
                }
            }
            out.push(summarize("EUKARYOME", "false known (absent taxa)", fk));
            out.push(summarize("EUKARYOME", "known species kept", kept));
        }
        _ => bail!("unknown scoring design {design} (culture, species, rounds, longglodb, eukaryome)"),
    }
    Ok(out)
}

pub const SCORE_HEADERS: [&str; 8] = ["design", "outcome", "n", "species_clusters", "frozen", "arm", "difference_low_pts", "difference_high_pts"];

pub fn scores_tsv(scores: &[OutcomeScore]) -> Vec<u8> {
    let headers: Vec<String> = SCORE_HEADERS.iter().map(|h| (*h).to_owned()).collect();
    let rows: Vec<Row> = scores.iter().map(|s| {
        let vals = [s.design.clone(), s.outcome.clone(), s.n.to_string(), s.clusters.to_string(), s.frozen.to_string(), s.arm.to_string(),
                    format!("{:.2}", s.diff_low), format!("{:.2}", s.diff_high)];
        SCORE_HEADERS.iter().map(|h| (*h).to_owned()).zip(vals).collect()
    }).collect();
    serialize_table(&headers, &rows)
}
