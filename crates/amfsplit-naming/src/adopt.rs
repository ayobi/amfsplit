//! Protocol v3, section 8, applied mechanically (`amfsplit adopt`): each rule run alone on the strata (F plus
//! that rule), its adoption rule applied to the scores, S checked last with the adopted rules only, and the
//! report written. For the dress rehearsal the inputs are the development stand-ins: LongGloDB v1 and
//! EUKARYOME v2.1 (conditions A and W), the 15 genus rounds, R6's negatives and the secondary VT analyses.
use super::rules::{read_windows, run_rules_multi, run_rules_on_tables, score_results, R4Inputs};
use super::scoring::{summarize, OutcomeScore, Unit};
use super::{read_table_any, Row};
use anyhow::{anyhow, Result};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

pub struct AdoptInputs {
    pub benchmark: PathBuf,
    /// local LongGloDB windows folder (lgdb_test.*, test_evidence.json.gz)
    pub lgdb: PathBuf,
    /// local EUKARYOME windows folder (overlap.win.*)
    pub euk: PathBuf,
    pub out: PathBuf,
    /// run only these rules (r1, r2, r4, r5, r3, r7, r8, r6, v, s); empty = all
    pub only: Vec<String>,
}

pub struct Decision { pub rule: String, pub criterion: String, pub measured: String, pub verdict: String }

fn pct(x: usize, n: usize) -> f64 { if n == 0 { 0.0 } else { 100.0 * x as f64 / n as f64 } }
fn line(o: &OutcomeScore) -> String { format!("{} of {} (frozen {}; difference {:+.2} to {:+.2} points)", o.arm, o.n, o.frozen, o.diff_low, o.diff_high) }

struct Ctx<'a> {
    fx: &'a Path, tables: Vec<PathBuf>, windows: BTreeMap<String, BTreeMap<String, String>>, qps: Vec<PathBuf>,
    triage_pos: PathBuf, r4: R4Inputs, genus: Vec<PathBuf>, out: &'a Path, scores: Vec<OutcomeScore>,
}

impl Ctx<'_> {
    fn strata(&mut self, arm: &str, chain: &[&str]) -> Result<BTreeMap<String, OutcomeScore>> {
        let tmp = self.out.join("work").join(arm);
        fs::create_dir_all(&tmp)?;
        let pairs: Vec<_> = self.tables.iter().enumerate().map(|(i, p)| (p.clone(), tmp.join(format!("{i}.tsv")))).collect();
        let qps: Vec<&Path> = self.qps.iter().map(|p| p.as_path()).collect();
        let r4 = chain.contains(&"r4").then_some(&self.r4);
        let triage = chain.contains(&"r6").then_some(self.triage_pos.as_path());
        let r = run_rules_multi(self.fx, chain, &pairs, &self.windows, &qps, r4, triage)?;
        let mut out = BTreeMap::new();
        for mut s in score_results(self.fx, "strata", &r.results, None)? { s.design = arm.to_owned(); out.insert(s.outcome.clone(), s.clone()); self.scores.push(s); }
        Ok(out)
    }
    fn genus_rounds(&mut self, arm: &str, chain: &[&str]) -> Result<BTreeMap<String, OutcomeScore>> {
        let tmp = self.out.join("work").join(format!("{arm}_genus"));
        fs::create_dir_all(&tmp)?;
        let pairs: Vec<_> = self.genus.iter().map(|p| {
            let key = p.parent().and_then(|d| d.parent()).and_then(|d| d.file_name()).map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
            (p.clone(), tmp.join(format!("{key}__{}", p.file_name().unwrap().to_string_lossy())))
        }).collect();
        let r = run_rules_on_tables(self.fx, chain, &pairs, &BTreeMap::new())?;
        let mut out = BTreeMap::new();
        for mut s in score_results(self.fx, "genus", &r.results, None)? { s.design = format!("{arm} (genus rounds)"); out.insert(s.outcome.clone(), s.clone()); self.scores.push(s); }
        Ok(out)
    }
}

fn get<'a>(m: &'a BTreeMap<String, OutcomeScore>, k: &str) -> Result<&'a OutcomeScore> { m.get(k).ok_or_else(|| anyhow!("no outcome {k}")) }

pub fn run_adopt(fixture_dir: &Path, inp: &AdoptInputs) -> Result<Vec<Decision>> {
    let b = inp.benchmark.join("benchmarks");
    let (m4, p3) = (b.join("m4_external_holdout"), b.join("protocol_v3"));
    let wk = p3.join("withheld_k");
    let mut windows: BTreeMap<String, BTreeMap<String, String>> = BTreeMap::new();
    for w in [inp.lgdb.join("lgdb_test"), inp.euk.join("overlap.win")] { for (m, s) in read_windows(&w)? { windows.entry(m).or_default().extend(s); } }
    let mut genus = Vec::new();
    for e in fs::read_dir(m4.join("rounds"))? {
        let d = e?.path();
        if d.file_name().is_some_and(|n| n.to_string_lossy().starts_with("genus_")) {
            for k in ["external_holdout_v2", "full_reference_v2"] { genus.push(d.join(format!("06_evaluation_v2/{k}.copies.tsv"))); }
        }
    }
    genus.sort();
    let r4 = R4Inputs::prepare(fixture_dir, "pooled", 97.5, &read_windows(&inp.lgdb.join("lgdb_test"))?,
                               &m4.join("runs_lgdb/06_scores/A_primary_v1.copies.tsv"), &p3.join("r4_name_checks.tsv"))?;
    let mut c = Ctx {
        fx: fixture_dir, windows, r4, genus, out: &inp.out, scores: Vec::new(), triage_pos: p3.join("r6_triage_positives.tsv"),
        tables: vec![m4.join("runs_lgdb/06_scores/A_primary_v2.copies.tsv"), m4.join("runs_euk/06_scores/A_primary_v2.copies.tsv.gz"),
                     wk.join("W_lgdb/W_primary_v2.copies.tsv.gz"), wk.join("W_euk/W_primary_v2.copies.tsv.gz")],
        qps: vec![m4.join("runs_lgdb/A/05_summary/lsu_placements.copies.tsv"), m4.join("runs_euk/A/05_summary/lsu_placements.copies.tsv.gz")],
    };
    let run = |r: &str| inp.only.is_empty() || inp.only.iter().any(|o| o == r);
    let mut d = Vec::new();
    let mut adopted: BTreeMap<&str, String> = BTreeMap::new();   // rule -> variant
    const KEPT: &str = "K, reference intact: known species kept";
    const FK: &str = "A, L and withheld K: false known";
    // R1a, R1b, R2: added demotion of K bounded; at least one targeted error removed
    let mut r1_pass = Vec::new();
    for (rule, ceiling) in [("r1a", 1.0), ("r1b", 1.0), ("r2", 3.0)] {
        if !run(&rule[..2]) { continue; }
        let s = c.strata(rule, &[rule])?;
        let (k, f) = (get(&s, KEPT)?, get(&s, FK)?);
        let bound = -k.diff_low;
        let benefit = f.arm < f.frozen;
        let pass = bound <= ceiling && benefit;
        if pass && rule != "r2" { r1_pass.push(rule); }
        if pass && rule == "r2" { adopted.insert("r2", "r2".into()); }
        d.push(Decision { rule: rule.to_uppercase(), criterion: format!("added demotion of K, reference intact, bounded at {ceiling} point(s); removes at least one targeted error"),
            measured: format!("K kept {}; demotion bound {bound:.2} points; false known in A, L and withheld K {} to {}", line(k), f.frozen, f.arm),
            verdict: if pass { "passes".into() } else { "fails".into() } });
    }
    if let Some(v) = r1_pass.first() { adopted.insert("r1", (*v).to_owned()); }
    // R4: demotion of K falls; false known in A, L and withheld K rises at most 1 point, bound under 3
    if run("r4") {
        let s = c.strata("r4", &["r4"])?;
        let (k, f) = (get(&s, KEPT)?, get(&s, FK)?);
        let rise = pct(f.arm, f.n) - pct(f.frozen, f.n);
        let pass = k.arm > k.frozen && rise <= 1.0 && f.diff_high < 3.0;
        if pass { adopted.insert("r4", "r4".into()); }
        d.push(Decision { rule: "R4".into(), criterion: "demotion of K falls; false known in A, L and withheld K rises at most 1 point, upper bound under 3".into(),
            measured: format!("K kept {}; false known {}", line(k), line(f)), verdict: if pass { "passes".into() } else { "fails".into() } });
    }
    // R5: recall at least 80%, other carriers reclassified at most 1%; the higher-recall variant meeting the cap
    if run("r5") {
        let mut best: Option<(&str, f64)> = None;
        for v in ["r5", "r5s"] {
            let s = c.strata(v, &[v])?;
            let (rh, ot) = (get(&s, "Rhizoglomus-labelled motif carriers called divergent (R5)")?, get(&s, "other motif carriers reclassified (R5)")?);
            let (recall, reclass) = (pct(rh.arm, rh.n), pct(ot.arm, ot.n));
            let pass = recall >= 80.0 && reclass <= 1.0;
            if pass && best.is_none_or(|(_, r)| recall > r) { best = Some((v, recall)); }
            d.push(Decision { rule: v.to_uppercase(), criterion: "at least 80% of Rhizoglomus-labelled carriers called divergent; at most 1% of other carriers reclassified".into(),
                measured: format!("{recall:.1}% ({} of {}); {reclass:.2}% ({} of {})", rh.arm, rh.n, ot.arm, ot.n), verdict: if pass { "passes".into() } else { "fails".into() } });
        }
        if let Some((v, _)) = best { adopted.insert("r5", v.to_owned()); }
    }
    // R3: at least 50% of wrong genera caught (genus rounds and strata A and L), at most 2% of right genera flagged
    if run("r3") {
        let mut best: Option<(&str, f64)> = None;
        for v in ["r3a", "r3a_all", "r3b"] {
            let s = c.strata(v, &[v])?;
            let g = c.genus_rounds(v, &[v])?;
            let (mut caught, mut wrong, mut flagged, mut right) = (0, 0, 0, 0);
            for (k, o) in g.iter().chain(s.iter()) {
                if k.starts_with("wrong-genus calls flagged") || k.starts_with("A and L: wrong genera flagged") { caught += o.arm; wrong += o.n; }
                if k.starts_with("right-genus calls falsely flagged") || k.starts_with("K and A: right genera falsely flagged") { flagged += o.arm; right += o.n; }
            }
            let (recall, false_rate) = (pct(caught, wrong), pct(flagged, right));
            let pass = recall >= 50.0 && false_rate <= 2.0;
            if pass && best.is_none_or(|(_, r)| recall > r) { best = Some((v, recall)); }
            d.push(Decision { rule: v.to_uppercase(), criterion: "at least 50% of wrong genera flagged in the genus rounds and strata A and L; at most 2% of right genera flagged".into(),
                measured: format!("{recall:.1}% ({caught} of {wrong}); {false_rate:.2}% ({flagged} of {right})"), verdict: if pass { "passes".into() } else { "fails".into() } });
        }
        if let Some((v, _)) = best { adopted.insert("r3", v.to_owned()); }
    }
    // R7: with the chosen R3 only
    if run("r7") {
        let r3 = adopted.get("r3").cloned();
        let with = r3.clone().unwrap_or_else(|| "r3b".into());
        let s = c.strata("r7", &["r7", with.as_str()])?;
        let (right, flag) = (get(&s, "A: new lineages given the right genus (R7)")?, get(&s, "wrong genera added by R7: flagged")?);
        let pass = r3.is_some() && pct(right.arm, right.n) >= 50.0 && (flag.n == 0 || pct(flag.arm, flag.n) >= 80.0);
        if pass { adopted.insert("r7", "r7".into()); }
        d.push(Decision { rule: "R7".into(), criterion: "with an adopted R3: the right genus for at least 50% of F's genus-less calls whose genus another culture holds; at least 80% of its wrong genera flagged".into(),
            measured: format!("right genus {} of {}; wrong genera flagged {} of {} (with {with}{})", right.arm, right.n, flag.arm, flag.n, if r3.is_none() { ", no R3 adopted" } else { "" }),
            verdict: if pass { "passes".into() } else if r3.is_none() { "stays opt-in (no R3 adopted)".into() } else { "fails".into() } });
    }
    // R8: the flag is on by default; R8d only if single-marker new lineages turn up among K, reference intact
    if run("r8") {
        let s = c.strata("r8", &["r8"])?;
        let k = get(&s, "K, reference intact: single-marker new lineages (R8)")?;
        adopted.insert("r8", if k.arm > 0 { "r8d".into() } else { "r8".into() });
        d.push(Decision { rule: "R8".into(), criterion: "flag on by default; R8d adopted only if single-marker new lineages turn up among K with the reference intact".into(),
            measured: format!("{} of {} K sequences", k.arm, k.n), verdict: if k.arm > 0 { "R8d adopted".into() } else { "flag on; R8d not adopted".into() } });
    }
    // R6: at most 1% of triage negatives called new lineage, upper bound under 2%
    if run("r6") {
        let named = read_table_any(&p3.join("r6_negatives_named_rows.tsv.gz"))?.rows;
        let tmp = inp.out.join("work/r6_negatives");
        fs::create_dir_all(&tmp)?;
        let r = super::rules::run_rules_triage(fixture_dir, &["r6"], &[(p3.join("r6_negatives_named_rows.tsv.gz"), tmp.join("named.tsv"))],
                                              &BTreeMap::new(), None, None, Some(&p3.join("r6_negatives_triage.tsv.gz")))?;
        let after: BTreeMap<String, String> = r.results.iter().flat_map(|t| t.arm.iter()).map(|row| (row["accession"].clone(), row["status"].clone())).collect();
        let genus = |id: &str| id.split(';').find_map(|t| t.strip_prefix("g__")).unwrap_or("unknown").to_owned();
        let labels = read_table_any(&p3.join("r6_negatives_triage.tsv.gz"))?.rows;
        // the negatives with all three regions at 300 nt or more: those triaged, plus those without any hit
        let total = NEGATIVES_WITH_THREE_REGIONS;
        let mut units: Vec<Unit> = labels.iter().map(|l| {
            let id = l["read_id"].clone();
            let nl = after.get(&id).is_some_and(|s| s.starts_with("new lineage") || s == "novel lineage");
            Unit { cluster: genus(&id), frozen: false, arm: nl }
        }).collect();
        while units.len() < total { units.push(Unit { cluster: "no hit".into(), frozen: false, arm: false }); }
        let o = summarize("r6 (negatives)", "triage negatives called new lineage", units);
        let rate = pct(o.arm, o.n);
        let pass = rate <= 1.0 && o.diff_high < 2.0;
        if pass { adopted.insert("r6", "r6".into()); }
        d.push(Decision { rule: "R6".into(), criterion: "at most 1% of triage negatives called new lineage, upper bound under 2%".into(),
            measured: format!("{rate:.2}% ({} of {}; upper bound {:.2}%; {} named under R6)", o.arm, o.n, o.diff_high, named.len()),
            verdict: if pass { "labelling replaces filtering".into() } else { "filtering stays (each run reports what triage dropped)".into() } });
        c.scores.push(o);
    }
    // V: with an adopted R3, every wrong genus given through a VT alone flagged; no false known through a VT shared by species
    if run("v") {
        let r3 = adopted.get("r3").cloned();
        let with = r3.clone().unwrap_or_else(|| "r3b".into());
        let (mut vt_wrong, mut vt_flagged, mut shared_fk) = (0, 0, 0);
        for (sec, prim, w) in [(m4.join("runs_lgdb/06_scores/A_secondary_vt_v2.copies.tsv"), m4.join("runs_lgdb/06_scores/A_primary_v2.copies.tsv"), &c.windows),
                               (m4.join("runs_euk/06_scores/A_secondary_vt_v2.copies.tsv.gz"), m4.join("runs_euk/06_scores/A_primary_v2.copies.tsv.gz"), &c.windows)] {
            let tmp = inp.out.join("work/v");
            fs::create_dir_all(&tmp)?;
            let qps: Vec<&Path> = c.qps.iter().map(|p| p.as_path()).collect();
            let r = run_rules_multi(fixture_dir, &[with.as_str()], &[(sec.clone(), tmp.join("sec.tsv"))], w, &qps, None, None)?;
            let primary: BTreeMap<String, Row> = read_table_any(&prim)?.rows.into_iter().map(|row| (row["accession"].clone(), row)).collect();
            for row in r.results.iter().flat_map(|t| t.arm.iter()) {
                let og = super::normalise_genus(row["organism"].split_whitespace().next().unwrap_or(""));
                let genera: Vec<&str> = row["call_genera"].split(';').filter(|g| !g.is_empty()).collect();
                let vt_alone = row["VT_state"] == "matched" && row["SH_state"] != "matched" && row["LSU_state"] != "matched";
                if vt_alone && !genera.is_empty() && !genera.contains(&og.as_str()) && row["test_role"] != "genus only" {
                    vt_wrong += 1;
                    if row.get("r3_flag").map(String::as_str) == Some("G1") { vt_flagged += 1; }
                }
                let known = |r: &Row| matches!(r.get("status").map(String::as_str), Some("known species" | "known species group"));
                let absent = row["test_role"] == "novel species, genus in v4" || row["test_role"] == "novel lineage, genus not in v4";
                if absent && known(row) && !primary.get(&row["accession"]).is_some_and(known)
                    && row["VT_species"].split(';').filter(|s| !s.is_empty()).count() > 1 { shared_fk += 1; }
            }
        }
        let pass = r3.is_some() && vt_wrong == vt_flagged && shared_fk == 0;
        if pass { adopted.insert("v", "v".into()); }
        d.push(Decision { rule: "V".into(), criterion: "with an adopted R3: every wrong genus given through a VT alone flagged; no false known through a VT shared by species".into(),
            measured: format!("VT-alone wrong genera flagged {vt_flagged} of {vt_wrong}; new false knowns through a shared VT {shared_fk} (with {with})"),
            verdict: if pass { "passes".into() } else if r3.is_none() { "stays opt-in (no R3 adopted)".into() } else { "fails".into() } });
    }
    // S: the adopted rules together; demotion within the adopted rules' largest ceiling
    if run("s") {
        let order = ["r1", "r2", "r4", "r5", "r7"];
        let ceiling = |a: &BTreeMap<&str, String>| -> Option<f64> {
            let mut c: Option<f64> = None;
            if a.contains_key("r1") { c = Some(c.unwrap_or(0.0).max(1.0)); }
            if a.contains_key("r2") { c = Some(c.unwrap_or(0.0).max(3.0)); }
            c
        };
        let mut set = adopted.clone();
        let mut dropped = Vec::new();
        loop {
            let mut chain: Vec<&str> = Vec::new();
            for (k, v) in [("r5", set.get("r5")), ("r4", set.get("r4")), ("r1", set.get("r1")), ("r2", set.get("r2")), ("r6", set.get("r6")), ("r7", set.get("r7")), ("r8", set.get("r8")), ("r3", set.get("r3"))] {
                if let Some(v) = v { chain.push(if k == "r1" || k == "r5" || k == "r8" || k == "r3" { v.as_str() } else { k }); }
            }
            if chain.is_empty() { d.push(Decision { rule: "S".into(), criterion: "the adopted rules together".into(), measured: "no rule adopted".into(), verdict: "nothing to check".into() }); break; }
            let s = c.strata("s", &chain)?;
            let k = get(&s, KEPT)?;
            let bound = -k.diff_low;
            let within = ceiling(&set).is_none_or(|cl| bound <= cl);
            if within || !order.iter().any(|r| set.contains_key(r)) {
                d.push(Decision { rule: "S".into(), criterion: "demotion of K within the adopted rules' largest ceiling".into(),
                    measured: format!("{}: K kept {}; demotion bound {bound:.2} points{}", chain.join(", "), line(k),
                                      if dropped.is_empty() { String::new() } else { format!("; left opt-in: {}", dropped.join(", ")) }),
                    verdict: if within { "passes".into() } else { "fails with no rule left to drop".into() } });
                break;
            }
            let last = *order.iter().rev().find(|r| set.contains_key(*r)).unwrap();
            set.remove(last);
            dropped.push(last.to_uppercase());
        }
    }
    // the report
    fs::create_dir_all(&inp.out)?;
    let mut tsv = String::from("design\toutcome\tn\tclusters\tfrozen\tarm\tdiff_low\tdiff_high\n");
    for s in &c.scores { tsv.push_str(&format!("{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\n", s.design, s.outcome, s.n, s.clusters, s.frozen, s.arm, s.diff_low, s.diff_high)); }
    fs::write(inp.out.join("adopt_scores.tsv"), tsv)?;
    let mut md = String::from("# Protocol v3, section 8 applied\n\nDevelopment stand-ins (LongGloDB v1, EUKARYOME v2.1): a rehearsal, which decides nothing. The power check runs at the set lock.\n\n| Rule | Criterion | Measured | Verdict |\n| --- | --- | --- | --- |\n");
    for x in &d { md.push_str(&format!("| {} | {} | {} | {} |\n", x.rule, x.criterion, x.measured, x.verdict)); }
    fs::write(inp.out.join("adopt_report.md"), md)?;
    Ok(d)
}

/// Non-Glomeromycota EUKARYOME v2.1 reads with full ITS, LR0R-FLR2 and SSU flank windows of at least 300 nt
/// (protocol_v3/README.md): 19,263 with a validation hit, in the triage table, and 139 without one.
const NEGATIVES_WITH_THREE_REGIONS: usize = 19_402;
