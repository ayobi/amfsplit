//! Gate 0b, LongGloDB: the independent-culture test of amf-novelty-benchmark (`lgdb_test.py`,
//! PROTOCOL_LGDB.md with amendments 1 and 2) reproduced in Rust.
//!
//! The 280 test ASVs join the frozen data as extra copies of extra cultures, every test culture is
//! held out of every call, and each ASV is named with the species-clade arm and then the ITS distance
//! rule (N6). Condition A places the test ASVs on the full V18 reference; each of the 19 condition B
//! rounds re-places everything with one novel species' V18 records removed. Both are scored under the
//! primary VT evidence and the secondary one (amendment 2: coverage over the spanned part of the
//! reference), before and after N6: eight tables, compared byte for byte with `runs_lgdb/06_scores/`.
//!
//! The Stage 1 windows of the test ASVs are LongGloDB sequences and are never redistributed, so this
//! gate reads them, with the test evidence, from a local folder checked against the hashes the
//! published run recorded. Everything else comes from an amf-novelty-benchmark checkout.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};
use serde_json::Value;

use super::*;

const TEST_ROLES: [&str; 4] = ["known species", "novel species, genus in v4", "novel lineage, genus not in v4", "genus only"];
const WINDOWS: [(&str, &str); 4] = [
    ("ssu_flank", "d76e02fd1ee62d649de794e0519c8fae15c5d39229c441b4b7fc3d0538be4d9e"),
    ("full_its", "bc1ae16bf4d7491e6c0be7f0c138f1200ca38282dab726a2bcb9f7c8be6d6465"),
    ("s58", "46da4e66960256f86066a7d7670c4316a3d78b380395dabc25bb05dfb20886f7"),
    ("LSU_LROR-FLR2", "82b147d60e796c4d172702c9f3f7979f0b660a5015bca8d827e7bd151f9ca190"),
];
const EVIDENCE_SHA: &str = "543f61883c4116ea1107495da67588bc47725900c0d5c6889e183762a02cf5fe";
const OVERLAP_SHA: &str = "e39919cc5118abd0a24ed738028b8032fdc4aba497be78bdc6deb2ebe8096501";
const PLAN_SHA: &str = "821b38074cf045bb9eb099b33a48215c16cc98c482d396febc1d372bf42e44cb";
const A_PLACEMENT_SHA: &str = "6eae1912f400af0954c09beb421ebdc561157fef657d02511c61e90d7e889805";
/// SHA-256 of the sorted lines "<key>\t<SHA-256 of the round's placement table>\n" of the 19 B rounds.
const B_PLACEMENTS_DIGEST: &str = "1b6bddb263ddca0e9f7b5a5ca17c693219eb80f31bc2698609b0e41db4fe5b5e";
const TABLES: [(&str, &str); 8] = [
    ("A_primary_v1", "32e8c0dc5bd88fe89f31e4fb39bcaae4f71c5cc55a6585a29b219953a231c047"),
    ("A_primary_v2", "f9fab9e8ec2a64042b96db3e6a387946a0b6d6b04abdde6892d3b3676149b4a0"),
    ("A_secondary_vt_v1", "83732ef1f5126714e446a041a7487a9b8c8c1ed2d95773be91eab344df76c4c0"),
    ("A_secondary_vt_v2", "2010f464b93cd8a67ef783e6ea3d136f40288a3cf43e67a6036d598512cac26e"),
    ("B_primary_v1", "c44d6bebcb19b4a3942f424d960a856af433905f9e3373328dfa2b9d88b8d3fa"),
    ("B_primary_v2", "b6b1dec756b3d6c99b89c77bab8ce60c8ffe1f728bddc02f9269d1408ae3225c"),
    ("B_secondary_vt_v1", "b84248ca32c8564248597b5e57311a4cc1dfeec0237058e1883792de491297cb"),
    ("B_secondary_vt_v2", "ced1c3cf61ace3c3fb3f685a86b483d848d400e427e318b9614452f5d08a936b"),
];

struct TestCulture {
    organism: String,
    code: String,
    role: String,
    labs: StringSet,
    collections: StringSet,
}

struct TestCopy {
    culture: String,
    role: String,
    organism: String,
}

#[derive(Clone, Debug)]
pub struct LgdbReport {
    pub test_copies: usize,
    pub test_cultures: usize,
    pub b_rounds: usize,
    pub its_threshold_pct: Option<f64>,
    pub its_threshold_n: usize,
    pub tables: Vec<TableReport>,
}

impl LgdbReport {
    pub fn passed(&self) -> bool {
        self.tables.len() == TABLES.len() && self.tables.iter().all(|t| t.byte_identical && t.differing_cells == 0)
    }
    pub fn summary(&self) -> String {
        let mut s = format!("LongGloDB test: {} ASVs in {} cultures; condition B {} rounds; ITS limit {} from {} distances\n",
            self.test_copies, self.test_cultures, self.b_rounds,
            self.its_threshold_pct.map(py_float).unwrap_or_default(), self.its_threshold_n);
        for t in &self.tables {
            s.push_str(&format!("{}: {} rows; byte-identical={}; differing cells={}\n", t.table, t.rows, t.byte_identical, t.differing_cells));
            for d in &t.first_differences { s.push_str("  "); s.push_str(d); s.push('\n'); }
        }
        s.push_str(if self.passed() {
            "GATE 0b, LongGloDB conditions A and B, both VT analyses (species clades plus N6): PASS (exact archived output parity)\n"
        } else {
            "GATE 0b, LongGloDB: FAIL (see differences above)\n"
        });
        s
    }
}

fn pinned(path: &Path, expected: &str) -> Result<()> {
    let actual = sha256_file(path).with_context(|| format!("reading {}", path.display()))?;
    if actual != expected { bail!("{} hash mismatch: {actual}", path.display()); }
    Ok(())
}

/// lgdb_test.py `augmented_data`: the frozen data plus the test copies, which never act as comparators
/// (they stay out of `mains` and `main58`) and are held out of every call.
fn augment(base: &FrozenData, cultures: &BTreeMap<String, TestCulture>, copies: &[(String, TestCopy)],
           windows: &BTreeMap<String, SeqMap>, evidence: &Value, vt_key: &str) -> Result<FrozenData> {
    let mut d = base.clone();
    for (cid, c) in cultures {
        if d.cultures.contains_key(cid) { bail!("test culture ID collides with a reference culture: {cid}"); }
        let genus = normalise_genus(c.organism.split_whitespace().next().unwrap_or(""));
        let label = format!("LGDB:{}", c.code);
        if is_named(&c.organism) {
            d.species.insert(cid.clone(), c.organism.clone());
            d.by_species.entry(c.organism.clone()).or_default().insert(cid.clone());
        }
        d.genus.insert(cid.clone(), genus.clone());
        d.by_label.entry(label.clone()).or_default().insert(cid.clone());
        d.cultures.insert(cid.clone(), CultureMeta { label, organism: c.organism.clone(), genus });
    }
    for (x, t) in copies {
        if d.copies.contains_key(x) { bail!("test copy ID collides with a reference copy: {x}"); }
        let s58 = windows["s58"].get(x).ok_or_else(|| anyhow!("{x} has no s58 window"))?;
        d.copies.insert(x.clone(), CopyMeta { culture: t.culture.clone(), organism: t.organism.clone(), copy_class: "main".into(),
                                              div_kind: String::new(), s58_len: Some(s58.len()) });
        d.copy_order.push(x.clone());
        d.by_culture.entry(t.culture.clone()).or_default().push(x.clone());
        for (m, _) in WINDOWS {
            let seq = windows[m].get(x).ok_or_else(|| anyhow!("{x} has no {m} window"))?;
            d.fa.get_mut(m).ok_or_else(|| anyhow!("frozen data has no {m} sequences"))?.insert(x.clone(), seq.clone());
        }
        for (arm, key) in [("VT", vt_key), ("SH", "SH")] {
            let ev = evidence.get(key).and_then(|a| a.get(x)).ok_or_else(|| anyhow!("no {key} evidence for {x}"))?;
            let unit_pident: BTreeMap<String, f64> = ev.get("unit_pident").and_then(Value::as_object)
                .ok_or_else(|| anyhow!("{key} evidence of {x} has no unit_pident"))?
                .iter().map(|(u, p)| p.as_f64().map(|p| (u.clone(), p)).ok_or_else(|| anyhow!("{key}/{x}/{u}: not a number")))
                .collect::<Result<_>>()?;
            let anchor_pident = ev.get("anchor_pident").and_then(Value::as_f64);
            let best = ev.get("best").and_then(Value::as_array).filter(|b| b.len() == 2)
                .ok_or_else(|| anyhow!("{key} evidence of {x} has no best unit"))?;
            let unit = best[0].as_str().ok_or_else(|| anyhow!("{key}/{x}: best unit is not text"))?.to_owned();
            let pident = best[1].as_f64().ok_or_else(|| anyhow!("{key}/{x}: best identity is not a number"))?;
            d.evidence.get_mut(arm).ok_or_else(|| anyhow!("frozen data has no {arm} evidence"))?
                .queries.insert(x.clone(), EvidenceQuery { anchor_pident, unit_pident });
            d.best.get_mut(arm).ok_or_else(|| anyhow!("frozen data has no {arm} best hits"))?.insert(x.clone(), (unit, pident));
        }
    }
    Ok(d)
}

/// lgdb_test.py `placements_for`: the first table that places a copy wins; every copy must be placed.
fn placements_for(data: &FrozenData, tables: &[&Path]) -> Result<BTreeMap<String, Row>> {
    let mut out: BTreeMap<String, Row> = BTreeMap::new();
    for t in tables {
        for r in read_table(t)?.rows {
            let x = row_value(&r, "accession")?.to_owned();
            out.entry(x).or_insert(r);
        }
    }
    let missing = data.copies.keys().filter(|x| !out.contains_key(*x)).count();
    if missing > 0 { bail!("{missing} copies have no placement"); }
    out.retain(|x, _| data.copies.contains_key(x));
    Ok(out)
}

fn annotate(rows: &mut [Row], condition: &str, analysis: &str, fold: &str,
            cultures: &BTreeMap<String, TestCulture>, copies: &BTreeMap<&str, &TestCopy>) -> Result<()> {
    for r in rows.iter_mut() {
        let x = row_value(r, "accession")?.to_owned();
        let t = copies.get(x.as_str()).ok_or_else(|| anyhow!("{x} is not a test ASV"))?;
        let c = &cultures[&t.culture];
        r.insert("condition".into(), condition.into());
        r.insert("analysis".into(), analysis.into());
        r.insert("fold".into(), fold.into());
        r.insert("test_role".into(), t.role.clone());
        r.insert("test_code".into(), c.code.clone());
        r.insert("test_lab".into(), join_set(&c.labs, ";"));
        r.insert("test_collection".into(), join_set(&c.collections, ";"));
    }
    Ok(())
}

/// Gate 0b, LongGloDB. `benchmark_dir` is an amf-novelty-benchmark checkout (release 0.3.1 or later);
/// `lgdb_dir` holds the four Stage 1 windows of the test ASVs (`lgdb_test.{ssu_flank,full_its,s58,
/// LSU_LROR-FLR2}.fasta`) and `test_evidence.json.gz` from the run's `01_evidence` stage.
pub fn run_gate0b_lgdb(fixture_dir: impl AsRef<Path>, benchmark_dir: impl AsRef<Path>, lgdb_dir: impl AsRef<Path>) -> Result<LgdbReport> {
    let (fixture_dir, bench, lg) = (fixture_dir.as_ref(), benchmark_dir.as_ref(), lgdb_dir.as_ref());
    let runs = bench.join("benchmarks/m4_external_holdout/runs_lgdb");
    for (m, sha) in WINDOWS { pinned(&lg.join(format!("lgdb_test.{m}.fasta")), sha)?; }
    let evidence_path = lg.join("test_evidence.json.gz");
    pinned(&evidence_path, EVIDENCE_SHA)?;
    let overlap_path = bench.join("benchmarks/longglodb/overlap.asvs.tsv");
    pinned(&overlap_path, OVERLAP_SHA)?;
    let plan_path = runs.join("plan.json");
    pinned(&plan_path, PLAN_SHA)?;
    let a_path = runs.join("A/05_summary/lsu_placements.copies.tsv");
    pinned(&a_path, A_PLACEMENT_SHA)?;
    let frozen_path = fixture_dir.join("placement/lsu_placements.copies.tsv");
    pinned(&frozen_path, "1de103b6ae45f0372dda9ab0799230c4d1243cb98add943e5bab67785e5e910d")?;

    let plan: Value = serde_json::from_str(&fs::read_to_string(&plan_path)?)?;
    let mut b_rounds: Vec<(String, Vec<String>, PathBuf)> = Vec::new();
    for f in plan.get("b_folds").and_then(Value::as_array).ok_or_else(|| anyhow!("plan.json has no b_folds"))? {
        let key = get_string(f, "key")?.to_owned();
        let ids = f.get("query_ids").and_then(Value::as_array).ok_or_else(|| anyhow!("{key}: no query_ids"))?
            .iter().map(|v| v.as_str().map(str::to_owned).ok_or_else(|| anyhow!("{key}: query_ids"))).collect::<Result<Vec<_>>>()?;
        let path = runs.join("B").join(&key).join("05_summary/lsu_placements.copies.tsv");
        b_rounds.push((key, ids, path));
    }
    let mut lines = b_rounds.iter().map(|(k, _, p)| Ok(format!("{k}\t{}\n", sha256_file(p)?))).collect::<Result<Vec<_>>>()?;
    lines.sort();
    let digest = sha256_bytes(lines.concat().as_bytes());
    if digest != B_PLACEMENTS_DIGEST { bail!("condition B placement tables differ from the published run: {digest}"); }

    // test set (lgdb_test.py test_set): a culture is an organism and a culture code
    let mut cultures: BTreeMap<String, TestCulture> = BTreeMap::new();
    let mut copies: Vec<(String, TestCopy)> = Vec::new();
    for r in &read_table(&overlap_path)?.rows {
        let role = row_value(r, "role")?;
        if !TEST_ROLES.contains(&role) { continue; }
        let (organism, code) = (row_value(r, "organism")?, row_value(r, "code")?);
        let cid = format!("{organism} | LGDB {code}");
        let c = cultures.entry(cid.clone()).or_insert_with(|| TestCulture {
            organism: organism.into(), code: code.into(), role: role.into(), labs: StringSet::new(), collections: StringSet::new() });
        if c.role != role { bail!("culture {cid} has ASVs in two roles"); }
        c.labs.insert(row_value(r, "lab")?.into());
        c.collections.insert(row_value(r, "collection")?.into());
        copies.push((row_value(r, "asv")?.into(), TestCopy { culture: cid, role: role.into(), organism: organism.into() }));
    }
    let by_id: BTreeMap<&str, &TestCopy> = copies.iter().map(|(x, t)| (x.as_str(), t)).collect();
    if by_id.len() != copies.len() { bail!("duplicate test ASV IDs"); }
    let windows: BTreeMap<String, SeqMap> = WINDOWS.iter()
        .map(|(m, _)| Ok(((*m).to_owned(), read_fasta(&lg.join(format!("lgdb_test.{m}.fasta")))?))).collect::<Result<_>>()?;
    for (m, seqs) in &windows {
        if seqs.len() != by_id.len() || seqs.keys().any(|x| !by_id.contains_key(x.as_str())) {
            bail!("Stage 1 window {m} does not hold exactly the {} test ASVs", by_id.len());
        }
    }
    let evidence: Value = serde_json::from_str(&read_text_maybe_gzip(&evidence_path)?)?;

    let base = FrozenData::load(fixture_dir.join("inputs/benchmark"))?;
    let primary = augment(&base, &cultures, &copies, &windows, &evidence, "VT")?;
    let secondary = augment(&base, &cultures, &copies, &windows, &evidence, "VT_spanned")?;
    drop(base);
    // both analyses share the ITS sequences, so one distance cache serves both
    let mut rule = n6_its::ItsRule::new(&primary)?;
    let held: StringSet = cultures.keys().cloned().collect();
    let training: StringSet = primary.cultures.keys().filter(|c| !held.contains(*c)).cloned().collect();
    let ids_all: Vec<String> = by_id.keys().map(|x| (*x).to_owned()).collect();

    let mut tables: BTreeMap<String, Vec<Row>> = BTreeMap::new();
    for (analysis, d) in [("primary", &primary), ("secondary_vt", &secondary)] {
        let trainers: BTreeMap<String, Trainer> = ["VT", "SH"].iter()
            .map(|arm| Ok(((*arm).to_owned(), Trainer::new(d, arm)?))).collect::<Result<_>>()?;
        let mut run = |placements: &BTreeMap<String, Row>, ids: &[String], condition: &str, fold: &str| -> Result<()> {
            let mut v1 = name_round(d, &trainers, placements, &held, &training, ids)?;
            let before: BTreeMap<String, Row> = v1.iter().map(|r| Ok((row_value(r, "accession")?.to_owned(), r.clone()))).collect::<Result<_>>()?;
            let mut v2 = v1.clone();
            rule.apply(d, &mut v2, &before, &training, n6_its::Pool::All, 99)?;
            for r in v2.iter_mut() { genus_fields(d, r)?; }
            annotate(&mut v1, condition, analysis, fold, &cultures, &by_id)?;
            annotate(&mut v2, condition, analysis, fold, &cultures, &by_id)?;
            tables.entry(format!("{condition}_{analysis}_v1")).or_default().extend(v1);
            tables.entry(format!("{condition}_{analysis}_v2")).or_default().extend(v2);
            Ok(())
        };
        let pa = placements_for(d, &[&frozen_path, &a_path])?;
        run(&pa, &ids_all, "A", "A")?;
        for (key, ids, path) in &b_rounds {
            let pb = placements_for(d, &[path.as_path()])?;
            run(&pb, ids, "B", key)?;
        }
    }

    let mut reports = Vec::new();
    for (name, sha) in TABLES {
        let path = runs.join("06_scores").join(format!("{name}.copies.tsv"));
        pinned(&path, sha)?;
        let mut rows = tables.remove(name).ok_or_else(|| anyhow!("no {name} rows were made"))?;
        rows.sort_by(|a, b| (a.get("fold"), a.get("accession")).cmp(&(b.get("fold"), b.get("accession"))));
        reports.push(compare_table(&path, &rows, name)?);
    }
    let (t, n) = rule.threshold(&primary, &training, 99)?;
    Ok(LgdbReport { test_copies: copies.len(), test_cultures: cultures.len(), b_rounds: b_rounds.len(),
                    its_threshold_pct: t, its_threshold_n: n, tables: reports })
}
