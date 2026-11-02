//! Gate 0b, the two external tests of amf-novelty-benchmark reproduced in Rust: LongGloDB
//! (`lgdb_test.py`, PROTOCOL_LGDB.md with amendments 1 and 2) and EUKARYOME (`euk_test.py`, the same
//! runner on another test set, PROTOCOL_EUK.md). What follows describes LongGloDB; EUKARYOME differs
//! only in its test set (5,511 reads, one culture per labelled taxon), its 51 condition B rounds
//! (30 species, 21 genus), and in publishing its evidence, placements and tables gzip-compressed.
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
const WINDOW_NAMES: [&str; 4] = ["ssu_flank", "full_its", "s58", "LSU_LROR-FLR2"];

/// What names one external test and pins its inputs to the published run.
struct Spec {
    name: &'static str,
    /// test set, relative to the benchmark checkout
    overlap: &'static str,
    overlap_sha: &'static str,
    /// run folder, relative to benchmarks/m4_external_holdout
    runs: &'static str,
    plan_sha: &'static str,
    /// culture IDs: "organism | LGDB code" (LongGloDB) or "organism | EUKARYOME"
    one_culture_per_taxon: bool,
    /// window file name in the local folder, with {m} for the region
    window_file: &'static str,
    window_shas: [&'static str; 4],
    /// test evidence: in the local folder (LongGloDB) or in the run's 01_evidence (EUKARYOME)
    evidence_in_runs: bool,
    evidence_sha: &'static str,
    /// placements and tables published as .tsv.gz
    gz: bool,
    a_placement_sha: &'static str,
    /// SHA-256 of the sorted lines "<key>\t<SHA-256 of the round's placement file>\n"
    b_placements_digest: &'static str,
    table_shas: [&'static str; 8],
}

const TABLE_NAMES: [&str; 8] = ["A_primary_v1", "A_primary_v2", "A_secondary_vt_v1", "A_secondary_vt_v2",
                                "B_primary_v1", "B_primary_v2", "B_secondary_vt_v1", "B_secondary_vt_v2"];

const LGDB: Spec = Spec {
    name: "LongGloDB",
    overlap: "benchmarks/longglodb/overlap.asvs.tsv",
    overlap_sha: "e39919cc5118abd0a24ed738028b8032fdc4aba497be78bdc6deb2ebe8096501",
    runs: "runs_lgdb",
    plan_sha: "821b38074cf045bb9eb099b33a48215c16cc98c482d396febc1d372bf42e44cb",
    one_culture_per_taxon: false,
    window_file: "lgdb_test.{m}.fasta",
    window_shas: ["d76e02fd1ee62d649de794e0519c8fae15c5d39229c441b4b7fc3d0538be4d9e",
                  "bc1ae16bf4d7491e6c0be7f0c138f1200ca38282dab726a2bcb9f7c8be6d6465",
                  "46da4e66960256f86066a7d7670c4316a3d78b380395dabc25bb05dfb20886f7",
                  "82b147d60e796c4d172702c9f3f7979f0b660a5015bca8d827e7bd151f9ca190"],
    evidence_in_runs: false,
    evidence_sha: "543f61883c4116ea1107495da67588bc47725900c0d5c6889e183762a02cf5fe",
    gz: false,
    a_placement_sha: "6eae1912f400af0954c09beb421ebdc561157fef657d02511c61e90d7e889805",
    b_placements_digest: "1b6bddb263ddca0e9f7b5a5ca17c693219eb80f31bc2698609b0e41db4fe5b5e",
    table_shas: ["32e8c0dc5bd88fe89f31e4fb39bcaae4f71c5cc55a6585a29b219953a231c047",
                 "f9fab9e8ec2a64042b96db3e6a387946a0b6d6b04abdde6892d3b3676149b4a0",
                 "83732ef1f5126714e446a041a7487a9b8c8c1ed2d95773be91eab344df76c4c0",
                 "2010f464b93cd8a67ef783e6ea3d136f40288a3cf43e67a6036d598512cac26e",
                 "c44d6bebcb19b4a3942f424d960a856af433905f9e3373328dfa2b9d88b8d3fa",
                 "b6b1dec756b3d6c99b89c77bab8ce60c8ffe1f728bddc02f9269d1408ae3225c",
                 "b84248ca32c8564248597b5e57311a4cc1dfeec0237058e1883792de491297cb",
                 "ced1c3cf61ace3c3fb3f685a86b483d848d400e427e318b9614452f5d08a936b"],
};

const EUK: Spec = Spec {
    name: "EUKARYOME",
    overlap: "benchmarks/eukaryome/overlap.asvs.tsv",
    overlap_sha: "ada9f7fa212c30c2d975eca3bf483413d9d1bb7c290b99038a89e7dc28ec479f",
    runs: "runs_euk",
    plan_sha: "373d6929ab794dd32881b714e28226c810bc2aa449aac45fc7c9de11a39c555e",
    one_culture_per_taxon: true,
    window_file: "overlap.win.{m}.fasta",
    window_shas: ["e67f99c47974feec819975d6e01589ab2b63387022682c573c9954d8bbfad58c",
                  "b52392c1aa0d6d4160ef8c5ebf428c3a4c8c0ecb10945dd3497b5df045bdc3e9",
                  "134a12265dfb1ef85631ba82203670aacac11839e9b7669ec6441f80e5742e87",
                  "b6ef33a9890527dd25f420b9209375b258055181c6ff6c9fd8656f0abcda277a"],
    evidence_in_runs: true,
    evidence_sha: "ae1d66e284975cdb5dd56cb80eb393bc1a4f8905aedc39d6a95adc838b41cbd0",
    gz: true,
    a_placement_sha: "04f7066d193fa706fa7868d055fd1dae27d7accf721a0429ca639e1e5fbad40f",
    b_placements_digest: "4c138cc4d341f3afc6175c7e80c4f6f4cfefa1d397cd36d612767fecfcbb1101",
    table_shas: ["56ee075859c007774e0cc8685dded31bd6330417ac648764179346ec0f5e43c6",
                 "15a76e4d9b97d016ff4cb6040c4910b64fc572c3062f13c1f22a2d65f5d756a6",
                 "d44a5ec66c5782560e493f11148d95ad819e5a4a04a85f956e5623f5eb547d11",
                 "082295d29a83dc40efa6dcf051877db13a9755346bb42490a411d1d4fc5e7b60",
                 "1ab427cf027dd3771f9f9abd405620e08b47b5ea5be385e9769c3744fd86bdd3",
                 "eec47fd28db6b35ec32adb1b858c9655db309b11cb69b3ee98c8833cdfaae93c",
                 "b0da94c0c2da2376babaed9c8a69f9a481c3518be81e077f71c3449c6227715f",
                 "9728271f313522aedf70ec40fcccf6c7ebff5b13daef7173dfbd5a566cd46266"],
};

pub(crate) struct TestCulture {
    pub(crate) organism: String,
    pub(crate) code: String,
    pub(crate) role: String,
    pub(crate) labs: StringSet,
    pub(crate) collections: StringSet,
}

pub(crate) struct TestCopy {
    pub(crate) culture: String,
    pub(crate) role: String,
    pub(crate) organism: String,
}

#[derive(Clone, Debug)]
pub struct ExternalReport {
    pub test: String,
    pub test_copies: usize,
    pub test_cultures: usize,
    pub b_rounds: usize,
    pub its_threshold_pct: Option<f64>,
    pub its_threshold_n: usize,
    pub tables: Vec<TableReport>,
}

impl ExternalReport {
    pub fn passed(&self) -> bool {
        self.tables.len() == TABLE_NAMES.len() && self.tables.iter().all(|t| t.byte_identical && t.differing_cells == 0)
    }
    pub fn summary(&self) -> String {
        let mut s = format!("{} test: {} sequences in {} cultures; condition B {} rounds; ITS limit {} from {} distances\n",
            self.test, self.test_copies, self.test_cultures, self.b_rounds,
            self.its_threshold_pct.map(py_float).unwrap_or_default(), self.its_threshold_n);
        for t in &self.tables {
            s.push_str(&format!("{}: {} rows; byte-identical={}; differing cells={}\n", t.table, t.rows, t.byte_identical, t.differing_cells));
            for d in &t.first_differences { s.push_str("  "); s.push_str(d); s.push('\n'); }
        }
        s.push_str(if self.passed() {
            format!("GATE 0b, {} conditions A and B, both VT analyses (species clades plus N6): PASS (exact archived output parity)\n", self.test)
        } else {
            format!("GATE 0b, {}: FAIL (see differences above)\n", self.test)
        }.as_str());
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
pub(crate) fn augment(base: &FrozenData, cultures: &BTreeMap<String, TestCulture>, copies: &[(String, TestCopy)],
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
        for m in WINDOW_NAMES {
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
pub(crate) fn placements_for(data: &FrozenData, tables: &[&Path]) -> Result<BTreeMap<String, Row>> {
    let mut out: BTreeMap<String, Row> = BTreeMap::new();
    for t in tables {
        for r in read_table_any(t)?.rows {
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
pub fn run_gate0b_lgdb(fixture_dir: impl AsRef<Path>, benchmark_dir: impl AsRef<Path>, lgdb_dir: impl AsRef<Path>) -> Result<ExternalReport> {
    run_external(&LGDB, fixture_dir.as_ref(), benchmark_dir.as_ref(), lgdb_dir.as_ref())
}

/// Gate 0b, EUKARYOME. `euk_dir` holds the four test windows of the 5,511 reads
/// (`overlap.win.{ssu_flank,full_its,s58,LSU_LROR-FLR2}.fasta`); everything else, the test evidence
/// included, comes from the benchmark checkout.
pub fn run_gate0b_euk(fixture_dir: impl AsRef<Path>, benchmark_dir: impl AsRef<Path>, euk_dir: impl AsRef<Path>) -> Result<ExternalReport> {
    run_external(&EUK, fixture_dir.as_ref(), benchmark_dir.as_ref(), euk_dir.as_ref())
}

fn run_external(spec: &Spec, fixture_dir: &Path, bench: &Path, lg: &Path) -> Result<ExternalReport> {
    let gz = if spec.gz { ".gz" } else { "" };
    let runs = bench.join("benchmarks/m4_external_holdout").join(spec.runs);
    for (m, sha) in WINDOW_NAMES.iter().zip(spec.window_shas) { pinned(&lg.join(spec.window_file.replace("{m}", m)), sha)?; }
    let evidence_path = if spec.evidence_in_runs { runs.join("01_evidence/test_evidence.json.gz") } else { lg.join("test_evidence.json.gz") };
    pinned(&evidence_path, spec.evidence_sha)?;
    let overlap_path = bench.join(spec.overlap);
    pinned(&overlap_path, spec.overlap_sha)?;
    let plan_path = runs.join("plan.json");
    pinned(&plan_path, spec.plan_sha)?;
    let a_path = runs.join(format!("A/05_summary/lsu_placements.copies.tsv{gz}"));
    pinned(&a_path, spec.a_placement_sha)?;
    let frozen_path = fixture_dir.join("placement/lsu_placements.copies.tsv");
    pinned(&frozen_path, "1de103b6ae45f0372dda9ab0799230c4d1243cb98add943e5bab67785e5e910d")?;

    let plan: Value = serde_json::from_str(&fs::read_to_string(&plan_path)?)?;
    let mut b_rounds: Vec<(String, Vec<String>, PathBuf)> = Vec::new();
    for f in plan.get("b_folds").and_then(Value::as_array).ok_or_else(|| anyhow!("plan.json has no b_folds"))? {
        let key = get_string(f, "key")?.to_owned();
        let ids = f.get("query_ids").and_then(Value::as_array).ok_or_else(|| anyhow!("{key}: no query_ids"))?
            .iter().map(|v| v.as_str().map(str::to_owned).ok_or_else(|| anyhow!("{key}: query_ids"))).collect::<Result<Vec<_>>>()?;
        let path = runs.join("B").join(&key).join(format!("05_summary/lsu_placements.copies.tsv{gz}"));
        b_rounds.push((key, ids, path));
    }
    let mut lines = b_rounds.iter().map(|(k, _, p)| Ok(format!("{k}\t{}\n", sha256_file(p)?))).collect::<Result<Vec<_>>>()?;
    lines.sort();
    let digest = sha256_bytes(lines.concat().as_bytes());
    if digest != spec.b_placements_digest { bail!("condition B placement tables differ from the published run: {digest}"); }

    // test set (lgdb_test.py test_set): a culture is an organism and a culture code
    let mut cultures: BTreeMap<String, TestCulture> = BTreeMap::new();
    let mut copies: Vec<(String, TestCopy)> = Vec::new();
    for r in &read_table(&overlap_path)?.rows {
        let role = row_value(r, "role")?;
        if !TEST_ROLES.contains(&role) { continue; }
        let (organism, code) = (row_value(r, "organism")?, row_value(r, "code")?);
        let cid = if spec.one_culture_per_taxon { format!("{organism} | EUKARYOME") } else { format!("{organism} | LGDB {code}") };
        let c = cultures.entry(cid.clone()).or_insert_with(|| TestCulture {
            organism: organism.into(), code: code.into(), role: role.into(), labs: StringSet::new(), collections: StringSet::new() });
        if c.role != role { bail!("culture {cid} has ASVs in two roles"); }
        c.labs.insert(row_value(r, "lab")?.into());
        c.collections.insert(row_value(r, "collection")?.into());
        copies.push((row_value(r, "asv")?.into(), TestCopy { culture: cid, role: role.into(), organism: organism.into() }));
    }
    let by_id: BTreeMap<&str, &TestCopy> = copies.iter().map(|(x, t)| (x.as_str(), t)).collect();
    if by_id.len() != copies.len() { bail!("duplicate test ASV IDs"); }
    let windows: BTreeMap<String, SeqMap> = WINDOW_NAMES.iter()
        .map(|m| Ok(((*m).to_owned(), read_fasta(&lg.join(spec.window_file.replace("{m}", m)))?))).collect::<Result<_>>()?;
    for (m, seqs) in &windows {
        if seqs.len() != by_id.len() || seqs.keys().any(|x| !by_id.contains_key(x.as_str())) {
            bail!("test window {m} does not hold exactly the {} test sequences", by_id.len());
        }
    }
    let evidence: Value = serde_json::from_str(&read_text_maybe_gzip(&evidence_path)?)?;

    let base = FrozenData::load(fixture_dir.join("inputs/benchmark"))?;
    let primary = augment(&base, &cultures, &copies, &windows, &evidence, "VT")?;
    let secondary = augment(&base, &cultures, &copies, &windows, &evidence, "VT_spanned")?;
    drop(base);
    // both analyses share the ITS sequences, so one distance cache serves both; pairs of two test
    // sequences are skipped, since every call holds all test cultures out
    let held: StringSet = cultures.keys().cloned().collect();
    let mut rule = n6_its::ItsRule::new_skipping_held_pairs(&primary, &held)?;
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
    for (name, sha) in TABLE_NAMES.iter().zip(spec.table_shas) {
        let path = runs.join("06_scores").join(format!("{name}.copies.tsv{gz}"));
        pinned(&path, sha)?;
        let mut rows = tables.remove(*name).ok_or_else(|| anyhow!("no {name} rows were made"))?;
        rows.sort_by(|a, b| (a.get("fold"), a.get("accession")).cmp(&(b.get("fold"), b.get("accession"))));
        reports.push(compare_table(&path, &rows, name)?);
    }
    let (t, n) = rule.threshold(&primary, &training, 99)?;
    Ok(ExternalReport { test: spec.name.into(), test_copies: copies.len(), test_cultures: cultures.len(), b_rounds: b_rounds.len(),
                    its_threshold_pct: t, its_threshold_n: n, tables: reports })
}
