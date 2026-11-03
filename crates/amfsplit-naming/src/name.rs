//! `db` and `name`: the paper's method for new sequences.
//!
//! `db` checks MaarjAM, UNITE and the frozen V18 reference folder against their pins, builds the VT
//! and SH BLAST databases once (UNITE filtered to Glomeromycota and renamed exactly as the frozen
//! loader does, so searches match the frozen ones), copies the V18 reference, and writes a manifest
//! with each tool's path, version and hash.
//!
//! `name` names query sequences the way the paper's external tests named theirs (lgdb_test.py,
//! condition A, primary analysis): every query joins the frozen culture data as a main-type copy of
//! its own held-out culture; VT and SH evidence come from BLAST+ with the culture copies' own
//! accessions dropped and each query's best unit by the external tests' rule; LSU placements come from
//! MAFFT and EPA-ng on the frozen V18 reference (the 903 culture copies keep their frozen placements);
//! the call is the species-clade arm plus the ITS distance rule (N6) against all 151 reference
//! cultures. Queries need all four windows from `amfsplit regions`: ssu_flank, full_its, s58 and
//! LSU_LROR-FLR2; any other query is listed in skipped.tsv with what it lacks.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};
use serde_json::Value;

use super::evidence::{best_unit, culture_accessions, load_reference, reduce_hits_with, run_blast_db, RefRecord, DATABASE_PINS, MAX_EVALUE};
use super::external::{augment, placements_for, TestCopy, TestCulture};
use super::placement::read_fasta_ordered;
use super::*;

const WINDOWS: [&str; 4] = ["ssu_flank", "full_its", "s58", "LSU_LROR-FLR2"];
/// VT and SH as the external tests searched them: arm, database kind, region, BLAST task, min identity, min coverage.
const NAME_ARMS: [(&str, &str, &str, &str, f64, f64); 2] = [
    ("VT", "maarjam", "ssu_flank", "megablast", 97.0, 0.95),
    ("SH", "unite", "full_its", "blastn", 98.5, 0.90),
];
/// The columns of a named row, in the frozen order, without those that only make sense when the
/// query's true culture is known (culture, organism, genus, group, correct, held_cultures and the
/// genus checks) or that are the same for every query (copy_class, div_kind).
pub const NAME_COLUMNS: [&str; 52] = [
    "accession", "VT_state", "VT_species", "VT_genera", "VT_units", "SH_state", "SH_species", "SH_genera", "SH_units",
    "LSU_state", "LSU_species", "LSU_genera", "LSU_units", "rank", "call_species", "call_genera", "unmatched_arms",
    "status", "novelty_claim", "artifact_claim", "placement", "reasons", "flags", "s58", "nearest_culture", "nearest_ssu_pct",
    "nearest_its_pct", "nearest_ssu_copy_count", "nearest_ssu_cultures", "artifact_comparator_its_pct", "best_vt", "best_vt_pident",
    "reconcile_rank", "sequence_type_signal", "d1_comparator_cultures", "d1_comparator_copies", "d1_short", "d2_paralog",
    "artifact_screen", "placement_rank_before_type_filter", "placement_species_before_type_filter", "lsu_best_edge_lwr",
    "lsu_credible_edge_count", "lsu_best_pendant_length", "lsu_retained_fraction", "lsu_external_species_clade",
    "lsu_external_species_clade_mass", "lsu_external_genus_clade", "lsu_external_genus_clade_mass", "its_nearest_training_pct",
    "its_threshold_pct", "its_rule_fired",
];
/// The short table: what was called, how sure, and why.
pub const CALL_COLUMNS: [&str; 16] = [
    "accession", "status", "placement", "novelty_claim", "artifact_claim", "reasons", "flags", "best_vt", "best_vt_pident",
    "VT_units", "SH_units", "LSU_units", "lsu_external_species_clade", "lsu_external_species_clade_mass",
    "its_nearest_training_pct", "its_rule_fired",
];

fn tool_record(tool: &str) -> Result<Value> {
    let path = which(tool)?;
    let version = std::process::Command::new(&path).arg(if tool.ends_with("blastn") || tool.ends_with("makeblastdb") { "-version" } else { "--version" })
        .output().ok().map(|o| format!("{}{}", String::from_utf8_lossy(&o.stdout), String::from_utf8_lossy(&o.stderr)))
        .and_then(|t| t.lines().map(str::trim).find(|l| !l.is_empty()).map(str::to_owned)).unwrap_or_default();
    Ok(serde_json::json!({"path": path.display().to_string(), "version": version, "sha256": sha256_file(&path)?}))
}

fn which(name: &str) -> Result<PathBuf> {
    let p = Path::new(name);
    if p.components().count() > 1 { return Ok(p.to_path_buf()); }
    std::env::var_os("PATH").and_then(|paths| std::env::split_paths(&paths).map(|d| d.join(name)).find(|c| c.is_file()))
        .ok_or_else(|| anyhow!("{name} is not on PATH"))
}

pub struct DbOptions {
    pub out: PathBuf,
    pub maarjam: PathBuf,
    pub unite: PathBuf,
    /// an existing frozen V18 reference folder, or
    pub reference_dir: Option<PathBuf>,
    /// the V18 full FASTA, from which the reference alignment is rebuilt (with the shipped tree,
    /// model and metadata in `shipped_v18`)
    pub v18_full: Option<PathBuf>,
    pub shipped_v18: PathBuf,
    pub mafft_threads: usize,
    pub blastn: String,
    pub makeblastdb: String,
    pub mafft: String,
    pub epa_ng: String,
    /// accept files that differ from the pinned ones (the run is then not the paper's method)
    pub allow_unpinned: bool,
}

/// Build the db folder; returns the manifest written to `out/manifest.json`.
pub fn build_db(o: &DbOptions) -> Result<Value> {
    if o.out.exists() && fs::read_dir(&o.out)?.next().is_some() { bail!("{} must be new or empty", o.out.display()); }
    let mut pinned = true;
    let mut check = |path: &Path, want: &str| -> Result<String> {
        let got = sha256_file(path)?;
        if got != want {
            if !o.allow_unpinned { bail!("{} is not the pinned file ({got}); pass --allow-unpinned to use it anyway", path.display()); }
            pinned = false;
        }
        Ok(got)
    };
    let pin = |arm: &str| DATABASE_PINS.iter().find(|(a, _, _)| *a == arm).copied().unwrap();
    let (_, vt_sha, vt_n) = pin("VT");
    let (_, sh_sha, sh_n) = pin("SH");
    let maarjam_sha = check(&o.maarjam, vt_sha)?;
    let unite_sha = check(&o.unite, sh_sha)?;
    let mut reference = serde_json::Map::new();
    let ref_out = o.out.join("v18_reference");
    fs::create_dir_all(&ref_out)?;
    match (&o.reference_dir, &o.v18_full) {
        (Some(dir), None) => {
            for (name, sha) in REFERENCE_PINS {
                reference.insert(name.into(), check(&dir.join(name), sha)?.into());
                fs::copy(dir.join(name), ref_out.join(name))?;
            }
        }
        (None, Some(full)) => {
            // lsu_placement.py build: R0001... in the release's order, uppercase, then MAFFT L-INS-i
            for (name, sha) in REFERENCE_PINS.iter().filter(|(n, _)| *n != "reference.aligned.fasta") {
                reference.insert((*name).into(), check(&o.shipped_v18.join(name), sha)?.into());
                fs::copy(o.shipped_v18.join(name), ref_out.join(name))?;
            }
            let meta = read_table(&o.shipped_v18.join("reference_metadata.tsv"))?;
            let seqs = read_fasta_ordered(full)?;
            let same = seqs.len() == meta.rows.len() && seqs.iter().zip(&meta.rows).all(|((x, _), r)| r.get("original_id") == Some(x));
            if !same { bail!("{} is not the V18 release the reference was built from (sequence IDs or order differ)", full.display()); }
            let mut unaligned = String::new();
            for ((_, seq), r) in seqs.iter().zip(&meta.rows) {
                unaligned.push('>'); unaligned.push_str(row_value(r, "ref_id")?); unaligned.push('\n'); unaligned.push_str(seq); unaligned.push('\n');
            }
            let un = ref_out.join("reference.unaligned.fasta");
            fs::write(&un, unaligned)?;
            let status = std::process::Command::new(which(&o.mafft)?)
                .args(["--localpair", "--maxiterate", "1000", "--thread", &o.mafft_threads.to_string(), "--threadit", "0"]).arg(&un)
                .stdout(fs::File::create(ref_out.join("reference.aligned.fasta"))?).stderr(fs::File::create(ref_out.join("reference_mafft.log"))?)
                .status()?;
            if !status.success() { bail!("MAFFT failed while rebuilding the reference alignment; see {}", ref_out.join("reference_mafft.log").display()); }
            let pin = REFERENCE_PINS.iter().find(|(n, _)| *n == "reference.aligned.fasta").unwrap().1;
            reference.insert("reference.aligned.fasta".into(), check(&ref_out.join("reference.aligned.fasta"), pin)?.into());
            reference.insert("built_from".into(), serde_json::json!({"v18_full": full.display().to_string(), "sha256": sha256_file(full)?}));
        }
        _ => bail!("give the V18 reference as --v18-full (the release FASTA) or --reference (a frozen reference folder)"),
    }
    let blast = o.out.join("blast");
    fs::create_dir_all(&blast)?;
    let mut arms = serde_json::Map::new();
    for ((arm, kind, region, task, min_id, min_cov), (path, records_pin)) in NAME_ARMS.iter().zip([(&o.maarjam, vt_n), (&o.unite, sh_n)]) {
        let (refs, fasta, read) = load_reference(path, kind, (*kind == "unite").then_some("Glomeromycota"))?;
        if pinned && read != records_pin { bail!("{arm}: {read} records read from {}, the pinned file has {records_pin}", path.display()); }
        let stem = blast.join(arm.to_lowercase());
        fs::write(stem.with_extension("fasta"), &fasta)?;
        let mut tsv = String::from("rid\tlabel\tacc\tn_N\n");
        for (rid, r) in &refs { tsv.push_str(&format!("{rid}\t{}\t{}\t{}\n", r.label, r.acc, r.n_n)); }
        fs::write(stem.with_extension("refs.tsv"), tsv)?;
        let status = std::process::Command::new(which(&o.makeblastdb)?).arg("-in").arg(stem.with_extension("fasta"))
            .args(["-dbtype", "nucl", "-out"]).arg(&stem).stdout(std::process::Stdio::null()).status()?;
        if !status.success() { bail!("makeblastdb failed for {arm} ({status})"); }
        arms.insert((*arm).into(), serde_json::json!({"database": path.display().to_string(), "records_read": read, "references": refs.len(),
            "region": region, "task": task, "min_id": min_id, "min_cov": min_cov, "max_evalue": MAX_EVALUE}));
    }
    let manifest = serde_json::json!({
        "schema": "amfsplit_db_v1", "pinned": pinned,
        "inputs": {"maarjam": {"path": o.maarjam.display().to_string(), "sha256": maarjam_sha},
                   "unite": {"path": o.unite.display().to_string(), "sha256": unite_sha}, "v18_reference": reference},
        "arms": arms,
        "tools": {"blastn": tool_record(&o.blastn)?, "makeblastdb": tool_record(&o.makeblastdb)?,
                  "mafft": tool_record(&o.mafft)?, "epa_ng": tool_record(&o.epa_ng)?},
    });
    fs::write(o.out.join("manifest.json"), serde_json::to_string_pretty(&manifest)? + "\n")?;
    Ok(manifest)
}

pub struct NameOptions {
    /// prefix of the `amfsplit regions` output: <prefix>.ssu_flank.fasta etc.
    pub windows: PathBuf,
    pub db: PathBuf,
    pub out: PathBuf,
    pub fixture_dir: PathBuf,
    pub blast_threads: usize,
    pub epa_threads: usize,
    /// tool overrides; default to the db manifest's
    pub blastn: Option<String>,
    pub mafft: Option<String>,
    pub epa_ng: Option<String>,
    /// a named-row table to compare with, on the shared columns
    pub compare_with: Option<PathBuf>,
    /// V (Protocol v3): VT coverage over the part of each MaarjAM reference the sequence spans
    pub vt_spanned: bool,
}

#[derive(Clone, Debug)]
pub struct NameReport {
    pub named: usize,
    pub skipped: usize,
    pub evidence_reused: bool,
    pub placement_reused: bool,
    pub status_counts: BTreeMap<String, usize>,
    /// (rows compared, differing cells by column)
    pub comparison: Option<(usize, BTreeMap<String, usize>)>,
}

fn windows_path(prefix: &Path, m: &str) -> PathBuf {
    let name = prefix.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
    prefix.with_file_name(format!("{name}.{m}.fasta"))
}

pub fn run_name(o: &NameOptions) -> Result<NameReport> {
    let bench = o.fixture_dir.join("inputs/benchmark");
    let work = o.out.join("work");
    fs::create_dir_all(&work)?;
    let manifest = || -> Result<Value> {
        let m: Value = serde_json::from_str(&fs::read_to_string(o.db.join("manifest.json")).with_context(|| format!("reading {}/manifest.json", o.db.display()))?)?;
        if m.get("schema").and_then(Value::as_str) != Some("amfsplit_db_v1") { bail!("{} is not an amfsplit db", o.db.display()); }
        Ok(m)
    };
    let tool = |m: &Value, key: &str, over: &Option<String>| -> Result<String> {
        over.clone().or_else(|| m.pointer(&format!("/tools/{key}/path")).and_then(Value::as_str).map(str::to_owned))
            .ok_or_else(|| anyhow!("no {key} in the db manifest"))
    };

    // windows: queries that have all four
    let windows: BTreeMap<String, SeqMap> = WINDOWS.iter().map(|m| Ok(((*m).to_owned(), read_fasta(&windows_path(&o.windows, m))?))).collect::<Result<_>>()?;
    let all: BTreeSet<&String> = windows.values().flat_map(|w| w.keys()).collect();
    let culture_copies = culture_accessions(&o.fixture_dir)?;
    let mut ids = Vec::new();
    let mut skipped = String::from("accession\tmissing_windows\n");
    let mut n_skipped = 0;
    for x in &all {
        let missing: Vec<&str> = WINDOWS.iter().copied().filter(|m| !windows[*m].contains_key(*x)).collect();
        if missing.is_empty() { ids.push((*x).clone()); } else { skipped.push_str(&format!("{x}\t{}\n", missing.join(","))); n_skipped += 1; }
    }
    if ids.is_empty() { bail!("no query has all four windows ({})", WINDOWS.join(", ")); }
    fs::write(o.out.join("skipped.tsv"), skipped)?;
    let write_queries = |m: &str, path: &Path| -> Result<()> {
        let mut s = String::new();
        for x in &ids { s.push('>'); s.push_str(x); s.push('\n'); s.push_str(&windows[m][x]); s.push('\n'); }
        fs::write(path, s).map_err(Into::into)
    };

    // 1. VT and SH evidence, cached per arm (work/evidence.json; work/evidence_vt_spanned.json under V)
    let cache = work.join(if o.vt_spanned { "evidence_vt_spanned.json" } else { "evidence.json" });
    let mut arms: serde_json::Map<String, Value> = if cache.is_file() {
        match serde_json::from_str::<Value>(&fs::read_to_string(&cache)?)? { Value::Object(m) => m, _ => bail!("{} is not an evidence cache", cache.display()) }
    } else { serde_json::Map::new() };
    let evidence_reused = NAME_ARMS.iter().all(|(arm, ..)| arms.contains_key(*arm));
    if !evidence_reused {
        let m = manifest()?;
        let blastn = tool(&m, "blastn", &o.blastn)?;
        let drop = culture_copies.clone();
        for (arm, _, region, task, min_id, min_cov) in NAME_ARMS {
            if arms.contains_key(arm) { continue; }
            let spanned = o.vt_spanned && arm == "VT";
            let stem = o.db.join("blast").join(arm.to_lowercase());
            let refs: BTreeMap<String, RefRecord> = read_table(&stem.with_extension("refs.tsv"))?.rows.iter().map(|r| Ok((row_value(r, "rid")?.to_owned(),
                RefRecord { label: row_value(r, "label")?.to_owned(), acc: row_value(r, "acc")?.to_owned(), n_n: row_value(r, "n_N")?.parse()? }))).collect::<Result<_>>()?;
            let qpath = work.join(format!("{}.queries.fasta", arm.to_lowercase()));
            write_queries(region, &qpath)?;
            let hits = run_blast_db(&qpath, &stem, task, o.blast_threads, &blastn, spanned)?;
            let mut qs = serde_json::Map::new();
            for x in &ids {
                let h = hits.get(x).map(Vec::as_slice).unwrap_or(&[]);
                let (units, anchor) = reduce_hits_with(h, &refs, min_cov, MAX_EVALUE, min_id, &drop, spanned)?;
                let (bu, bp) = best_unit(&units, h, &refs);
                let unit_pident: BTreeMap<String, f64> = units.into_iter().collect();
                qs.insert(x.clone(), serde_json::json!({"unit_pident": unit_pident, "anchor_pident": anchor, "best": [bu, bp]}));
            }
            arms.insert(arm.into(), Value::Object(qs));
        }
        fs::write(&cache, serde_json::to_string(&Value::Object(arms.clone()))?)?;
    }
    let evidence: Value = serde_json::from_str(&fs::read_to_string(&cache)?)?;

    // 2. LSU placements of the queries (cached in work/place)
    let place_dir = work.join("place");
    let table = place_dir.join("lsu_placements.copies.tsv");
    let placement_reused = table.is_file();
    if !placement_reused {
        let m = manifest()?;
        let q = work.join("lsu.queries.fasta");
        write_queries("LSU_LROR-FLR2", &q)?;
        run_place(&PlaceOptions {
            reference_dir: o.db.join("v18_reference"), queries: q, out: place_dir.clone(),
            mafft: tool(&m, "mafft", &o.mafft)?, epa_ng: tool(&m, "epa_ng", &o.epa_ng)?, mafft_threads: 4, epa_threads: o.epa_threads,
            mass: 0.95, allow_other_reference: !m.get("pinned").and_then(Value::as_bool).unwrap_or(false), compare_with: None,
        })?;
    }

    // 3. the call: each query its own held-out culture, named against all reference cultures
    let base = FrozenData::load(&bench)?;
    let mut cultures = BTreeMap::new();
    let mut copies = Vec::new();
    for x in &ids {
        let cid = format!("{x} | query");
        cultures.insert(cid.clone(), TestCulture { organism: String::new(), code: x.clone(), role: String::new(), labs: StringSet::new(), collections: StringSet::new() });
        copies.push((x.clone(), TestCopy { culture: cid, role: String::new(), organism: String::new() }));
    }
    let data = augment(&base, &cultures, &copies, &windows, &evidence, "VT")?;
    drop(base);
    let held: StringSet = cultures.keys().cloned().collect();
    let training: StringSet = data.cultures.keys().filter(|c| !held.contains(*c)).cloned().collect();
    let frozen_placements = o.fixture_dir.join("placement/lsu_placements.copies.tsv");
    let placements = placements_for(&data, &[frozen_placements.as_path(), table.as_path()])?;
    let trainers: BTreeMap<String, Trainer> = ["VT", "SH"].iter().map(|arm| Ok(((*arm).to_owned(), Trainer::new(&data, arm)?))).collect::<Result<_>>()?;
    let v1 = name_round(&data, &trainers, &placements, &held, &training, &ids)?;
    let before: BTreeMap<String, Row> = v1.iter().map(|r| Ok((row_value(r, "accession")?.to_owned(), r.clone()))).collect::<Result<_>>()?;
    let mut rows = v1;
    let mut rule = n6_its::ItsRule::new_skipping_held_pairs(&data, &held)?;
    rule.apply(&data, &mut rows, &before, &training, n6_its::Pool::All, 99)?;

    let headers: Vec<String> = NAME_COLUMNS.iter().map(|h| (*h).to_owned()).collect();
    fs::write(o.out.join("rows.tsv"), serialize_table(&headers, &rows))?;
    let call_headers: Vec<String> = CALL_COLUMNS.iter().map(|h| (*h).to_owned()).collect();
    fs::write(o.out.join("calls.tsv"), serialize_table(&call_headers, &rows))?;
    let mut status_counts = BTreeMap::new();
    for r in &rows { *status_counts.entry(r.get("status").cloned().unwrap_or_default()).or_insert(0) += 1; }

    let comparison = match &o.compare_with {
        None => None,
        Some(path) => {
            let other = read_row_map(&read_table_any(path)?);
            let (mut n, mut diff) = (0usize, BTreeMap::new());
            for r in &rows {
                let x = row_value(r, "accession")?;
                let Some(e) = other.get(x) else { continue };
                n += 1;
                for h in NAME_COLUMNS.iter().filter(|h| e.contains_key(**h)) {
                    if r.get(*h) != e.get(*h) { *diff.entry((*h).to_owned()).or_insert(0) += 1; }
                }
            }
            Some((n, diff))
        }
    };
    let run = serde_json::json!({
        "windows": WINDOWS.iter().map(|m| (m.to_string(), sha256_file(&windows_path(&o.windows, m)).unwrap_or_default())).collect::<BTreeMap<_, _>>(),
        "db": o.db.display().to_string(), "db_manifest_sha256": sha256_file(&o.db.join("manifest.json")).ok(),
        "named": rows.len(), "skipped": n_skipped, "status_counts": status_counts,
        "evidence_reused": evidence_reused, "placement_reused": placement_reused,
        "method": "species clades plus the ITS distance rule (N6), each query its own held-out culture, as the external tests (condition A, primary VT)",
    });
    fs::write(o.out.join("run.json"), serde_json::to_string_pretty(&run)? + "\n")?;
    Ok(NameReport { named: rows.len(), skipped: n_skipped, evidence_reused, placement_reused, status_counts, comparison })
}
