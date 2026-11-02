//! VT, SH and LSU evidence: the frozen `arms.py` (`load_reference`, `run_blast`) and `arm_sets.py`
//! (`vt_hits`) in Rust. Each query window is searched with BLAST+ against MaarjAM (VT, SSU),
//! UNITE (SH, full ITS, Glomeromycota records only) or the V18 AMF-only LSU database, and its hits
//! become, per unit, the best identity among hits that pass coverage and e-value (`unit_pident`),
//! plus the identity of the highest-bit-score hit that also passes the arm's minimum identity
//! (`anchor_pident`). Hits on records whose accession is one of the culture reference's own are
//! dropped, so a culture never supports itself.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};
use flate2::read::GzDecoder;
use serde_json::Value;

use super::*;

/// The three arms as the frozen ARMS table defines them: database kind, query region, BLAST task,
/// minimum identity, minimum coverage.
pub const EVIDENCE_ARMS: [(&str, &str, &str, &str, f64, f64); 3] = [
    ("VT", "maarjam", "ssu_flank", "megablast", 97.0, 0.95),
    ("SH", "unite", "full_its", "blastn", 98.5, 0.90),
    ("LSU", "delavaux", "LSU_LROR-FLR2", "blastn", 97.0, 0.90),
];
pub const MAX_EVALUE: f64 = 1e-50;
/// The databases behind the frozen evidence (crosswalk_v4.run.json): SHA-256 and records read.
pub const DATABASE_PINS: [(&str, &str, usize); 3] = [
    ("VT", "2b2c57e664e763473a23c8eba12de07d448cff5b45b1fb75794ca5701930969a", 3081),
    ("SH", "abe98f43e16c179f2c129f137989cad48e6855ce4be9ae0f365defd1f856369b", 2069189),
    ("LSU", "e6794cbbe166735a94e44b146aae2ecffb35ebfee61cf6ae8ed2e48ba9227198", 460),
];

#[derive(Clone, Debug)]
pub struct RefRecord {
    pub label: String,
    pub acc: String,
    pub n_n: i64,
}

#[derive(Clone, Debug)]
pub struct Hit {
    pub reference: String,
    pub pident: f64,
    pub length: i64,
    pub qlen: i64,
    pub slen: i64,
    pub evalue: f64,
    pub bitscore: f64,
}

/// `ACC.match`: 1 or 2 capitals, an optional underscore and at least five digits at the start.
fn acc_prefix(f: &str) -> bool {
    let b = f.as_bytes();
    for k in [2usize, 1] {
        if b.len() < k || !b[..k].iter().all(u8::is_ascii_uppercase) { continue; }
        let mut j = k;
        if b.get(j) == Some(&b'_') { j += 1; }
        if b[j..].iter().take_while(|c| c.is_ascii_digit()).count() >= 5 { return true; }
    }
    false
}

/// First `SH\d+\.\d+FU` in the text.
fn sh_code(s: &str) -> Option<&str> {
    let b = s.as_bytes();
    for i in 0..b.len().saturating_sub(1) {
        if &b[i..i + 2] != b"SH" { continue; }
        let mut j = i + 2;
        let d1 = b[j..].iter().take_while(|c| c.is_ascii_digit()).count();
        if d1 == 0 { continue; }
        j += d1;
        if b.get(j) != Some(&b'.') { continue; }
        j += 1;
        let d2 = b[j..].iter().take_while(|c| c.is_ascii_digit()).count();
        if d2 == 0 { continue; }
        j += d2;
        if b.get(j..j + 2) == Some(b"FU") { return Some(&s[i..j + 2]); }
    }
    None
}

fn is_vtx(s: &str) -> bool { s.strip_prefix("VTX").is_some_and(|d| !d.is_empty() && d.bytes().all(|c| c.is_ascii_digit())) }

fn parse_maarjam(h: &str) -> Result<(String, String)> {
    // 'VTX00223__AB555676__Glomeraceae_Glomus_sp.'
    if let Some(rest) = h.strip_prefix("VTX") {
        let digits = rest.bytes().take_while(u8::is_ascii_digit).count();
        if digits > 0 {
            if let Some(after) = rest[digits..].strip_prefix("__") {
                let acc_len = after.find('_').unwrap_or(after.len());
                if acc_len > 0 {
                    if let Some(tail) = after[acc_len..].strip_prefix("__") {
                        if tail.chars().next().is_some_and(|c| !c.is_whitespace()) {
                            return Ok((h[..3 + digits].to_owned(), after[..acc_len].to_owned()));
                        }
                    }
                }
            }
        }
    }
    // 'gi|MOA01016|gb|AB555676| Glomeraceae Glomus sp. VTX00223'
    let (idpart, rest) = h.split_once(' ').unwrap_or((h, ""));
    let p: Vec<&str> = idpart.split('|').collect();
    let toks: Vec<&str> = rest.split_whitespace().collect();
    match toks.last() {
        Some(last) if p.len() >= 4 && is_vtx(last) => Ok(((*last).to_owned(), p[3].to_owned())),
        _ => bail!("unexpected MaarjAM header: {h}"),
    }
}

fn parse_unite(h: &str) -> Result<(String, String)> {
    let sh = sh_code(h).ok_or_else(|| anyhow!("no SH code in UNITE header: {h}"))?;
    let acc = h.split(|c: char| c == '|' || c.is_whitespace())
        .find(|f| (acc_prefix(f) || f.starts_with("UDB")) && sh_code(f) != Some(*f)).unwrap_or("");
    Ok((sh.to_owned(), acc.to_owned()))
}

fn parse_delavaux(h: &str) -> Result<(String, String)> {
    let first = h.split_whitespace().next().unwrap_or("");
    let b = first.as_bytes();
    let mut found = None;
    'outer: for k in [2usize, 1] {
        if b.len() < k || !b[..k].iter().all(u8::is_ascii_uppercase) { continue; }
        let mut j = k;
        if b.get(j) == Some(&b'_') { j += 1; }
        let d = b[j..].iter().take_while(|c| c.is_ascii_digit()).count();
        if d == 0 { continue; }
        j += d;
        if b.get(j) == Some(&b'.') {
            let d2 = b[j + 1..].iter().take_while(|c| c.is_ascii_digit()).count();
            if d2 > 0 { j += 1 + d2; }
        }
        if b.get(j) == Some(&b'_') && j + 1 < b.len() { found = Some(j); break 'outer; }
    }
    let j = found.ok_or_else(|| anyhow!("unexpected Delavaux header: {h}"))?;
    let toks: Vec<&str> = first[j + 1..].split('_').collect();
    let rest = if toks[0].ends_with("aceae") { &toks[1..] } else { &toks[..] };
    if rest.is_empty() { bail!("no genus in Delavaux header: {h}"); }
    Ok((rest.join(" "), first[..j].to_owned()))
}

/// arms.py `load_reference`: records renamed r0000001, ... in file order (only those whose header
/// contains `keep`, when given); returns the records, the renamed FASTA and the number read.
pub fn load_reference(path: &Path, kind: &str, keep: Option<&str>) -> Result<(BTreeMap<String, RefRecord>, String, usize)> {
    let file = fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let reader: Box<dyn BufRead> = if path.extension().and_then(|x| x.to_str()) == Some("gz") {
        Box::new(std::io::BufReader::new(GzDecoder::new(file)))
    } else { Box::new(std::io::BufReader::new(file)) };
    let (mut refs, mut out, mut n) = (BTreeMap::new(), String::new(), 0usize);
    let mut push = |header: &str, seq: &str, refs: &mut BTreeMap<String, RefRecord>| -> Result<()> {
        n += 1;
        if keep.is_some_and(|k| !header.contains(k)) { return Ok(()); }
        let rid = format!("r{:07}", refs.len() + 1);
        let (label, acc) = match kind { "maarjam" => parse_maarjam(header)?, "unite" => parse_unite(header)?,
                                        "delavaux" => parse_delavaux(header)?, _ => bail!("unknown reference kind {kind}") };
        let n_n = seq.bytes().filter(|c| c.eq_ignore_ascii_case(&b'N')).count() as i64;
        out.push('>'); out.push_str(&rid); out.push('\n'); out.push_str(seq); out.push('\n');
        refs.insert(rid, RefRecord { label, acc, n_n });
        Ok(())
    };
    let (mut header, mut seq): (Option<String>, String) = (None, String::new());
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if let Some(h) = line.strip_prefix('>') {
            if let Some(prev) = header.take() { push(&prev, &seq, &mut refs)?; }
            header = Some(h.to_owned());
            seq.clear();
        } else if !line.is_empty() {
            seq.push_str(line);
        }
    }
    if let Some(prev) = header { push(&prev, &seq, &mut refs)?; }
    if refs.is_empty() { bail!("no sequences in {}{}", path.display(), keep.map(|k| format!(" match {k:?}")).unwrap_or_default()); }
    Ok((refs, out, n))
}

/// arms.py `run_blast`: makeblastdb, then blastn with one HSP per subject.
pub fn run_blast(query: &Path, ref_fasta: &str, workdir: &Path, name: &str, task: &str, threads: usize,
                 makeblastdb: &str, blastn: &str) -> Result<BTreeMap<String, Vec<Hit>>> {
    fs::create_dir_all(workdir)?;
    let db = workdir.join(name);
    let fasta = workdir.join(format!("{name}.fasta"));
    fs::write(&fasta, ref_fasta)?;
    let status = std::process::Command::new(makeblastdb).args(["-in"]).arg(&fasta).args(["-dbtype", "nucl", "-out"]).arg(&db)
        .stdout(std::process::Stdio::null()).status().with_context(|| format!("starting {makeblastdb}"))?;
    if !status.success() { bail!("makeblastdb failed ({status})"); }
    let out = std::process::Command::new(blastn)
        .args(["-task", task, "-query"]).arg(query).arg("-db").arg(&db)
        .args(["-dust", "no", "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore", "-evalue", "1e-10",
               "-max_target_seqs", "5000", "-max_hsps", "1", "-num_threads", &threads.to_string()])
        .output().with_context(|| format!("starting {blastn}"))?;
    if !out.status.success() { bail!("blastn failed: {}", String::from_utf8_lossy(&out.stderr)); }
    let mut hits: BTreeMap<String, Vec<Hit>> = BTreeMap::new();
    for line in String::from_utf8(out.stdout)?.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 8 { bail!("unexpected blastn line: {line}"); }
        hits.entry(f[0].to_owned()).or_default().push(Hit {
            reference: f[1].to_owned(), pident: f[2].parse()?, length: f[3].parse()?, qlen: f[4].parse()?,
            slen: f[5].parse()?, evalue: f[6].parse()?, bitscore: f[7].parse()?,
        });
    }
    Ok(hits)
}

fn base_acc(acc: &str) -> &str { acc.split('.').next().unwrap_or("") }

/// arm_sets.py `vt_hits`.
pub fn reduce_hits(hits: &[Hit], refs: &BTreeMap<String, RefRecord>, min_cov: f64, max_ev: f64, min_id: f64,
                   drop: &StringSet) -> Result<(BTreeMap<String, f64>, Option<f64>)> {
    let mut ids: BTreeMap<String, f64> = BTreeMap::new();
    let mut best: Option<(f64, f64)> = None;
    for h in hits {
        let info = refs.get(&h.reference).ok_or_else(|| anyhow!("BLAST hit on unknown record {}", h.reference))?;
        if drop.contains(base_acc(&info.acc)) { continue; }
        let cov = h.length as f64 / (h.qlen.min(1.max(h.slen - info.n_n))) as f64;
        if cov < min_cov || h.evalue > max_ev { continue; }
        if h.pident > *ids.get(&info.label).unwrap_or(&-1.0) { ids.insert(info.label.clone(), h.pident); }
        if h.pident >= min_id && best.is_none_or(|b| (h.bitscore, h.pident) > b) { best = Some((h.bitscore, h.pident)); }
    }
    Ok((ids, best.map(|b| b.1)))
}

/// Accessions of the culture reference copies, every status, as the frozen evidence drops them.
pub fn culture_accessions(fixture_dir: &Path) -> Result<StringSet> {
    read_table(&fixture_dir.join("inputs/benchmark/ref/copies.tsv"))?.rows.iter()
        .map(|r| Ok(base_acc(row_value(r, "accession")?).to_owned())).collect()
}

#[derive(Clone, Debug)]
pub struct ArmEvidenceReport {
    pub arm: String,
    pub references: usize,
    pub frozen_references: usize,
    pub copies: usize,
    pub reproduced: usize,
    pub first_differences: Vec<String>,
}

#[derive(Clone, Debug)]
pub struct EvidenceReport { pub arms: Vec<ArmEvidenceReport> }

impl EvidenceReport {
    pub fn passed(&self) -> bool {
        !self.arms.is_empty() && self.arms.iter().all(|a| a.references == a.frozen_references && a.reproduced == a.copies)
    }
    pub fn summary(&self) -> String {
        let mut s = String::new();
        for a in &self.arms {
            s.push_str(&format!("{}: {} references kept (frozen {}); {}/{} copies reproduce the frozen evidence\n",
                                a.arm, a.references, a.frozen_references, a.reproduced, a.copies));
            for d in &a.first_differences { s.push_str("  "); s.push_str(d); s.push('\n'); }
        }
        let arms: Vec<&str> = self.arms.iter().map(|a| a.arm.as_str()).collect();
        let what = if arms.len() == 3 { "VT, SH and LSU".to_owned() } else { arms.join(", ") };
        s.push_str(&if self.passed() { format!("EVIDENCE ({what}, BLAST+ from the culture windows): PASS (exact frozen evidence)\n") }
                   else { "EVIDENCE: FAIL (see above)\n".to_owned() });
        s
    }
}

/// Rebuild the frozen evidence of the 903 culture copies from their windows and the three databases
/// (`databases`: VT MaarjAM, SH UNITE, LSU V18 AMF-only), and compare every query with the frozen
/// evidence files in the fixture.
/// `arms` limits the check to some arms (empty: all three); `sample` searches only every k-th retained
/// copy (k = ceil(903 / sample)). Each query's BLAST hits depend only on the query and the database,
/// so a sample checks exactly the queries it contains.
pub fn run_evidence_parity(fixture_dir: impl AsRef<Path>, databases: &BTreeMap<String, PathBuf>, workdir: &Path,
                           threads: usize, makeblastdb: &str, blastn: &str, arms_only: &[String], sample: Option<usize>) -> Result<EvidenceReport> {
    let fixture_dir = fixture_dir.as_ref();
    let bench = fixture_dir.join("inputs/benchmark");
    let drop = culture_accessions(fixture_dir)?;
    let retained: BTreeSet<String> = read_table(&bench.join("ref/copies.tsv"))?.rows.iter()
        .filter(|r| r.get("status").map(String::as_str) == Some("retained"))
        .map(|r| Ok(row_value(r, "accession")?.to_owned())).collect::<Result<_>>()?;
    let retained: BTreeSet<String> = match sample {
        Some(n) if n > 0 && n < retained.len() => {
            let k = retained.len().div_ceil(n);
            retained.into_iter().enumerate().filter(|(i, _)| i % k == 0).map(|(_, x)| x).collect()
        }
        _ => retained,
    };
    let mut arms = Vec::new();
    for (arm, kind, region, task, min_id, min_cov) in EVIDENCE_ARMS {
        if !arms_only.is_empty() && !arms_only.iter().any(|a| a.eq_ignore_ascii_case(arm)) { continue; }
        let path = databases.get(arm).ok_or_else(|| anyhow!("no database for the {arm} arm"))?;
        let (refs, fasta, _) = load_reference(path, kind, (kind == "unite").then_some("Glomeromycota"))?;
        // the region's FASTA records of the copies checked, copied verbatim
        let region_path = bench.join(format!("ref/fasta/{region}.fasta"));
        let query = workdir.join(format!("{}.queries.fasta", arm.to_lowercase()));
        fs::create_dir_all(workdir)?;
        let text = fs::read_to_string(&region_path)?;
        let mut kept = String::new();
        for rec in text.split('>').filter(|r| !r.trim().is_empty()) {
            let id = rec.split_whitespace().next().unwrap_or("");
            if retained.contains(id) { kept.push('>'); kept.push_str(rec); }
        }
        fs::write(&query, kept)?;
        let hits = run_blast(&query, &fasta, workdir, &arm.to_lowercase(), task, threads, makeblastdb, blastn)?;
        let frozen: Value = serde_json::from_str(&read_text_maybe_gzip(&bench.join(format!("evidence/{arm}.evidence.json.gz")))?)?;
        let queries = frozen.get("queries").and_then(Value::as_object).ok_or_else(|| anyhow!("{arm} evidence has no queries"))?;
        let (mut reproduced, mut first) = (0, Vec::new());
        for x in &retained {
            let (ids, anchor) = reduce_hits(hits.get(x).map(Vec::as_slice).unwrap_or(&[]), &refs, min_cov, MAX_EVALUE, min_id, &drop)?;
            let q = queries.get(x).ok_or_else(|| anyhow!("{arm} evidence has no {x}"))?;
            let f_ids: BTreeMap<String, f64> = q.get("unit_pident").and_then(Value::as_object).ok_or_else(|| anyhow!("{arm}/{x}: no unit_pident"))?
                .iter().map(|(u, p)| (u.clone(), p.as_f64().unwrap_or(f64::NAN))).collect();
            let f_anchor = q.get("anchor_pident").and_then(Value::as_f64);
            if ids == f_ids && anchor == f_anchor { reproduced += 1; }
            else if first.len() < 5 { first.push(format!("{x}: {} units (frozen {}), anchor {anchor:?} (frozen {f_anchor:?})", ids.len(), f_ids.len())); }
        }
        let frozen_references = frozen.get("reference_records").and_then(Value::as_u64).unwrap_or(0) as usize;
        arms.push(ArmEvidenceReport { arm: arm.into(), references: refs.len(), frozen_references, copies: retained.len(), reproduced, first_differences: first });
    }
    Ok(EvidenceReport { arms })
}
