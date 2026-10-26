//! Frozen culture-bridge naming and novelty decisions for AMF linked long reads.
//!
//! This crate ports the M2b baseline used by the reproducibility gate.  It consumes only the
//! checksum-pinned culture reference and numeric per-unit identity evidence; it does not need
//! to redistribute or download the external marker databases.  The public API is deliberately
//! narrow until the production reference-search interface is wired into `amfsplit name`.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, bail, Context, Result};
use flate2::read::GzDecoder;
use sha2::{Digest, Sha256};
use serde_json::Value;
use amfsplit_edlib::{align, Mode, Task};
use amfsplit_regions::pycompat::{csv_row, float_repr, py_round, read_tsv_rows};

pub type Row = BTreeMap<String, String>;

type SeqMap = BTreeMap<String, String>;
type StringSet = BTreeSet<String>;

#[derive(Clone, Debug)]
struct Table {
    headers: Vec<String>,
    rows: Vec<Row>,
}

fn read_table(path: &Path) -> Result<Table> {
    let text = fs::read_to_string(path).with_context(|| format!("reading TSV {}", path.display()))?;
    let parsed = read_tsv_rows(&text);
    if parsed.is_empty() {
        bail!("TSV has no header: {}", path.display());
    }
    let headers = parsed[0].clone();
    let mut rows = Vec::with_capacity(parsed.len().saturating_sub(1));
    for (i, raw) in parsed.into_iter().enumerate().skip(1) {
        if raw.iter().all(|s| s.is_empty()) { continue; }
        let mut row = Row::new();
        for (j, key) in headers.iter().enumerate() {
            row.insert(key.clone(), raw.get(j).cloned().unwrap_or_default());
        }
        if raw.len() > headers.len() {
            bail!("{} row {} has {} fields; header has {}", path.display(), i + 1, raw.len(), headers.len());
        }
        rows.push(row);
    }
    Ok(Table { headers, rows })
}

fn read_text_maybe_gzip(path: &Path) -> Result<String> {
    if path.extension().and_then(|x| x.to_str()) == Some("gz") {
        let f = fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
        let mut decoder = GzDecoder::new(f);
        let mut text = String::new();
        std::io::Read::read_to_string(&mut decoder, &mut text)
            .with_context(|| format!("decompressing {}", path.display()))?;
        Ok(text)
    } else {
        fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))
    }
}

fn read_fasta(path: &Path) -> Result<SeqMap> {
    let text = read_text_maybe_gzip(path)?;
    let mut out = SeqMap::new();
    let mut name: Option<String> = None;
    let mut seq = String::new();
    for line in text.lines() {
        let line = line.trim();
        if let Some(head) = line.strip_prefix('>') {
            if let Some(old) = name.take() { out.insert(old, std::mem::take(&mut seq)); }
            name = head.split_whitespace().next().map(str::to_owned);
            if name.is_none() { bail!("empty FASTA header in {}", path.display()); }
        } else if name.is_some() {
            seq.push_str(line);
        }
    }
    if let Some(old) = name { out.insert(old, seq); }
    for s in out.values_mut() { s.make_ascii_uppercase(); }
    Ok(out)
}

fn sha256_bytes(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn sha256_file(path: &Path) -> Result<String> {
    let bytes = fs::read(path).with_context(|| format!("reading for SHA-256: {}", path.display()))?;
    Ok(sha256_bytes(&bytes))
}

fn get_string<'a>(v: &'a Value, key: &str) -> Result<&'a str> {
    v.get(key).and_then(Value::as_str).ok_or_else(|| anyhow!("missing or non-string JSON field {key}"))
}

fn get_f64(v: &Value, key: &str) -> Result<f64> {
    v.get(key).and_then(Value::as_f64).ok_or_else(|| anyhow!("missing or non-numeric JSON field {key}"))
}

fn row_value<'a>(row: &'a Row, key: &str) -> Result<&'a str> {
    row.get(key).map(String::as_str).ok_or_else(|| anyhow!("missing TSV column {key}"))
}

fn parse_f64_field(row: &Row, key: &str) -> Result<Option<f64>> {
    let s = row_value(row, key)?;
    if s.is_empty() { Ok(None) }
    else { s.parse::<f64>().map(Some).with_context(|| format!("parsing {key}={s:?}")) }
}

fn py_float(x: f64) -> String { float_repr(x) }

fn py_round_float(x: f64, digits: usize) -> f64 { py_round(x, digits) }

/// Infix edit distance as used by the frozen Python reference: the shorter sequence is the query;
/// equal lengths preserve the order of the arguments (important for edlib HW asymmetry).
fn percent_distance(a: &str, b: &str, k: i32) -> Result<Option<f64>> {
    let (q, t) = if a.len() <= b.len() { (a.as_bytes(), b.as_bytes()) } else { (b.as_bytes(), a.as_bytes()) };
    if q.is_empty() { bail!("empty sequence in infix distance"); }
    let d = align(q, t, Mode::Hw, Task::Distance, k, &[])
        .map_err(|e| anyhow!("edlib distance failed: {e}"))?.edit_distance;
    if d < 0 { Ok(None) } else { Ok(Some(100.0 * f64::from(d) / q.len() as f64)) }
}

fn all_intersections_nonempty(sets: &[StringSet]) -> bool {
    let Some(first) = sets.first() else { return false; };
    first.iter().any(|u| sets.iter().skip(1).all(|s| s.contains(u)))
}

fn set_intersection(sets: &[StringSet]) -> StringSet {
    let Some(first) = sets.first() else { return StringSet::new(); };
    first.iter().filter(|u| sets.iter().skip(1).all(|s| s.contains(*u))).cloned().collect()
}

fn join_set(items: &StringSet, separator: &str) -> String { items.iter().cloned().collect::<Vec<_>>().join(separator) }

fn split_set(s: &str, sep: char) -> StringSet {
    s.split(sep).filter(|x| !x.is_empty()).map(str::to_owned).collect()
}

fn is_named(organism: &str) -> bool {
    if organism.is_empty() { return false; }
    // Equivalent to the frozen Python regex r"\bsp\.?(\s|$)" (case-sensitive).
    let bytes = organism.as_bytes();
    for i in 0..bytes.len() {
        if bytes[i..].starts_with(b"sp") && (i == 0 || !bytes[i - 1].is_ascii_alphanumeric() && bytes[i - 1] != b'_') {
            let mut j = i + 2;
            if j < bytes.len() && bytes[j] == b'.' { j += 1; }
            if j == bytes.len() || bytes[j].is_ascii_whitespace() { return false; }
        }
    }
    true
}

fn normalise_genus(g: &str) -> String {
    let g = g.trim_matches(|c| c == '[' || c == ']');
    match g {
        "Rhizoglomus" => "Rhizophagus".into(),
        "Claroideoglomus" => "Entrophospora".into(),
        _ => g.to_owned(),
    }
}

#[derive(Clone, Debug)]
struct CopyMeta {
    culture: String,
    organism: String,
    copy_class: String,
    div_kind: String,
    s58_len: Option<usize>,
}

#[derive(Clone, Debug)]
struct CultureMeta {
    label: String,
    organism: String,
    genus: String,
}

#[derive(Clone, Debug)]
struct EvidenceQuery {
    anchor_pident: Option<f64>,
    unit_pident: BTreeMap<String, f64>,
}

#[derive(Clone, Debug)]
struct ArmEvidence {
    region: String,
    min_id: f64,
    default_spread_quantile: usize,
    queries: BTreeMap<String, EvidenceQuery>,
}

#[derive(Clone, Debug)]
struct FrozenData {
    copies: BTreeMap<String, CopyMeta>,
    copy_order: Vec<String>,
    cultures: BTreeMap<String, CultureMeta>,
    by_culture: BTreeMap<String, Vec<String>>,
    by_label: BTreeMap<String, StringSet>,
    species: BTreeMap<String, String>,
    by_species: BTreeMap<String, StringSet>,
    genus: BTreeMap<String, String>,
    fa: BTreeMap<String, SeqMap>,
    evidence: BTreeMap<String, ArmEvidence>,
    best: BTreeMap<String, BTreeMap<String, (String, f64)>>,
    mains: Vec<String>,
    main58: BTreeMap<String, Vec<String>>,
    ref_fingerprint: String,
}

const ARM_NAMES: [&str; 3] = ["VT", "SH", "LSU"];

impl FrozenData {
    /// `inputs` is the frozen bundle's `inputs/benchmark` directory.
    pub fn load(inputs: impl AsRef<Path>) -> Result<Self> {
        let root = inputs.as_ref().to_path_buf();
        let ref_dir = root.join("ref");
        let base_dir = root.join("baseline");
        let evidence_dir = root.join("evidence");
        let fixture_root = root.parent().and_then(|p| p.parent())
            .ok_or_else(|| anyhow!("fixture inputs must be under parity/naming/inputs/benchmark"))?;
        let frozen_source_dir = fixture_root.join("frozen_sources");
        if !frozen_source_dir.is_dir() {
            bail!("frozen source bundle is missing: {}", frozen_source_dir.display());
        }
        let m2b_source_dir = frozen_source_dir.join("m2b");
        let source_lock_path = m2b_source_dir.join("SOURCE_SHA256SUMS");
        let source_lock_sha = sha256_file(&source_lock_path)?;
        if source_lock_sha != "c43dcac9caf0cf178734a15890e3014cfeac8c93c8792859a1c3d0de5069372f" {
            bail!("frozen M2b source checksum lock hash mismatch: {source_lock_sha}");
        }
        let source_lock = fs::read_to_string(&source_lock_path)
            .with_context(|| format!("reading frozen M2b source checksum lock {}", source_lock_path.display()))?;
        for line in source_lock.lines().filter(|line| !line.trim().is_empty()) {
            let mut fields = line.split_whitespace();
            let expected = fields.next().ok_or_else(|| anyhow!("bad source checksum line: {line}"))?;
            let name = fields.next().ok_or_else(|| anyhow!("bad source checksum line: {line}"))?;
            if fields.next().is_some() { bail!("bad source checksum line: {line}"); }
            let actual = sha256_file(&m2b_source_dir.join(name))?;
            if actual != expected { bail!("frozen M2b source hash mismatch for {name}: expected {expected}, found {actual}"); }
        }

        let ref_manifest_path = ref_dir.join("MANIFEST.json");
        let ref_manifest_sha = sha256_file(&ref_manifest_path)?;
        if ref_manifest_sha != "0c87060958bd1f1af58fc5c769b8b0ca2de861751660fb5ba32c8000d27d632a" {
            bail!("frozen culture-reference manifest hash mismatch: {ref_manifest_sha}");
        }
        let ref_manifest: Value = serde_json::from_slice(&fs::read(&ref_manifest_path)?)?;
        let ref_fingerprint = get_string(&ref_manifest, "fingerprint")?.to_owned();
        if ref_fingerprint != "181938d695f3df40" {
            bail!("unexpected culture-reference fingerprint {ref_fingerprint}");
        }
        if let Some(outputs) = ref_manifest.get("outputs").and_then(Value::as_object) {
            for (rel, expected) in outputs {
                let expected = expected.as_str().ok_or_else(|| anyhow!("non-string digest for reference output {rel}"))?;
                let file = ref_dir.join(rel);
                let actual = sha256_file(&file)?;
                if actual != expected { bail!("frozen reference hash mismatch for {}: expected {}, found {}", rel, expected, actual); }
            }
        } else { bail!("reference manifest has no outputs map"); }

        let copy_table = read_table(&ref_dir.join("copies.tsv"))?;
        let culture_table = read_table(&ref_dir.join("cultures.tsv"))?;
        let mut copies = BTreeMap::new();
        let mut copy_order = Vec::new();
        for r in &copy_table.rows {
            if row_value(r, "status")? != "retained" { continue; }
            let x = row_value(r, "accession")?.to_owned();
            let n = row_value(r, "s58_len")?.parse::<usize>().ok();
            copies.insert(x.clone(), CopyMeta {
                culture: row_value(r, "culture")?.to_owned(),
                organism: row_value(r, "organism")?.to_owned(),
                copy_class: row_value(r, "copy_class")?.to_owned(),
                div_kind: r.get("div_kind").cloned().unwrap_or_default(),
                s58_len: n,
            });
            copy_order.push(x);
        }
        let mut cultures = BTreeMap::new();
        let mut by_label: BTreeMap<String, StringSet> = BTreeMap::new();
        let mut by_species: BTreeMap<String, StringSet> = BTreeMap::new();
        let mut species = BTreeMap::new();
        let mut genus = BTreeMap::new();
        for r in &culture_table.rows {
            let c = row_value(r, "culture")?.to_owned();
            let label = row_value(r, "culture_label")?.to_owned();
            let organism = row_value(r, "organism")?.to_owned();
            let g = normalise_genus(row_value(r, "genus")?);
            let name_status = r.get("name_status").cloned().unwrap_or_else(|| "ok".into());
            if is_named(&organism) && name_status == "ok" {
                species.insert(c.clone(), organism.clone());
                by_species.entry(organism.clone()).or_default().insert(c.clone());
            }
            genus.insert(c.clone(), g.clone());
            by_label.entry(label.clone()).or_default().insert(c.clone());
            cultures.insert(c.clone(), CultureMeta { label, organism, genus: g });
        }
        let mut by_culture: BTreeMap<String, Vec<String>> = BTreeMap::new();
        // Preserve the source table's copy order within each culture. Some equal-length edlib HW
        // distances are asymmetric, so the `combinations()` order is part of the frozen behavior.
        for x in &copy_order {
            let c = copies.get(x).expect("copy_order is drawn from copies").culture.clone();
            by_culture.entry(c).or_default().push(x.clone());
        }
        let mut fa = BTreeMap::new();
        for marker in ["ssu_flank", "full_its", "s58", "LSU_LROR-FLR2"] {
            let p = ref_dir.join("fasta").join(format!("{marker}.fasta"));
            fa.insert(marker.to_owned(), read_fasta(&p)?);
        }

        let evidence_manifest_path = evidence_dir.join("MANIFEST.json");
        let evidence_manifest_sha = sha256_file(&evidence_manifest_path)?;
        if evidence_manifest_sha != "969fa89314d4099cc8739be71a985bdc64fc7a174364e70b3231090e8fceab09" {
            bail!("frozen evidence manifest hash mismatch: {evidence_manifest_sha}");
        }
        let evidence_manifest: Value = serde_json::from_slice(&fs::read(&evidence_manifest_path)?)?;
        if get_string(&evidence_manifest, "schema")? != "amf_phase4_evidence_bundle_v1" { bail!("unexpected evidence manifest schema"); }
        if get_string(&evidence_manifest, "reference_fingerprint")? != ref_fingerprint { bail!("evidence manifest refers to another culture-reference fingerprint"); }
        let entries = evidence_manifest.get("files").and_then(Value::as_array).ok_or_else(|| anyhow!("evidence manifest files missing"))?;
        let locked_evidence = [
            ("VT.evidence.json.gz", "5bcd06453d32f8718f573368009b12d19b55c58c2712e439408d26e1a9217f56"),
            ("SH.evidence.json.gz", "de21b3fa16fe8df6ca43421260acfd912391ef513e6bf5ead86853be4b9e81f6"),
            ("LSU.evidence.json.gz", "1a259082f41ef44cdfca2d02b64dc7ef17c34473b693db56a292f5d35ecc9c4b"),
        ];
        if entries.len() != locked_evidence.len() { bail!("evidence manifest does not contain exactly the three frozen culture evidence files"); }
        for (name, expected_sha) in locked_evidence {
            let p = evidence_dir.join(name);
            let bytes = fs::read(&p).with_context(|| format!("reading evidence file {}", p.display()))?;
            let item = entries.iter().find(|item| item.get("path").and_then(Value::as_str) == Some(name))
                .ok_or_else(|| anyhow!("evidence manifest omits {name}"))?;
            let want_bytes = item.get("bytes").and_then(Value::as_u64).ok_or_else(|| anyhow!("missing byte count for {name}"))?;
            let want_sha = get_string(item, "sha256")?;
            let actual_sha = sha256_bytes(&bytes);
            if bytes.len() as u64 != want_bytes || actual_sha != want_sha || actual_sha != expected_sha {
                bail!("evidence file checksum/size mismatch for {name}: got {actual_sha}");
            }
        }
        let crosswalk: Value = serde_json::from_slice(&fs::read(base_dir.join("crosswalk_v4.run.json"))?)?;
        for (name, expected) in [
            ("crosswalk_v4.run.json", "538ac476ee5b5a0ce7d374b374638633139e0f8b564739f41b8b16b188eb40d5"),
            ("vt_sets_v4.copies.tsv", "3cb6e82e100e19f9312ad748bdc6543cecd51c46e3f80c62d615c59262cc1d7a"),
            ("sets_SH_v4.copies.tsv", "9320f140fc45fb85911a75da2b2b5bd1c9beb46991fc55aebd4bcd75bedef657"),
            ("sets_LSU_v4.copies.tsv", "94b300b42aba7a2a14002832fa88cb0212d84d6987358155a033098b4cb1b517"),
        ] {
            let actual = sha256_file(&base_dir.join(name))?;
            if actual != expected { bail!("frozen baseline table hash mismatch for {name}: expected {expected}, found {actual}"); }
        }
        let mut evidence = BTreeMap::new();
        let mut best: BTreeMap<String, BTreeMap<String, (String, f64)>> = BTreeMap::new();
        let specs = [
            ("VT", "ssu_flank", "megablast", 97.0, 0.95, 95usize, 3081usize, "vt_sets_v4.copies.tsv"),
            ("SH", "full_its", "blastn", 98.5, 0.90, 0usize, 16800usize, "sets_SH_v4.copies.tsv"),
            ("LSU", "LSU_LROR-FLR2", "blastn", 97.0, 0.90, 0usize, 460usize, "sets_LSU_v4.copies.tsv"),
        ];
        for (arm, region, task, min_id, min_cov, spread_q, reference_records, legacy_name) in specs {
            let p = evidence_dir.join(format!("{arm}.evidence.json.gz"));
            let mut dec = GzDecoder::new(fs::File::open(&p).with_context(|| format!("opening {}", p.display()))?);
            let mut bytes = Vec::new();
            std::io::Read::read_to_end(&mut dec, &mut bytes)?;
            let ev: Value = serde_json::from_slice(&bytes).with_context(|| format!("parsing {}", p.display()))?;
            if get_string(&ev, "schema")? != "amf_arm_unit_evidence_v1" || get_string(&ev, "arm")? != arm { bail!("wrong schema/arm in {}", p.display()); }
            let source_hashes = ev.get("source_script_sha256").and_then(Value::as_object).ok_or_else(|| anyhow!("{arm} source script hashes missing"))?;
            for src in ["arms.py", "arm_sets.py"] {
                let want = source_hashes.get(src).and_then(Value::as_str).ok_or_else(|| anyhow!("{arm} evidence has no source hash for {src}"))?;
                let got = sha256_file(&frozen_source_dir.join(src))?;
                if got != want { bail!("frozen source hash mismatch for {src}: expected {want}, found {got}"); }
            }
            if get_string(&ev, "region")? != region { bail!("{arm} evidence region changed"); }
            if get_string(&ev, "task")? != task { bail!("{arm} evidence BLAST task changed"); }
            if (get_f64(&ev, "min_id")? - min_id).abs() > f64::EPSILON || (get_f64(&ev, "min_cov")? - min_cov).abs() > f64::EPSILON { bail!("{arm} evidence thresholds differ from the frozen specification"); }
            if ev.get("query_records").and_then(Value::as_u64) != Some(copies.len() as u64) { bail!("{arm} evidence query-record count changed"); }
            if ev.get("reference_records").and_then(Value::as_u64) != Some(reference_records as u64) { bail!("{arm} evidence reference-record count changed"); }
            let got_qsha = get_string(&ev, "query_sha256")?;
            let query_path = ref_dir.join("fasta").join(format!("{region}.fasta"));
            if got_qsha != sha256_file(&query_path)? { bail!("{arm} evidence query FASTA hash mismatch"); }
            if get_string(&ev, "copies_sha256")? != sha256_file(&ref_dir.join("copies.tsv"))? { bail!("{arm} evidence copies.tsv hash mismatch"); }
            if get_string(&ev, "cultures_sha256")? != sha256_file(&ref_dir.join("cultures.tsv"))? { bail!("{arm} evidence cultures.tsv hash mismatch"); }
            if get_string(&ev, "culture_reference_fingerprint")? != ref_fingerprint { bail!("{arm} evidence reference fingerprint mismatch"); }
            let arm_sha = crosswalk.get("arms").and_then(|x| x.get(arm)).and_then(|x| x.get("sha256")).and_then(Value::as_str).ok_or_else(|| anyhow!("crosswalk has no SHA for {arm}"))?;
            if get_string(&ev, "reference_sha256")? != arm_sha { bail!("{arm} evidence external-reference digest mismatch"); }
            if get_f64(&ev, "max_evalue")? != 1e-50 { bail!("{arm} evidence max e-value changed"); }
            if ev.get("default_spread_quantile").and_then(Value::as_u64) != Some(spread_q as u64) { bail!("{arm} spread quantile changed"); }
            let query_values = ev.get("queries").and_then(Value::as_object).ok_or_else(|| anyhow!("{arm} evidence has no queries map"))?;
            if query_values.len() != copies.len() { bail!("{arm} evidence has {} queries for {} retained copies", query_values.len(), copies.len()); }
            let mut queries = BTreeMap::new();
            for (x, q) in query_values {
                let anchor = match q.get("anchor_pident") { Some(Value::Number(n)) => n.as_f64(), Some(Value::Null) | None => None, _ => bail!("bad anchor_pident for {arm}/{x}") };
                let mut unit_pident = BTreeMap::new();
                let ids = q.get("unit_pident").and_then(Value::as_object).ok_or_else(|| anyhow!("{arm}/{x} has no unit_pident map"))?;
                for (u, v) in ids {
                    let p = v.as_f64().ok_or_else(|| anyhow!("non-number identity for {arm}/{x}/{u}"))?;
                    if !p.is_finite() || !(0.0..=100.0).contains(&p) { bail!("invalid identity for {arm}/{x}/{u}"); }
                    unit_pident.insert(u.clone(), p);
                }
                queries.insert(x.clone(), EvidenceQuery { anchor_pident: anchor, unit_pident });
            }
            if !query_values.keys().all(|x| copies.contains_key(x)) { bail!("{arm} evidence query identifiers differ from the retained reference"); }
            let mut arm_best = BTreeMap::new();
            let old_table = read_table(&base_dir.join(legacy_name))?;
            if old_table.rows.len() != copies.len() { bail!("legacy {arm} table copy count differs"); }
            let col = if old_table.headers.iter().any(|x| x == "top_units") { "top_units" } else { "top_vts" };
            for r in &old_table.rows {
                let x = row_value(r, "accession")?;
                let q = queries.get(x).ok_or_else(|| anyhow!("legacy {arm} accession {x} absent from evidence"))?;
                let got_anchor = parse_f64_field(r, "anchor_pident")?;
                if q.anchor_pident != got_anchor { bail!("{arm} anchor identity changed for {x}: {:?} vs {:?}", q.anchor_pident, got_anchor); }
                let top = row_value(r, col)?;
                if !top.is_empty() {
                    let (unit, _) = top.split(',').next().unwrap_or("").rsplit_once(':').ok_or_else(|| anyhow!("bad {arm} top unit for {x}: {top:?}"))?;
                    let max = q.unit_pident.values().copied().fold(f64::NEG_INFINITY, f64::max);
                    if q.unit_pident.get(unit).copied() != Some(max) { bail!("recorded best is not full-precision maximum: {arm} {x}"); }
                    arm_best.insert(x.to_owned(), (unit.to_owned(), max));
                } else {
                    if !q.unit_pident.is_empty() { bail!("missing legacy best despite evidence: {arm} {x}"); }
                    arm_best.insert(x.to_owned(), (String::new(), 0.0));
                }
            }
            best.insert(arm.to_owned(), arm_best);
            evidence.insert(arm.to_owned(), ArmEvidence { region: region.to_owned(), min_id, default_spread_quantile: spread_q, queries });
        }

        let mains: Vec<String> = copy_order.iter().filter(|x| copies.get(*x).is_some_and(|c| c.copy_class == "main") && fa["ssu_flank"].contains_key(*x)).cloned().collect();
        let mut main58: BTreeMap<String, Vec<String>> = BTreeMap::new();
        for x in &mains {
            if copies.get(x).and_then(|r| r.s58_len).is_some() && fa["s58"].contains_key(x) {
                main58.entry(copies[x].culture.clone()).or_default().push(x.clone());
            }
        }
        for (x, r) in &copies {
            if let (Some(n), Some(seq)) = (r.s58_len, fa["s58"].get(x)) {
                if n != seq.len() { bail!("5.8S length metadata mismatch for {x}: {n} vs {}", seq.len()); }
            }
        }
        Ok(Self { copies, copy_order, cultures, by_culture, by_label, species, by_species, genus, fa, evidence, best, mains, main58, ref_fingerprint })
    }

    fn truth_group(&self, c: &str, drop: &StringSet) -> String {
        if let Some(sp) = self.species.get(c) {
            if self.by_species.get(sp).is_some_and(|all| all.iter().any(|other| !drop.contains(other))) {
                "named, species elsewhere".into()
            } else { "named, species held out".into() }
        } else { "unnamed ('sp.')".into() }
    }

    fn folds(&self, mode: &str) -> Result<Vec<(StringSet, Vec<String>)>> {
        if mode != "physical_culture" && mode != "physical_species" { bail!("unsupported holdout mode {mode:?}"); }
        let mut grouped: BTreeMap<Vec<String>, Vec<String>> = BTreeMap::new();
        for c in self.cultures.keys() {
            let mut drop = if mode.ends_with("species") {
                self.species.get(c).and_then(|sp| self.by_species.get(sp)).cloned().unwrap_or_else(|| [c.clone()].into_iter().collect())
            } else { [c.clone()].into_iter().collect() };
            // The published physical folds also exclude all cultures sharing the same culture label.
            let labels: Vec<String> = drop.iter().filter_map(|cc| self.cultures.get(cc).map(|meta| meta.label.clone())).collect();
            if mode.starts_with("physical_") {
                for label in labels {
                    if let Some(members) = self.by_label.get(&label) { drop.extend(members.iter().cloned()); }
                }
            }
            let key: Vec<String> = drop.iter().cloned().collect();
            let ids = grouped.entry(key).or_default();
            if let Some(xs) = self.by_culture.get(c) { ids.extend(xs.iter().cloned()); }
        }
        Ok(grouped.into_iter().map(|(d, mut xs)| { xs.sort(); (d.into_iter().collect(), xs) }).collect())
    }
}

#[derive(Clone, Debug)]
struct MergeProposal { culture: String, units: Vec<String> }

#[derive(Clone, Debug)]
struct Fit {
    sets: BTreeMap<String, StringSet>,
    holders: BTreeMap<String, StringSet>,
    anchors: BTreeMap<String, String>,
    best: BTreeMap<String, (String, f64)>,
}

#[derive(Clone, Debug)]
struct Trainer {
    arm: String,
    spread_quantile: usize,
    main_of: BTreeMap<String, Vec<String>>,
    raw: Vec<BTreeMap<String, StringSet>>,
    spreads: BTreeMap<String, Vec<f64>>,
    proposals: Vec<MergeProposal>,
}

fn set_at_delta(q: &EvidenceQuery, delta: f64, min_id: f64) -> StringSet {
    let Some(anchor) = q.anchor_pident else { return StringSet::new(); };
    let threshold = min_id.max(anchor - delta) - 1e-9;
    q.unit_pident.iter().filter(|(_, p)| **p >= threshold).map(|(u, _)| u.clone()).collect()
}

fn quantiles_exclusive(values: &[f64], n: usize) -> Vec<f64> {
    let mut d = values.to_vec();
    d.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let len = d.len();
    if len < 2 { return vec![0.0; n.saturating_sub(1)]; }
    let m = len + 1;
    let mut out = Vec::with_capacity(n.saturating_sub(1));
    for i in 1..n {
        let product = i * m;
        let j = (product / n).clamp(1, len - 1);
        // Python clamps the interpolation index, but retains the original product when it
        // calculates delta. For small samples this intentionally extrapolates at the ends.
        let delta = product as i64 - (j * n) as i64;
        // delta is signed: the exclusive method extrapolates outside the observed range for
        // very small samples, so an unsigned remainder/subtraction would wrap or panic.
        // Keep the same operation order as statistics.quantiles() on Python floats.
        out.push((d[j - 1] * (n as i64 - delta) as f64 + d[j] * delta as f64) / n as f64);
    }
    out
}

fn median(values: &[f64]) -> Result<f64> {
    if values.is_empty() { bail!("median of empty data"); }
    let mut d = values.to_vec(); d.sort_by(|a,b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    Ok(if d.len()%2 == 1 { d[d.len()/2] } else { (d[d.len()/2-1] + d[d.len()/2])/2.0 })
}

impl Trainer {
    fn new(data: &FrozenData, arm: &str) -> Result<Self> {
        let ev = data.evidence.get(arm).ok_or_else(|| anyhow!("missing {arm} evidence"))?;
        let qseq = &data.fa[&ev.region];
        let mut main_of = BTreeMap::new();
        let mut spreads = BTreeMap::new();
        for (c, xs) in &data.by_culture {
            let main: Vec<String> = xs.iter().filter(|x| data.copies.get(*x).is_some_and(|r| r.copy_class == "main") && qseq.contains_key(*x)).cloned().collect();
            let mut local = Vec::new();
            for i in 0..main.len() { for j in i+1..main.len() {
                local.push(percent_distance(&qseq[&main[i]], &qseq[&main[j]], -1)?.ok_or_else(|| anyhow!("unbounded distance returned no result"))?);
            }}
            main_of.insert(c.clone(), main);
            spreads.insert(c.clone(), local);
        }
        let mut raw = vec![BTreeMap::new(); 31];
        for d in 0..=30 {
            let delta = d as f64 * 0.1;
            for (x, q) in &ev.queries { raw[d].insert(x.clone(), set_at_delta(q, delta, ev.min_id)); }
        }
        let ssu = &data.fa["ssu_flank"];
        let mut proposals = Vec::new();
        for c in data.cultures.keys() {
            let Some(xs) = main_of.get(c) else { continue; };
            // Preserve first-seen unit-set order exactly as Python's insertion-ordered dict does.
            // This matters because equal-length edlib HW distances can be asymmetric.
            let mut side: Vec<(StringSet, Vec<String>)> = Vec::new();
            for x in xs {
                let set = raw[0].get(x).cloned().unwrap_or_default();
                if set.is_empty() { continue; }
                if let Some(i) = side.iter().position(|(seen, _)| seen == &set) {
                    side[i].1.push(x.clone());
                } else {
                    side.push((set, vec![x.clone()]));
                }
            }
            if side.len() < 2 { continue; }
            let sets: Vec<StringSet> = side.iter().map(|(set, _)| set.clone()).collect();
            if all_intersections_nonempty(&sets) { continue; }
            let sides: Vec<&Vec<String>> = side.iter().map(|(_, members)| members).collect();
            let mut gap = 0.0f64;
            for i in 0..sides.len() { for j in i+1..sides.len() {
                let mut closest = f64::INFINITY;
                for p in sides[i] { for q in sides[j] {
                    let d = if let (Some(a), Some(b)) = (ssu.get(p), ssu.get(q)) {
                        percent_distance(a, b, -1)?.unwrap_or(100.0)
                    } else { 100.0 };
                    if d < closest { closest = d; }
                }}
                if closest > gap { gap = closest; }
            }}
            if gap <= 1.0 {
                let units: Vec<String> = side.iter().flat_map(|(set, _)| set.iter().cloned()).collect::<StringSet>().into_iter().collect();
                proposals.push(MergeProposal { culture: c.clone(), units });
            }
        }
        Ok(Self { arm: arm.to_owned(), spread_quantile: ev.default_spread_quantile, main_of, raw, spreads, proposals })
    }

    fn fit(&self, data: &FrozenData, training: &StringSet) -> Result<Fit> {
        if training.is_empty() || !training.iter().all(|c| data.cultures.contains_key(c)) { bail!("invalid training culture set"); }
        let mut parent: BTreeMap<String, String> = BTreeMap::new();
        fn find(parent: &mut BTreeMap<String,String>, x: &str) -> String {
            parent.entry(x.to_owned()).or_insert_with(|| x.to_owned());
            let p=parent[x].clone();
            if p==x { return x.to_owned(); }
            let root=find(parent,&p); parent.insert(x.to_owned(),root.clone()); root
        }
        let mut proposal_list: Vec<&MergeProposal> = self.proposals.iter().filter(|p| training.contains(&p.culture)).collect();
        proposal_list.sort_by(|a,b| a.culture.cmp(&b.culture));
        for p in &proposal_list {
            let mut units=p.units.clone(); units.sort(); units.dedup();
            if let Some(first)=units.first() { for u in &units { find(&mut parent,u); }
                for u in units.iter().skip(1) { let root0=find(&mut parent,first); let root1=find(&mut parent,u); if root0!=root1 { parent.insert(root1,root0); } }
            }
        }
        let all_units: Vec<String> = parent.keys().cloned().collect();
        let mut comps: BTreeMap<String,StringSet> = BTreeMap::new();
        for u in &all_units { let root=find(&mut parent,u); comps.entry(root).or_default().insert(u.clone()); }
        let mut unit_map = BTreeMap::new();
        for component in comps.values() {
            if component.len()>1 { let composite=join_set(component,"+"); for u in component { unit_map.insert(u.clone(),composite.clone()); } }
        }

        let mut spread: Vec<f64> = Vec::new();
        for c in training { if let Some(v)=self.spreads.get(c) { spread.extend(v.iter().copied()); } }
        let qs = quantiles_exclusive(&spread,100);
        let d_spread = if self.spread_quantile==0 { 0.0 } else {
            let q=*qs.get(self.spread_quantile-1).unwrap_or(&0.0);
            let rounded = py_round_float(q / 0.1, 6).ceil() * 0.1;
            py_round_float(rounded.min(3.0),3)
        };
        let mut calibration = Vec::with_capacity(31);
        let mut d_cons: Option<f64> = None;
        for d in 0..=30 {
            let mut tested=0usize; let mut consistent=0usize;
            for c in training {
                let xs=self.main_of.get(c).cloned().unwrap_or_default();
                let sets: Vec<StringSet> = xs.iter().filter_map(|x| self.raw[d].get(x)).filter(|s| !s.is_empty()).map(|s| s.iter().map(|u|unit_map.get(u).cloned().unwrap_or_else(||u.clone())).collect()).collect();
                if sets.len()>=2 { tested+=1; if !set_intersection(&sets).is_empty() { consistent+=1; } }
            }
            let share=if tested==0 { f64::NAN } else { py_round_float(consistent as f64/tested as f64,4) };
            let delta=d as f64*0.1;
            if d_cons.is_none() && share.is_finite() && share>=0.99 { d_cons=Some(delta); }
            calibration.push((delta,tested,consistent,share));
        }
        let d_cons=d_cons.unwrap_or(3.0);
        let delta=d_cons.max(d_spread);
        let d_idx=(delta*10.0).round() as usize;
        let mut sets: BTreeMap<String,StringSet> = BTreeMap::new();
        let chosen_raw=&self.raw[d_idx.min(30)];
        for x in data.copies.keys() {
            let mapped: StringSet=chosen_raw.get(x).into_iter().flat_map(|s|s.iter()).map(|u|unit_map.get(u).cloned().unwrap_or_else(||u.clone())).collect();
            sets.insert(x.clone(),mapped);
        }
        let mut culture_sets = BTreeMap::new(); let mut holders: BTreeMap<String,StringSet> = BTreeMap::new();
        for c in training {
            let xs=self.main_of.get(c).cloned().unwrap_or_default();
            let per_copy: Vec<StringSet>=xs.iter().filter_map(|x| sets.get(x)).filter(|s|!s.is_empty()).cloned().collect();
            let inter=set_intersection(&per_copy);
            let cu_set=if !inter.is_empty(){inter}else{per_copy.iter().flat_map(|s|s.iter().cloned()).collect()};
            for u in &cu_set { holders.entry(u.clone()).or_default().insert(c.clone()); }
            culture_sets.insert(c.clone(),cu_set);
        }
        let best=data.best.get(&self.arm).cloned().unwrap_or_default();
        let mut anchors=BTreeMap::new();
        for x in data.copies.keys() {
            let best_u=best.get(x).map(|x|x.0.as_str()).unwrap_or("");
            let anchor=sets.get(x).into_iter().flat_map(|s|s.iter()).find(|unit|unit.split('+').any(|part|part==best_u)).cloned().unwrap_or_default();
            anchors.insert(x.clone(),anchor);
        }
        Ok(Fit { sets, holders, anchors, best })
    }
}

/// Per-copy LSU placement units for the species-clade policy of the frozen compare_holdouts.py:
/// the copy's `credible_species_clade_units` when its LR0R-FLR2 window was eligible for the bridge
/// alignment, otherwise none.
fn species_clade_sets(data: &FrozenData, placements: &BTreeMap<String, Row>) -> Result<BTreeMap<String, StringSet>> {
    let mut out = BTreeMap::new();
    for x in data.copies.keys() {
        let p = placements.get(x).ok_or_else(|| anyhow!("missing placement annotation for {x}"))?;
        let units: StringSet = if row_value(p, "bridge_alignment_eligible")? == "1" {
            row_value(p, "credible_species_clade_units")?.split(',').filter(|u| !u.is_empty()).map(str::to_owned).collect()
        } else { StringSet::new() };
        out.insert(x.clone(), units);
    }
    Ok(out)
}

/// LSU arm built from placement units (`placement_fit`, frozen compare_holdouts.py lines 24-31):
/// no merges or margins; a training culture holds the units its main-type copies share or, when
/// they share none, every unit any of them has. Anchors are empty, as in the frozen code.
fn placement_fit(data: &FrozenData, training: &StringSet, psets: &BTreeMap<String, StringSet>) -> Fit {
    let mut holders: BTreeMap<String, StringSet> = BTreeMap::new();
    for c in training {
        let per_copy: Vec<StringSet> = data.by_culture.get(c).into_iter().flatten()
            .filter(|x| data.copies.get(*x).is_some_and(|m| m.copy_class == "main"))
            .filter_map(|x| psets.get(x)).filter(|u| !u.is_empty()).cloned().collect();
        let shared = set_intersection(&per_copy);
        let units: StringSet = if !shared.is_empty() { shared } else { per_copy.iter().flat_map(|u| u.iter().cloned()).collect() };
        for u in units { holders.entry(u).or_default().insert(c.clone()); }
    }
    Fit { sets: psets.clone(), holders, anchors: BTreeMap::new(), best: data.best.get("LSU").cloned().unwrap_or_default() }
}

fn reconcile(data: &FrozenData, x: &str, fits: &BTreeMap<String, Fit>, drop: &StringSet) -> Result<Row> {
    let meta=data.copies.get(x).ok_or_else(|| anyhow!("unknown accession {x}"))?;
    let c=&meta.culture;
    let mut per: BTreeMap<String,(String,StringSet,StringSet,StringSet)> = BTreeMap::new();
    for arm in ARM_NAMES {
        let fit=fits.get(arm).ok_or_else(|| anyhow!("missing fit for arm {arm}"))?;
        let units=fit.sets.get(x).cloned().unwrap_or_default();
        let mut others=StringSet::new();
        for u in &units { if let Some(cultures)=fit.holders.get(u) { others.extend(cultures.iter().filter(|cu|!drop.contains(*cu)).cloned()); } }
        let sp: StringSet=others.iter().filter_map(|o| data.species.get(o).cloned()).collect();
        let gn: StringSet=others.iter().filter_map(|o| data.genus.get(o).cloned()).filter(|g|!g.is_empty()).collect();
        let state=(if units.is_empty(){"no call"}else if !gn.is_empty(){"matched"}else{"unmatched"}).to_owned();
        per.insert(arm.to_owned(),(state,sp,gn,units));
    }
    // Build intersections in the frozen VT, SH, LSU order.
    let matched_names: Vec<String> = ARM_NAMES.iter().copied().filter(|n| per.get(*n).is_some_and(|row| row.0 == "matched")).map(str::to_owned).collect();
    let (rank,call_sp,call_gn)=if matched_names.is_empty() {
        ("unplaced".to_owned(),StringSet::new(),StringSet::new())
    } else {
        let gn_sets: Vec<StringSet>=matched_names.iter().map(|n|per[n].2.clone()).collect();
        let call_gn=set_intersection(&gn_sets);
        let sp_arms: Vec<String> = matched_names.iter().filter(|n| per.get(n.as_str()).is_some_and(|row| !row.1.is_empty())).cloned().collect();
        let mut call_sp=if sp_arms.is_empty(){StringSet::new()}else{set_intersection(&sp_arms.iter().map(|n|per[n].1.clone()).collect::<Vec<_>>())};
        if !call_gn.is_empty() {
            call_sp.retain(|sp| data.by_species.get(sp).is_some_and(|holders| holders.iter().any(|o| !drop.contains(o) && call_gn.contains(data.genus.get(o).map(String::as_str).unwrap_or("")))));
        } else { call_sp.clear(); }
        let rank=if call_gn.is_empty(){"discordant"}else if call_sp.len()==1{"species"}else if !call_sp.is_empty(){"species group"}else{"genus"};
        (rank.to_owned(),call_sp,call_gn)
    };
    let true_sp=data.species.get(c).cloned().unwrap_or_default();
    let true_gn=data.genus.get(c).cloned().unwrap_or_default();
    let correct=if rank=="species" && call_sp.len()==1 && call_sp.contains(&true_sp) { "yes" }
      else if rank=="species group" && call_sp.contains(&true_sp) { "contains" }
      else if rank=="genus" && call_gn.contains(&true_gn) { "genus" }
      else if (rank=="species" || rank=="species group") && call_gn.contains(&true_gn) { "genus only" }
      else if rank=="species" || rank=="species group" || rank=="genus" { "no" }
      else { "" };
    let mut r=Row::new();
    r.insert("accession".into(),x.into()); r.insert("culture".into(),c.clone()); r.insert("organism".into(),meta.organism.clone());
    r.insert("genus".into(),true_gn); r.insert("copy_class".into(),meta.copy_class.clone()); r.insert("div_kind".into(),meta.div_kind.clone());
    r.insert("group".into(),data.truth_group(c,drop));
    let mut unmatched=Vec::new();
    for arm in ARM_NAMES {
        let (state,sp,gn,units)=&per[arm];
        r.insert(format!("{arm}_state"),state.clone());
        r.insert(format!("{arm}_species"),join_set(sp,";"));
        r.insert(format!("{arm}_genera"),join_set(gn,";"));
        r.insert(format!("{arm}_units"),join_set(units,","));
        if state=="unmatched" { unmatched.push(arm); }
    }
    r.insert("rank".into(),rank); r.insert("call_species".into(),join_set(&call_sp,";")); r.insert("call_genera".into(),join_set(&call_gn,";"));
    r.insert("correct".into(),correct.into()); r.insert("unmatched_arms".into(),unmatched.join(";"));
    Ok(r)
}

fn format_g(x: f64, precision: usize) -> String {
    if x == 0.0 { return "0".into(); }
    if !x.is_finite() { return if x.is_nan(){"nan".into()}else if x>0.0{"inf".into()}else{"-inf".into()}; }
    let exponent=x.abs().log10().floor() as i32;
    if exponent < -4 || exponent >= precision as i32 {
        let scientific=format!("{:.*e}",precision.saturating_sub(1),x);
        if let Some((mant,exp))=scientific.split_once('e') {
            let mut m=mant.trim_end_matches('0').trim_end_matches('.').to_owned();
            if m.is_empty() {m="0".into();}
            let e=exp.parse::<i32>().unwrap_or(0);
            format!("{}e{}{:02}",m,if e<0{'-'}else{'+'},e.abs())
        } else {
            scientific
        }
    } else {
        let decimals=(precision as i32-1-exponent).max(0) as usize;
        let s=format!("{:.*}",decimals,x);
        if s.contains('.') { s.trim_end_matches('0').trim_end_matches('.').to_owned() } else { s }
    }
}

fn format_g0(x: f64) -> String { format_g(x, 6) }

fn nearest(data: &FrozenData, x: &str, drop: &StringSet, vt_fit: &Fit) -> Result<(Option<String>,Option<f64>,Vec<String>)> {
    let Some(q)=data.fa["ssu_flank"].get(x) else { return Ok((None,None,Vec::new())); };
    let x_best=data.best["VT"].get(x).map(|t|t.0.as_str()).unwrap_or("");
    let mut candidates: Vec<String>=data.mains.iter().filter(|y|!drop.contains(&data.copies[*y].culture)).cloned().collect();
    candidates.sort_by(|a,b| {
        let ab=data.best["VT"].get(a).map(|t|t.0.as_str()).unwrap_or("") != x_best;
        let bb=data.best["VT"].get(b).map(|t|t.0.as_str()).unwrap_or("") != x_best;
        ab.cmp(&bb).then_with(||data.copies[a].culture.cmp(&data.copies[b].culture)).then_with(||a.cmp(b))
    });
    let mut near_d: Option<f64>=None; let mut near: Vec<String>=Vec::new();
    for y in candidates {
        let target=&data.fa["ssu_flank"][&y];
        let k=near_d.map_or(-1,|d|(d * q.len().min(target.len()) as f64 /100.0) as i32 +1);
        let Some(d)=percent_distance(q,target,k)? else { continue; };
        if near_d.is_none() || d < near_d.unwrap()-1e-12 { near_d=Some(d); near=vec![y]; }
        else if (d-near_d.unwrap()).abs()<=1e-12 { near.push(y); }
    }
    near.sort_by(|a,b|data.copies[a].culture.cmp(&data.copies[b].culture).then_with(||a.cmp(b)));
    let representative=near.first().cloned();
    let _=vt_fit;
    Ok((representative,near_d,near))
}

fn evaluate_novelty(data: &FrozenData, query_ids: &[String], fits: &BTreeMap<String,Fit>, drop: &StringSet, rec: &BTreeMap<String,Row>) -> Result<Vec<Row>> {
    let vt_fit=&fits["VT"];
    let vt_units=&vt_fit.sets; let vt_best=&vt_fit.best; let vt_anchor=&vt_fit.anchors; let vt_holders=&vt_fit.holders;
    let para_arms=["VT","SH"];
    let mut main_reach: BTreeMap<String,BTreeMap<String,StringSet>>=BTreeMap::new();
    let mut div_reach: BTreeMap<String,BTreeMap<String,StringSet>>=BTreeMap::new();
    for arm in para_arms {
        for (y,u) in &fits[arm].anchors {
            if u.is_empty() || drop.contains(&data.copies[y].culture) { continue; }
            let meta=&data.copies[y];
            if meta.copy_class=="main" { main_reach.entry(arm.into()).or_default().entry(u.clone()).or_default().insert(meta.culture.clone()); }
            else if meta.div_kind=="short_5.8S" { div_reach.entry(arm.into()).or_default().entry(u.clone()).or_default().insert(meta.culture.clone()); }
        }
    }
    let mut out=Vec::with_capacity(query_ids.len());
    for x in query_ids {
        let r=data.copies.get(x).ok_or_else(||anyhow!("unknown copy {x}"))?; let c=&r.culture;
        let mut reasons:Vec<String>=Vec::new(); let mut flags:Vec<String>=Vec::new();
        let (near,near_d,nearest_copies)=nearest(data,x,drop,vt_fit)?;
        let its_d=if let (Some(y),Some(q),Some(t))=(near.as_ref(),data.fa["full_its"].get(x),near.as_ref().and_then(|id|data.fa["full_its"].get(id))) { percent_distance(q,t,-1)? } else {None};

        let mut artifact=false; let mut artifact_its:Vec<f64>=Vec::new();
        if near_d.is_some_and(|d|d>=1.0) {
            let d=near_d.unwrap();
            artifact_its=nearest_copies.iter().filter_map(|y|Some((data.fa["full_its"].get(x)?,data.fa["full_its"].get(y)?))).map(|(q,t)|percent_distance(q,t,-1)).collect::<Result<Vec<_>>>()?.into_iter().flatten().collect();
            let votes:Vec<bool>=artifact_its.iter().map(|its|d>=0.5*its).collect();
            artifact=votes.len()==nearest_copies.len() && !votes.is_empty() && votes.iter().all(|v|*v);
            if votes.iter().any(|v|*v) && (!votes.iter().all(|v|*v) || votes.len()!=nearest_copies.len()) { flags.push("A2 artifact screen unresolved across SSU-tied nearest copies".into()); }
        }
        if artifact {
            let d=near_d.unwrap(); let min=artifact_its.iter().copied().fold(f64::INFINITY,f64::min); let max=artifact_its.iter().copied().fold(f64::NEG_INFINITY,f64::max);
            reasons.push(format!("A1 screening flag: SSU {:.1}% vs ITS {:.1}-{:.1}% across {} SSU-tied nearest copies; not proof of an artifact",d,min,max,nearest_copies.len()));
        }

        let mut short=false; let mut paralog=false;
        let ncul=near.as_ref().map(|y|data.copies[y].culture.clone()).unwrap_or_default();
        let mut comp:StringSet=vt_holders.get(vt_anchor.get(x).map(String::as_str).unwrap_or("")).cloned().unwrap_or_default();
        comp.retain(|cu|!drop.contains(cu) && data.main58.get(cu).is_some_and(|v|!v.is_empty()));
        let mut comp_src="cultures sharing its VT".to_owned();
        if comp.is_empty() && near_d.is_some_and(|d|d<=1.0) {
            comp=nearest_copies.iter().map(|y|data.copies[y].culture.clone()).filter(|cu|data.main58.get(cu).is_some_and(|v|!v.is_empty())).collect();
            comp_src=format!("SSU-tied nearest cultures ({:.2}% SSU)",near_d.unwrap());
        }
        let mut s58_info=String::new(); let mut compare_ids:Vec<String>=Vec::new();
        if r.s58_len.is_some() && !comp.is_empty() && data.fa["s58"].contains_key(x) {
            let mut culture_medians=Vec::new();
            for cu in &comp {
                let values:Vec<f64>=data.main58[cu].iter().filter_map(|y|data.copies[y].s58_len.map(|n|n as f64)).collect();
                culture_medians.push(median(&values)?);
            }
            let ref_len=median(&culture_medians)?;
            let delta=r.s58_len.unwrap() as f64-ref_len;
            for cu in &comp { compare_ids.extend(data.main58[cu].iter().cloned()); }
            if delta<=-3.0 { compare_ids.retain(|y|data.fa["s58"][y].len() as f64>=ref_len); }
            let mut d58=f64::INFINITY;
            for y in &compare_ids {
                let d=percent_distance(&data.fa["s58"][x],&data.fa["s58"][y],-1)?.ok_or_else(||anyhow!("unbounded 5.8S distance failed"))?;
                if d<d58 { d58=d; }
            }
            s58_info=format!("5.8S {} bp vs {} in {}, {:.1}% different",r.s58_len.unwrap(),format_g0(ref_len),comp_src,d58);
            if delta<=-3.0 && d58>=3.0 { short=true; reasons.push(format!("D1 {s58_info}")); }
            else if delta<=-3.0 { flags.push(format!("D4 {s58_info} (shorter, below divergence threshold: inspect boundaries)")); }
            else if delta.abs()<3.0 && d58>=3.0 { flags.push(format!("D3 {s58_info}")); }
        }
        for arm in para_arms {
            let u=fits[arm].anchors.get(x).map(String::as_str).unwrap_or("");
            if u.is_empty(){continue;}
            let m=main_reach.get(arm).and_then(|m|m.get(u)).cloned().unwrap_or_default();
            let d=div_reach.get(arm).and_then(|d|d.get(u)).cloned().unwrap_or_default();
            let m_count=m.iter().filter(|cu|!drop.contains(*cu)).count();
            let d_count=d.iter().filter(|cu|!drop.contains(*cu)).count();
            if d_count>=2 && m_count==0 { paralog=true; reasons.push(format!("D2 paralog unit {u} ({arm}; short-5.8S copies of {d_count} other cultures)")); }
        }
        let divergent=short||paralog;
        let rc=rec.get(x).ok_or_else(||anyhow!("no reconciliation row for {x}"))?;
        let vt_state=row_value(rc,"VT_state")?;
        let (best_vt,best_id)=vt_best.get(x).cloned().unwrap_or_default();
        let boundary=!best_vt.is_empty() && vt_state=="no call" && best_id>=96.0 && best_id<97.0;
        let off_vt=!best_vt.is_empty() && vt_state=="no call" && best_id<96.0;
        if vt_state=="no call" && !best_vt.is_empty() && best_id>=97.0 { flags.push("B2 VT identity passes threshold but no call: inspect coverage, filtering and provenance".into()); }
        if vt_state=="no call" && best_vt.is_empty() { flags.push("Q1 no recorded best VT hit; distance from all VTs is not established".into()); }
        let rank=row_value(rc,"rank")?;
        let um=split_set(row_value(rc,"unmatched_arms")?,';');
        let status=if artifact {"artifact".to_owned()}
          else if divergent {
            let place = if !row_value(rc,"VT_species")?.is_empty() { row_value(rc,"VT_species")? } else { row_value(rc,"VT_genera")? };
            reasons.push(format!("placed by VT: {}",if place.is_empty(){"no VT match"}else{place}));
            "divergent class".to_owned()
        }
          else if boundary { reasons.push(format!("B1 best VT {best_vt} {:.2}%",best_id)); "boundary".into() }
          else if rank=="species" || rank=="species group" {
            if um.contains("VT") || um.contains("LSU") {
                for arm in ["VT","LSU"] { if um.contains(arm) { reasons.push(format!("{} {arm} unmatched",if arm=="VT"{"N1"}else{"N2"})); } }
                "novel species, genus known".into()
            } else { reasons.push("K1".into()); format!("known {rank}") }
          }
          else if rank=="genus" { reasons.push("N3 no species shared by the arms".into()); "novel species, genus known".into() }
          else if rank=="discordant" { reasons.push("X1 arms share no genus".into()); "discordant".into() }
          else { reasons.push((if off_vt {"N5 recorded best VT identity below boundary"} else {"N4 no arm matched"}).into()); "novel lineage".into() };
        let placement=if status.starts_with("known") { row_value(rc,"call_species")?.to_owned() }
            else if status=="divergent class" { let s=row_value(rc,"VT_species")?; if !s.is_empty(){s.to_owned()}else{row_value(rc,"VT_genera")?.to_owned()} }
            else if status=="novel species, genus known" { row_value(rc,"call_genera")?.to_owned() }
            else { String::new() };
        let nearest_its_cultures: String = nearest_copies.iter().map(|y|data.copies[y].culture.clone()).collect::<StringSet>().iter().cloned().collect::<Vec<_>>().join(";");
        let artifact_comparator= {
            let mut d=artifact_its.clone(); d.sort_by(|a,b|a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            d.iter().map(|x|format_g(*x,8)).collect::<Vec<_>>().join(";")
        };
        let mut row=Row::new();
        row.insert("accession".into(),x.clone()); row.insert("culture".into(),c.clone()); row.insert("organism".into(),r.organism.clone()); row.insert("copy_class".into(),r.copy_class.clone()); row.insert("div_kind".into(),r.div_kind.clone()); row.insert("group".into(),row_value(rc,"group")?.to_owned());
        row.insert("status".into(),status.clone());
        row.insert("novelty_claim".into(),(if status.starts_with("novel") {"candidate_only"} else {"not_asserted"}).into());
        row.insert("artifact_claim".into(),(if status=="artifact" {"screen_flag_only"} else {"not_asserted"}).into());
        row.insert("placement".into(),placement);
        row.insert("reasons".into(),reasons.join("; ")); row.insert("flags".into(),flags.join("; ")); row.insert("s58".into(),s58_info);
        row.insert("nearest_culture".into(),ncul);
        row.insert("nearest_ssu_pct".into(),near_d.map(|d|py_float(py_round_float(d,2))).unwrap_or_default());
        row.insert("nearest_its_pct".into(),its_d.map(|d|py_float(py_round_float(d,2))).unwrap_or_default());
        row.insert("nearest_ssu_copy_count".into(),nearest_copies.len().to_string()); row.insert("nearest_ssu_cultures".into(),nearest_its_cultures); row.insert("artifact_comparator_its_pct".into(),artifact_comparator);
        row.insert("best_vt".into(),best_vt.clone()); row.insert("best_vt_pident".into(),py_float(best_id)); row.insert("reconcile_rank".into(),rank.into()); row.insert("unmatched_arms".into(),row_value(rc,"unmatched_arms")?.into());
        row.insert("sequence_type_signal".into(),(if short {"short_58S"} else if paralog {"paralog_unit"} else {"unresolved"}).into());
        row.insert("d1_comparator_cultures".into(),join_set(&comp,";"));
        let show_compare= r.s58_len.is_some() && !comp.is_empty() && data.fa["s58"].contains_key(x);
        row.insert("d1_comparator_copies".into(),if show_compare {compare_ids.join(";")} else {String::new()});
        row.insert("d1_short".into(),(if short {"1"} else {"0"}).into());
        row.insert("d2_paralog".into(),(if paralog {"1"} else {"0"}).into());
        row.insert("artifact_screen".into(),(if artifact {"1"} else {"0"}).into());
        row.insert("placement_rank_before_type_filter".into(),rank.into()); row.insert("placement_species_before_type_filter".into(),row_value(rc,"call_species")?.into());
        out.push(row);
    }
    Ok(out)
}

fn read_row_map(table: &Table) -> BTreeMap<String, Row> {
    table.rows.iter().filter_map(|r|r.get("accession").map(|x|(x.clone(),r.clone()))).collect()
}

fn serialize_table(headers: &[String], rows: &[Row]) -> Vec<u8> {
    let mut out=String::new(); out.push_str(&csv_row(headers));
    for r in rows {
        let fields:Vec<String>=headers.iter().map(|h|r.get(h).cloned().unwrap_or_default()).collect();
        out.push_str(&csv_row(&fields));
    }
    out.into_bytes()
}

#[derive(Clone, Debug)]
pub struct ModeReport {
    pub mode: String,
    pub folds: usize,
    pub rows: usize,
    pub sha256: String,
    pub byte_identical: bool,
    pub differing_cells: usize,
    pub first_differences: Vec<String>,
}

#[derive(Clone, Debug)]
pub struct Gate0Report {
    pub policy: String,
    pub reference_fingerprint: String,
    pub evidence_hashes: BTreeMap<String,String>,
    pub modes: Vec<ModeReport>,
}

impl Gate0Report {
    pub fn passed(&self) -> bool { self.modes.iter().all(|m|m.byte_identical && m.differing_cells==0) }
    pub fn summary(&self) -> String {
        let mut s=format!("policy {}; reference fingerprint {}\n",self.policy,self.reference_fingerprint);
        for m in &self.modes {
            s.push_str(&format!("{}: {} holdout rounds; {}/903 rows; byte-identical={}; generated SHA-256 {}; differing cells={}\n",m.mode,m.folds,m.rows,m.byte_identical,m.sha256,m.differing_cells));
            for d in &m.first_differences { s.push_str("  "); s.push_str(d); s.push('\n'); }
        }
        let gate = match self.policy.as_str() {
            "baseline" => "GATE 0a, nearest-reference LSU arm",
            "epa_species_clades" => "GATE 0b, LSU species-clade arm",
            _ => "GATE 0",
        };
        s.push_str(&if self.passed(){format!("{gate}: PASS (exact archived per-copy output parity)\n")}else{format!("{gate}: FAIL (see per-copy differences above)\n")});
        s
    }
}

/// Execute both frozen culture-only holdouts against the saved M2b assignment tables.
/// `fixture_dir` is `parity/naming`; its inputs are checksum-pinned and exclude external raw reads.
pub fn run_gate0(fixture_dir: impl AsRef<Path>) -> Result<Gate0Report> {
    run_policy(fixture_dir.as_ref(), "baseline")
}

/// Gate 0b, step 1: the same holdouts with the LSU arm built from placement species-clade units,
/// against the archived `epa_species_clades_*` tables of the frozen M2b comparison.
pub fn run_gate0b_species_clades(fixture_dir: impl AsRef<Path>) -> Result<Gate0Report> {
    run_policy(fixture_dir.as_ref(), "epa_species_clades")
}

fn run_policy(fixture_dir: &Path, policy: &str) -> Result<Gate0Report> {
    if !matches!(policy, "baseline" | "epa_species_clades") { bail!("unknown naming policy {policy}"); }
    let input_dir=fixture_dir.join("inputs/benchmark");
    let data=FrozenData::load(&input_dir)?;
    let trainers:BTreeMap<String,Trainer>=ARM_NAMES.iter().map(|arm|Ok(((*arm).to_owned(),Trainer::new(&data,arm)?))).collect::<Result<_>>()?;
    let placement_path=fixture_dir.join("placement/lsu_placements.copies.tsv");
    let placement_sha=sha256_file(&placement_path)?;
    if placement_sha != "1de103b6ae45f0372dda9ab0799230c4d1243cb98add943e5bab67785e5e910d" {
        bail!("frozen LSU placement annotation hash mismatch: {placement_sha}");
    }
    let placement_table=read_table(&placement_path)?;
    let placements=read_row_map(&placement_table);
    if placements.len()!=data.copies.len() || placements.keys().any(|x|!data.copies.contains_key(x)) { bail!("LSU placement IDs differ from the 903-copy culture reference"); }
    let psets = if policy == "epa_species_clades" { Some(species_clade_sets(&data, &placements)?) } else { None };
    let mut modes=Vec::new();
    for mode in ["physical_culture","physical_species"] {
        let mut results=Vec::with_capacity(data.copies.len());
        let folds=data.folds(mode)?;
        let fold_count=folds.len();
        for (drop,query_ids) in &folds {
            let training:StringSet=data.cultures.keys().filter(|c|!drop.contains(*c)).cloned().collect();
            let mut fits:BTreeMap<String,Fit>=ARM_NAMES.iter().map(|arm|Ok(((*arm).to_owned(),trainers[*arm].fit(&data,&training)?))).collect::<Result<_>>()?;
            if let Some(ps)=&psets { fits.insert("LSU".into(),placement_fit(&data,&training,ps)); }
            if fits["LSU"].holders.values().any(|h|h.iter().any(|c|drop.contains(c))) { bail!("held-out culture leaked into LSU training holders in {mode}"); }
            let mut reconciliation=BTreeMap::new();
            for x in query_ids { reconciliation.insert(x.clone(),reconcile(&data,x,&fits,drop)?); }
            let novelties=evaluate_novelty(&data,query_ids,&fits,drop,&reconciliation)?;
            for novelty in novelties {
                let x=row_value(&novelty,"accession")?.to_owned();
                let mut row=reconciliation.get(&x).cloned().ok_or_else(||anyhow!("missing reconciliation for {x}"))?;
                row.extend(novelty);
                row.insert("policy".into(),policy.into());
                row.insert("mode".into(),mode.into());
                row.insert("held_cultures".into(),join_set(drop,";"));
                let p=placements.get(&x).ok_or_else(||anyhow!("missing placement annotation for {x}"))?;
                for k in ["best_edge_lwr","credible_edge_count","best_pendant_length","retained_fraction","external_species_clade","external_species_clade_mass","external_genus_clade","external_genus_clade_mass"] {
                    row.insert(format!("lsu_{k}"),row_value(p,k)?.to_owned());
                }
                results.push(row);
            }
        }
        results.sort_by(|a,b|a.get("accession").cmp(&b.get("accession")));
        if results.len()!=data.copies.len() { bail!("{mode} produced {} per-copy results, expected {}",results.len(),data.copies.len()); }
        let expected_path=fixture_dir.join("golden").join(format!("{policy}_{mode}.copies.tsv"));
        let expected_sha = match (policy, mode) {
            ("baseline", "physical_culture") => "6cc5412acf147ba7644d767851872517f78551b372e614c3f07297227d539f34",
            ("baseline", "physical_species") => "906420d70e9d8d24ab5b115ecd63e72451c396b4f1e7c60cfccf7b8d6a8a0800",
            ("epa_species_clades", "physical_culture") => "a357b563a5def31932085dabe4697f3b49e0f586d56a6c91acefe8bb1cb7cfd4",
            ("epa_species_clades", "physical_species") => "9b0d4ba422214a763cbf0c5949f8f7b2b163fdfcd8e2fd5ebf81a55db521a53a",
            _ => unreachable!("policy and mode are checked above"),
        };
        let actual_golden_sha = sha256_file(&expected_path)?;
        if actual_golden_sha != expected_sha { bail!("golden output hash mismatch for {mode}: {actual_golden_sha}"); }
        let expected_bytes=fs::read(&expected_path)?;
        let expected=read_table(&expected_path)?;
        let expected_by_id=read_row_map(&expected);
        let actual_by_id=read_row_map(&Table{headers:expected.headers.clone(),rows:results.clone()});
        let mut diffs=Vec::new(); let mut differing_cells=0usize;
        if expected_by_id.len()!=actual_by_id.len() { differing_cells+=expected_by_id.len().abs_diff(actual_by_id.len()); }
        for (x, exp) in &expected_by_id {
            let Some(act)=actual_by_id.get(x) else { differing_cells+=1; if diffs.len()<12 {diffs.push(format!("{x}: missing output row"));} continue; };
            for h in &expected.headers {
                let a=act.get(h).map(String::as_str).unwrap_or(""); let e=exp.get(h).map(String::as_str).unwrap_or("");
                if a!=e { differing_cells+=1; if diffs.len()<12 { diffs.push(format!("{x} [{h}]: got {:?}; expected {:?}",a,e)); } }
            }
        }
        let serialized=serialize_table(&expected.headers,&results);
        let byte_identical=serialized==expected_bytes;
        if !byte_identical && differing_cells==0 && diffs.len()<12 { diffs.push("per-cell output matched but TSV serialization did not (header/order/escaping difference)".into()); }
        modes.push(ModeReport {mode:mode.into(),folds:fold_count,rows:results.len(),sha256:sha256_bytes(&serialized),byte_identical,differing_cells,first_differences:diffs});
    }
    let evidence_dir=input_dir.join("evidence");
    let evidence_hashes=["MANIFEST.json","VT.evidence.json.gz","SH.evidence.json.gz","LSU.evidence.json.gz"].iter().map(|n|Ok(((*n).to_owned(),sha256_file(&evidence_dir.join(n))?))).collect::<Result<BTreeMap<_,_>>>()?;
    Ok(Gate0Report { policy:policy.into(),reference_fingerprint:data.ref_fingerprint,evidence_hashes,modes })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../parity/naming").canonicalize().unwrap()
    }

    #[test]
    fn frozen_evidence_hashes_match_the_lock() {
        let f=fixture().join("inputs/benchmark/evidence");
        assert_eq!(sha256_file(&f.join("MANIFEST.json")).unwrap(),"969fa89314d4099cc8739be71a985bdc64fc7a174364e70b3231090e8fceab09");
        assert_eq!(sha256_file(&f.join("VT.evidence.json.gz")).unwrap(),"5bcd06453d32f8718f573368009b12d19b55c58c2712e439408d26e1a9217f56");
        assert_eq!(sha256_file(&f.join("SH.evidence.json.gz")).unwrap(),"de21b3fa16fe8df6ca43421260acfd912391ef513e6bf5ead86853be4b9e81f6");
        assert_eq!(sha256_file(&f.join("LSU.evidence.json.gz")).unwrap(),"1a259082f41ef44cdfca2d02b64dc7ef17c34473b693db56a292f5d35ecc9c4b");
    }

    #[test]
    fn python_exclusive_quantiles_match_reference_interpolation() {
        assert_eq!(quantiles_exclusive(&[1.0, 3.0, 5.0, 7.0], 4), vec![1.5, 4.0, 6.5]);
        let q = quantiles_exclusive(&[1.0, 2.0], 100);
        assert_eq!(q.len(), 99);
        assert!((q[0] - 0.03).abs() < 1e-12);
        assert!((q[32] - 0.99).abs() < 1e-12);
        assert!((q[33] - 1.02).abs() < 1e-12);
        assert!((q[98] - 2.97).abs() < 1e-12);
    }

    #[test]
    fn python_general_float_format_matches_frozen_novelty_strings() {
        assert_eq!(format_g(248.75, 6), "248.75");
        assert_eq!(format_g(100.0, 8), "100");
        assert_eq!(format_g(0.00001, 8), "1e-05");
        assert_eq!(format_g0(248.75), "248.75");
    }

    #[test]
    fn equal_length_hw_distances_keep_argument_direction() {
        let forward = percent_distance("AAC", "AGA", -1).unwrap().unwrap();
        let reverse = percent_distance("AGA", "AAC", -1).unwrap().unwrap();
        assert!((forward - 200.0 / 3.0).abs() < 1e-12, "forward={forward}");
        assert!((reverse - 100.0 / 3.0).abs() < 1e-12, "reverse={reverse}");
    }

    #[test]
    fn frozen_m2b_assignment_rows_match_every_golden_cell() {
        let report=run_gate0(fixture()).unwrap_or_else(|e|panic!("loading/running the frozen Gate 0 data failed: {e:#}"));
        assert!(report.passed(),"{}",report.summary());
    }

    #[test]
    fn frozen_m2b_species_clade_rows_match_every_golden_cell() {
        let report=run_gate0b_species_clades(fixture()).unwrap_or_else(|e|panic!("loading/running the Gate 0b species-clade data failed: {e:#}"));
        assert!(report.passed(),"{}",report.summary());
    }
}
