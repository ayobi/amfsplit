//! AMF triage of Stage 1 reads: is each read Glomeromycota, from the best hit of its LR0R–FLR2
//! window against an LSU reference that also holds non-AMF lineages?
//!
//! Two steps, ported from the Stage 1 prototype (amfsplit v0.3.1) with byte-identical outputs:
//!   validate   `validate_windows.py`: vsearch global search of a window FASTA against a
//!              reference, the hits table and a summary of best hits;
//!   triage     `triage.py`: a call per read (AMF, unresolved, non_AMF, no_hit, no_window), the
//!              AMF read IDs, and AMF-only copies of every Stage 1 FASTA.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use anyhow::{Context, Result, bail};
use regex::Regex;

use crate::pycompat::{
    Json, float_repr, median_f64, most_common, py_first_token, py_round, read_tsv_rows,
};

/// The v0.3.1 prototype's default non-AMF genera: the four outgroups of an older Delavaux LSU
/// release. V18 (V18_LSUDB_052025.fasta) holds 35 non-AMF records, of which only Mortierella is
/// listed here, so the port has no default: the caller says which targets are AMF.
pub const V031_DEFAULT_OUTGROUPS: &str = "Mortierella,Rhodotorula,Exophiala,Citrus";

/// vsearch options of `validate_windows.py`, which the hits table depends on.
pub const VSEARCH_FIXED: [&str; 8] = [
    "--strand",
    "both",
    "--maxaccepts",
    "16",
    "--maxrejects",
    "64",
    "--userfields",
    "query+target+id+qcov+tcov+alnlen",
];

#[derive(Clone, Debug)]
pub struct ValidateParams {
    pub windows: PathBuf,
    pub reference: PathBuf,
    pub out: String,
    pub id: f64,
    pub threads: i64,
    pub vsearch: PathBuf,
}

/// Best hit of one query in a vsearch hits table.
#[derive(Clone, Debug)]
struct BestHit {
    target: String,
    id: f64,
    qcov: f64,
    tcov: f64,
}

fn count_fasta_headers(path: &Path) -> Result<i64> {
    let f = fs::File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let mut n = 0;
    for line in BufReader::new(f).split(b'\n') {
        if line?.first() == Some(&b'>') {
            n += 1;
        }
    }
    Ok(n)
}

/// The exact vsearch command line of `validate_windows.py`.
pub fn vsearch_command(p: &ValidateParams, hits_path: &Path) -> Command {
    let mut cmd = Command::new(&p.vsearch);
    cmd.arg("--usearch_global")
        .arg(&p.windows)
        .arg("--db")
        .arg(&p.reference)
        .arg("--id")
        .arg(float_repr(p.id))
        .args(&VSEARCH_FIXED[..6])
        .arg("--userout")
        .arg(hits_path)
        .args(&VSEARCH_FIXED[6..])
        .arg("--threads")
        .arg(p.threads.to_string())
        .arg("--quiet");
    cmd
}

/// `validate_windows.py`: run vsearch, write `{out}.hits.tsv` and `{out}.summary.json`; returns
/// the text the prototype prints on stdout.
pub fn validate(p: &ValidateParams) -> Result<String> {
    let hits_path = PathBuf::from(format!("{}.hits.tsv", p.out));
    let status = vsearch_command(p, &hits_path)
        .status()
        .with_context(|| format!("cannot run vsearch ({})", p.vsearch.display()))?;
    if !status.success() {
        bail!("vsearch failed ({status})");
    }
    let n_query = count_fasta_headers(&p.windows)?;
    let mut order: Vec<String> = Vec::new();
    let mut best: HashMap<String, BestHit> = HashMap::new();
    let text = fs::read_to_string(&hits_path)
        .with_context(|| format!("cannot read {}", hits_path.display()))?;
    for line in text.split_inclusive('\n') {
        let f: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        if f.len() != 6 {
            bail!("{}: expected 6 fields: {line:?}", hits_path.display());
        }
        let num = |s: &str| {
            s.parse::<f64>()
                .with_context(|| format!("not a number: {s:?}"))
        };
        let hit = BestHit {
            target: f[1].to_string(),
            id: num(f[2])?,
            qcov: num(f[3])?,
            tcov: num(f[4])?,
        };
        match best.get_mut(f[0]) {
            None => {
                order.push(f[0].to_string());
                best.insert(f[0].to_string(), hit);
            }
            Some(b) if hit.id > b.id => *b = hit,
            Some(_) => {}
        }
    }
    let vals =
        |get: fn(&BestHit) -> f64| -> Vec<f64> { order.iter().map(|q| get(&best[q])).collect() };
    let med = |get: fn(&BestHit) -> f64| -> Json {
        let v = vals(get);
        if v.is_empty() {
            Json::Null
        } else {
            Json::Float(py_round(median_f64(&v), 1))
        }
    };
    let denom = n_query.max(1) as f64;
    let ge97 = order.iter().filter(|q| best[*q].id >= 97.0).count() as i64;
    let top = most_common(order.iter().map(|q| best[q].target.as_str()), Some(5));
    let summary = Json::obj(vec![
        ("windows", Json::Int(n_query)),
        ("with_hit", Json::Int(order.len() as i64)),
        (
            "with_hit_pct",
            Json::Float(py_round((100 * order.len() as i64) as f64 / denom, 1)),
        ),
        ("best_identity_median", med(|b| b.id)),
        ("query_coverage_median", med(|b| b.qcov)),
        ("target_coverage_median", med(|b| b.tcov)),
        (
            "identity_ge_97_pct",
            Json::Float(py_round((100 * ge97) as f64 / denom, 1)),
        ),
        (
            "top_targets",
            Json::Arr(
                top.into_iter()
                    .map(|(t, n)| Json::Arr(vec![Json::Str(t), Json::Int(n)]))
                    .collect(),
            ),
        ),
    ]);
    let text = summary.dump_indent1();
    let sum_path = format!("{}.summary.json", p.out);
    fs::write(&sum_path, &text).with_context(|| format!("cannot write {sum_path}"))?;
    Ok(text + "\n")
}

#[derive(Clone, Debug)]
pub struct TriageParams {
    pub hits: PathBuf,
    pub out: String,
    pub records: Option<PathBuf>,
    pub fasta_prefix: Option<String>,
    /// Non-AMF genera, comma-separated (the prototype's way).
    pub outgroups: Option<String>,
    /// Regex marking AMF targets; overrides `outgroups`, as in the prototype.
    pub amf_pattern: Option<String>,
    /// FASTA whose record IDs are the AMF targets (e.g. V18_LSUDB_052025_AMFONLY.fasta).
    pub amf_reference: Option<PathBuf>,
    pub min_id: f64,
    pub keep_unresolved: bool,
}

impl TriageParams {
    pub fn new(hits: PathBuf, out: String) -> Self {
        TriageParams {
            hits,
            out,
            records: None,
            fasta_prefix: None,
            outgroups: None,
            amf_pattern: None,
            amf_reference: None,
            min_id: 90.0,
            keep_unresolved: false,
        }
    }
}

#[derive(Clone, Debug)]
struct Call {
    call: &'static str,
    target: String,
    identity: Option<f64>,
}

/// What a triage run produced: the stdout text of the prototype, and notes for the console.
#[derive(Clone, Debug)]
pub struct TriageReport {
    pub stdout: String,
    pub notes: Vec<String>,
}

/// How a hit target is judged AMF.
enum AmfRule {
    /// Target ID listed in an AMF-only reference.
    Reference {
        ids: HashSet<String>,
        path: String,
        sha256: String,
    },
    /// Regex search on the target name (prototype: `--amf-pattern`).
    Pattern(Regex),
    /// Second underscore-separated field not among these genera (prototype: `--outgroups`).
    Outgroups(Vec<String>),
}

/// Record IDs of a FASTA (plain or gzipped): the first whitespace-separated token of each header.
fn reference_ids(path: &Path) -> Result<HashSet<String>> {
    let mut ids = HashSet::new();
    let mut r = crate::seqio::open_maybe_gz(path)?;
    let mut line = Vec::new();
    loop {
        line.clear();
        if r.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.first() == Some(&b'>') {
            let id = py_first_token(&line[1..])
                .with_context(|| format!("{}: empty header", path.display()))?;
            ids.insert(String::from_utf8_lossy(id).into_owned());
        }
    }
    Ok(ids)
}

fn sha256_file(path: &Path) -> Result<String> {
    use sha2::{Digest, Sha256};
    let mut h = Sha256::new();
    let mut f = fs::File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    std::io::copy(&mut f, &mut h)?;
    Ok(format!("{:x}", h.finalize()))
}

fn amf_rule(p: &TriageParams) -> Result<AmfRule> {
    match (&p.amf_reference, &p.amf_pattern, &p.outgroups) {
        (Some(_), Some(_), _) | (Some(_), _, Some(_)) => {
            bail!("--amf-reference cannot be combined with --amf-pattern or --outgroups")
        }
        (Some(path), None, None) => {
            let ids = reference_ids(path)?;
            if ids.is_empty() {
                bail!("--amf-reference {} holds no FASTA records", path.display());
            }
            Ok(AmfRule::Reference {
                ids,
                path: path.display().to_string(),
                sha256: sha256_file(path)?,
            })
        }
        (None, Some(pat), _) => {
            Ok(AmfRule::Pattern(Regex::new(pat).with_context(|| {
                format!("--amf-pattern is not a valid regex: {pat:?}")
            })?))
        }
        (None, None, Some(list)) => {
            let mut set: Vec<String> = Vec::new();
            for g in list.split(',') {
                let g = g.trim();
                if !g.is_empty() && !set.iter().any(|x| x == g) {
                    set.push(g.to_string());
                }
            }
            set.sort();
            Ok(AmfRule::Outgroups(set))
        }
        (None, None, None) => bail!(
            "say which hit targets are AMF:\n  \
             --amf-reference FASTA   the AMF records of the searched reference; with the Delavaux V18 LSU \
             database, search V18_LSUDB_052025.fasta and give V18_LSUDB_052025_AMFONLY.fasta here (recommended)\n  \
             --amf-pattern REGEX     a regex on the target name, e.g. 'Glomeromycota' for a reference with lineages in its headers\n  \
             --outgroups GENERA      non-AMF genera; amfsplit v0.3.1 defaulted to {V031_DEFAULT_OUTGROUPS}, \
             the outgroups of an older Delavaux release (only Mortierella is in V18)"
        ),
    }
}

/// Names of the files `glob.glob(f"{prefix}.*.fasta")` returns, sorted as Python sorts them.
fn glob_prefix_fastas(prefix: &str) -> Result<Vec<PathBuf>> {
    let p = Path::new(prefix);
    let dir = p
        .parent()
        .filter(|d| !d.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let base = p
        .file_name()
        .map(|b| b.to_string_lossy().into_owned())
        .unwrap_or_default();
    let head = format!("{base}.");
    let mut out: Vec<String> = Vec::new();
    for entry in fs::read_dir(dir).with_context(|| format!("cannot list {}", dir.display()))? {
        let name = entry?.file_name().to_string_lossy().into_owned();
        if name.len() >= head.len() + ".fasta".len()
            && name.starts_with(&head)
            && name.ends_with(".fasta")
        {
            out.push(if p.parent().is_some_and(|d| !d.as_os_str().is_empty()) {
                format!("{}/{}", dir.display(), name)
            } else {
                name
            });
        }
    }
    out.sort();
    Ok(out.into_iter().map(PathBuf::from).collect())
}

fn fasta_ids(path: &Path) -> Result<Vec<String>> {
    let f = fs::File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let mut ids = Vec::new();
    for line in BufReader::new(f).split(b'\n') {
        let line = line?;
        if line.first() == Some(&b'>') {
            let id = py_first_token(&line[1..])
                .with_context(|| format!("{}: empty header", path.display()))?;
            ids.push(String::from_utf8_lossy(id).into_owned());
        }
    }
    Ok(ids)
}

/// Read IDs named AMF on a best hit that V18 marks as non-AMF: its header ends in a fungal
/// lineage tag (Asco, Basid, Chyt, Zygo) or names Homo, Arabidopsis, Oryza or Geosiphon. Only a
/// console hint; the outputs stay those of the prototype.
fn nonamf_hint(calls: &[(String, Call)]) -> Option<String> {
    let tagged = |t: &str| {
        ["_Asco", "_Basid", "_Chyt", "_Zygo"]
            .iter()
            .any(|s| t.ends_with(s))
            || ["_Homo_", "_Arabidopsis_", "_Oryza_", "_Geosiphon_"]
                .iter()
                .any(|s| t.contains(s))
    };
    let flagged: Vec<&str> = calls
        .iter()
        .filter(|(_, c)| matches!(c.call, "AMF" | "unresolved") && tagged(&c.target))
        .map(|(_, c)| c.target.as_str())
        .collect();
    if flagged.is_empty() {
        return None;
    }
    let examples = most_common(flagged.iter().copied(), Some(3))
        .into_iter()
        .map(|(t, _)| t)
        .collect::<Vec<_>>()
        .join(", ");
    Some(format!(
        "note: {} read(s) were called AMF or unresolved on a best hit that the reference labels non-AMF \
         (e.g. {examples}). --outgroups names only some of the reference's non-AMF genera. With the \
         Delavaux V18 LSU database, use --amf-reference V18_LSUDB_052025_AMFONLY.fasta instead.",
        flagged.len()
    ))
}

/// `triage.py`: calls per read, `{out}.triage.tsv`, `{out}.amf.ids`, `{out}.triage.summary.json`,
/// and with a FASTA prefix, `<file>.amf.fasta` beside every `{prefix}.*.fasta`.
pub fn triage(p: &TriageParams) -> Result<TriageReport> {
    let rule = amf_rule(p)?;
    let is_amf = |target: &str| -> bool {
        match &rule {
            AmfRule::Reference { ids, .. } => ids.contains(target),
            AmfRule::Pattern(re) => re.is_match(target),
            AmfRule::Outgroups(outgroups) => {
                let parts: Vec<&str> = target.split('_').collect();
                parts.len() < 2 || !outgroups.iter().any(|g| g == parts[1])
            }
        }
    };

    // best hit per query, first seen kept on ties
    let mut best_order: Vec<String> = Vec::new();
    let mut best: HashMap<String, (String, f64)> = HashMap::new();
    let text =
        fs::read_to_string(&p.hits).with_context(|| format!("cannot read {}", p.hits.display()))?;
    for line in text.split_inclusive('\n') {
        let f: Vec<&str> = line.trim_end_matches('\n').split('\t').collect();
        if f.len() < 3 {
            bail!("{}: expected at least 3 fields: {line:?}", p.hits.display());
        }
        let ident: f64 = f[2]
            .parse()
            .with_context(|| format!("identity is not a number: {:?}", f[2]))?;
        match best.get_mut(f[0]) {
            None => {
                best_order.push(f[0].to_string());
                best.insert(f[0].to_string(), (f[1].to_string(), ident));
            }
            Some(b) if ident > b.1 => *b = (f[1].to_string(), ident),
            Some(_) => {}
        }
    }

    let mut notes = Vec::new();
    if let AmfRule::Reference { ids, path, .. } = &rule {
        let targets: HashSet<&str> = text
            .split_inclusive('\n')
            .filter_map(|l| l.split('\t').nth(1))
            .collect();
        let n_amf = targets.iter().filter(|t| ids.contains(**t)).count();
        if !targets.is_empty() && n_amf == 0 {
            bail!(
                "none of the {} hit targets is in --amf-reference {path}: the search reference and the \
                 AMF reference must come from the same database release",
                targets.len()
            );
        }
        if !targets.is_empty() && n_amf == targets.len() {
            notes.push(format!(
                "note: all {} hit targets are in the AMF reference, so no read can be called non_AMF. \
                 If `amfsplit validate` searched the AMF-only file, search the full reference instead \
                 (with Delavaux V18: V18_LSUDB_052025.fasta).",
                targets.len()
            ));
        }
    }

    // calls, in the prototype's insertion order: hits first, then the rest of the universe
    let mut calls: Vec<(String, Call)> = Vec::new();
    let mut index: HashMap<String, usize> = HashMap::new();
    for q in &best_order {
        let (t, ident) = &best[q];
        let call = if !is_amf(t) {
            "non_AMF"
        } else if *ident >= p.min_id {
            "AMF"
        } else {
            "unresolved"
        };
        index.insert(q.clone(), calls.len());
        calls.push((
            q.clone(),
            Call {
                call,
                target: t.clone(),
                identity: Some(*ident),
            },
        ));
    }
    let universe: Vec<String> = match &p.records {
        None => best_order.clone(),
        Some(r) => {
            let text =
                fs::read_to_string(r).with_context(|| format!("cannot read {}", r.display()))?;
            let rows = read_tsv_rows(&text);
            let Some(header) = rows.first() else {
                bail!("{} is empty", r.display())
            };
            let col = header
                .iter()
                .position(|h| h == "read_id")
                .with_context(|| format!("{} has no read_id column", r.display()))?;
            rows[1..]
                .iter()
                .map(|row| row.get(col).cloned().unwrap_or_default())
                .collect()
        }
    };
    let mut windowed: HashSet<String> = best.keys().cloned().collect();
    if let Some(prefix) = &p.fasta_prefix {
        let lsu = PathBuf::from(format!("{prefix}.LSU_LROR-FLR2.fasta"));
        if lsu.exists() {
            windowed.extend(fasta_ids(&lsu)?);
        }
    }
    for rid in &universe {
        if !index.contains_key(rid) {
            let call = if windowed.contains(rid) {
                "no_hit"
            } else {
                "no_window"
            };
            index.insert(rid.clone(), calls.len());
            calls.push((
                rid.clone(),
                Call {
                    call,
                    target: String::new(),
                    identity: None,
                },
            ));
        }
    }

    let keep: &[&str] = if p.keep_unresolved {
        &["AMF", "unresolved"]
    } else {
        &["AMF"]
    };
    let mut tsv = String::from("read_id\tcall\tbest_target\tidentity\n");
    let row_ids: Vec<String> = if p.records.is_some() {
        universe.clone()
    } else {
        let mut ids: Vec<String> = calls.iter().map(|(r, _)| r.clone()).collect();
        ids.sort();
        ids
    };
    for rid in &row_ids {
        let c = &calls[index[rid]].1;
        let ident = c.identity.map_or(String::new(), float_repr);
        tsv.push_str(&format!("{rid}\t{}\t{}\t{ident}\n", c.call, c.target));
    }
    let tsv_path = format!("{}.triage.tsv", p.out);
    fs::write(&tsv_path, tsv).with_context(|| format!("cannot write {tsv_path}"))?;

    let mut amf_ids: Vec<&str> = calls
        .iter()
        .filter(|(_, c)| keep.contains(&c.call))
        .map(|(r, _)| r.as_str())
        .collect();
    amf_ids.sort();
    amf_ids.dedup();
    let ids_path = format!("{}.amf.ids", p.out);
    fs::write(
        &ids_path,
        amf_ids.iter().map(|r| format!("{r}\n")).collect::<String>(),
    )
    .with_context(|| format!("cannot write {ids_path}"))?;
    let amf_set: HashSet<&str> = amf_ids.iter().copied().collect();

    let n = calls.len() as i64;
    let counts = most_common(calls.iter().map(|(_, c)| c.call), None);
    let non_amf = most_common(
        calls
            .iter()
            .filter(|(_, c)| c.call == "non_AMF")
            .map(|(_, c)| {
                // t.split("_")[1] if "_" in t else t
                c.target.split('_').nth(1).unwrap_or(c.target.as_str())
            }),
        Some(6),
    );
    let pairs = |v: &[(String, i64)]| {
        Json::Arr(
            v.iter()
                .map(|(t, k)| Json::Arr(vec![Json::Str(t.clone()), Json::Int(*k)]))
                .collect(),
        )
    };
    let calls_json = Json::Obj(
        counts
            .iter()
            .map(|(k, v)| {
                let pct = py_round((100 * v) as f64 / n.max(1) as f64, 1);
                (
                    k.clone(),
                    Json::obj(vec![("n", Json::Int(*v)), ("pct", Json::Float(pct))]),
                )
            })
            .collect(),
    );
    let mut summary = vec![
        ("reads".to_string(), Json::Int(n)),
        ("calls".to_string(), calls_json),
        ("non_AMF_best_hits".to_string(), pairs(&non_amf)),
        ("params".to_string(), {
            let mut params = vec![
                ("min_id".to_string(), Json::Float(p.min_id)),
                (
                    "amf_pattern".to_string(),
                    p.amf_pattern.clone().map_or(Json::Null, Json::Str),
                ),
                (
                    "outgroups".to_string(),
                    match &rule {
                        AmfRule::Outgroups(o) => {
                            Json::Arr(o.iter().map(|g| Json::Str(g.clone())).collect())
                        }
                        _ => Json::Null,
                    },
                ),
            ];
            // only with --amf-reference, so the prototype's modes keep its exact summary
            if let AmfRule::Reference { ids, path, sha256 } = &rule {
                params.push((
                    "amf_reference".to_string(),
                    Json::obj(vec![
                        ("path", Json::Str(path.clone())),
                        ("sha256", Json::Str(sha256.clone())),
                        ("records", Json::Int(ids.len() as i64)),
                    ]),
                ));
            }
            params.push(("keep_unresolved".to_string(), Json::Bool(p.keep_unresolved)));
            Json::Obj(params)
        }),
    ];

    if let Some(prefix) = &p.fasta_prefix {
        let mut written: Vec<(String, Json)> = Vec::new();
        for path in glob_prefix_fastas(prefix)? {
            let ps = path.to_string_lossy().into_owned();
            if ps.ends_with(".amf.fasta") {
                continue;
            }
            let out_path = format!("{}.amf.fasta", &ps[..ps.len() - ".fasta".len()]);
            let mut src =
                BufReader::new(fs::File::open(&path).with_context(|| format!("cannot open {ps}"))?);
            let mut dst = std::io::BufWriter::new(
                fs::File::create(&out_path).with_context(|| format!("cannot write {out_path}"))?,
            );
            let (mut n_out, mut keep_rec) = (0i64, false);
            let mut line = Vec::new();
            loop {
                line.clear();
                if src.read_until(b'\n', &mut line)? == 0 {
                    break;
                }
                // Python's text mode reads CRLF as LF and writes LF
                if line.ends_with(b"\r\n") {
                    line.truncate(line.len() - 2);
                    line.push(b'\n');
                }
                if line.first() == Some(&b'>') {
                    let id = py_first_token(&line[1..])
                        .map(|t| String::from_utf8_lossy(t).into_owned())
                        .unwrap_or_default();
                    keep_rec = amf_set.contains(id.as_str());
                    n_out += i64::from(keep_rec);
                }
                if keep_rec {
                    dst.write_all(&line)?;
                }
            }
            dst.flush()?;
            let name = out_path.rsplit('/').next().unwrap_or(&out_path).to_string();
            written.push((name, Json::Int(n_out)));
        }
        summary.push(("amf_fastas".to_string(), Json::Obj(written)));
    }
    let sum_path = format!("{}.triage.summary.json", p.out);
    fs::write(&sum_path, Json::Obj(summary.clone()).dump_indent1())
        .with_context(|| format!("cannot write {sum_path}"))?;
    let stdout = Json::Obj(summary.into_iter().take(3).collect()).dump_indent1() + "\n";
    Ok(TriageReport {
        stdout,
        notes: {
            if matches!(rule, AmfRule::Outgroups(_)) {
                notes.extend(nonamf_hint(&calls));
            }
            notes
        },
    })
}
