//! Withholding rounds (Protocol v3, section 5): V18 with a taxon's sequences removed, re-placed exactly as
//! the frozen LSU-removal rounds did (benchmarks/m4_external_holdout/external_holdout.py, `run_fold`):
//! 01_reference prunes the withheld tips from backbone.topology.newick (Bio.Phylo `prune`) and re-aligns
//! the kept reference with MAFFT L-INS-i; 02_model refits branch lengths and GTR+G4 with raxml-ng
//! --evaluate (seed 20260929, at most 4 threads) on the fixed topology; 03_queries adds the queries with
//! MAFFT --addfragments --keeplength; 04_epa places them with EPA-ng; 05_summary is the frozen summary.
use super::newick::Tree;
use super::placement::{compare_up_to_ties, placement_tsv, read_fasta_ordered, summarize_with_ties, TableMatch};
use anyhow::{anyhow, bail, Context, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

pub struct WithholdOptions {
    /// frozen V18 reference folder: reference.unaligned.fasta, reference_metadata.tsv, backbone.topology.newick
    pub frozen_reference: PathBuf,
    pub removed: Vec<String>,
    pub queries: PathBuf,
    pub out: PathBuf,
    pub threads: usize,
    pub mafft: String,
    pub raxml_ng: String,
    pub epa_ng: String,
    /// a frozen round folder (runs/<key>) to compare every stage with
    pub compare_with: Option<PathBuf>,
    /// stage 06, naming the fold's queries: (naming fixture, holdout_plan.json, fold key)
    pub naming: Option<(PathBuf, PathBuf, String)>,
    /// a finished round of the same fold whose 01_reference and 02_model are reused, so only the new
    /// queries are aligned and placed (withheld K: a species' K sequences in its own round)
    pub reuse_reference: Option<PathBuf>,
}

pub struct WithholdReport { pub stages: Vec<(String, String)>, pub placement: Option<TableMatch>, pub table: PathBuf, pub named: Vec<crate::TableReport> }

fn tool(name: &str) -> Result<PathBuf> {
    let p = Path::new(name);
    if p.components().count() > 1 { return Ok(p.to_path_buf()); }
    std::env::var_os("PATH").and_then(|ps| std::env::split_paths(&ps).map(|d| d.join(name)).find(|c| c.is_file()))
        .ok_or_else(|| anyhow!("{name} is not on PATH"))
}

fn run(cmd: &[String], log: &Path, stdout: Option<&Path>) -> Result<()> {
    let err = fs::File::create(log)?;
    let mut c = std::process::Command::new(&cmd[0]);
    c.args(&cmd[1..]);
    match stdout { Some(p) => { c.stdout(fs::File::create(p)?).stderr(err); } None => { c.stdout(err.try_clone()?).stderr(err); } }
    let st = c.status().with_context(|| format!("starting {}", cmd[0]))?;
    if !st.success() { bail!("{} failed ({st}); see {}", cmd[0], log.display()); }
    Ok(())
}

fn write_fasta(path: &Path, seqs: &[(String, String)]) -> Result<()> {
    let mut s = String::new();
    for (k, v) in seqs { s.push('>'); s.push_str(k); s.push('\n'); s.push_str(v); s.push('\n'); }
    Ok(fs::write(path, s)?)
}

/// external_holdout.py `splits`: each node's tip set against the rest, both sides above one tip,
/// written as the smaller of the two sorted tip lists.
fn splits(t: &Tree) -> BTreeSet<Vec<String>> {
    let all: BTreeSet<String> = t.tips().into_iter().collect();
    let mut out = BTreeSet::new();
    for side in t.clade_tip_sets() {
        let other: BTreeSet<String> = all.difference(&side).cloned().collect();
        if side.len() > 1 && other.len() > 1 {
            let a: Vec<String> = side.into_iter().collect();
            let b: Vec<String> = other.into_iter().collect();
            out.insert(if a < b { a } else { b });
        }
    }
    out
}

pub fn run_withhold(o: &WithholdOptions) -> Result<WithholdReport> {
    let s = |p: &Path| p.display().to_string();
    let fr = &o.frozen_reference;
    let removed: BTreeSet<String> = o.removed.iter().cloned().collect();
    let [r1, r2, r3, r4, r5] = ["01_reference", "02_model", "03_queries", "04_epa", "05_summary"].map(|d| o.out.join(d));
    for d in [&r1, &r2, &r3, &r4, &r5] { fs::create_dir_all(d)?; }
    let mafft = tool(&o.mafft)?;
    if let Some(src) = &o.reuse_reference {
        for d in ["01_reference", "02_model"] {
            for e in fs::read_dir(src.join(d)).with_context(|| format!("reading {}", src.join(d).display()))? {
                let p = e?.path();
                if p.is_file() { fs::copy(&p, o.out.join(d).join(p.file_name().unwrap()))?; }
            }
        }
        // the reused round must be this fold: same withheld tips, same kept reference
        let kept_ids: BTreeSet<String> = read_fasta_ordered(&o.out.join("01_reference/aligned.fasta"))?.into_iter().map(|(x, _)| x).collect();
        let raw_ids: BTreeSet<String> = read_fasta_ordered(&fr.join("reference.unaligned.fasta"))?.into_iter().map(|(x, _)| x).collect();
        let expect: BTreeSet<String> = raw_ids.difference(&removed).cloned().collect();
        if kept_ids != expect { bail!("{} is not a round of this fold (its reference holds other taxa)", src.display()); }
    }
    // 01_reference
    let raw = read_fasta_ordered(&fr.join("reference.unaligned.fasta"))?;
    let meta_text = fs::read_to_string(fr.join("reference_metadata.tsv"))?;
    let mut lines = meta_text.lines();
    let header = lines.next().ok_or_else(|| anyhow!("empty reference_metadata.tsv"))?;
    let cols: Vec<&str> = header.split('\t').collect();
    let (id_col, amf_col) = (cols.iter().position(|c| *c == "ref_id").ok_or_else(|| anyhow!("no ref_id column"))?,
                             cols.iter().position(|c| *c == "is_amf").ok_or_else(|| anyhow!("no is_amf column"))?);
    let raw_ids: BTreeSet<&String> = raw.iter().map(|(x, _)| x).collect();
    if removed.iter().any(|x| !raw_ids.contains(x)) { bail!("removal includes taxa missing from the reference"); }
    for l in meta_text.lines().skip(1) {
        let f: Vec<&str> = l.split('\t').collect();
        if removed.contains(f[id_col]) && f[amf_col] != "1" { bail!("removal includes an outgroup ({})", f[id_col]); }
    }
    let kept: Vec<(String, String)> = raw.iter().filter(|(x, _)| !removed.contains(x)).cloned().collect();
    let aligned_path = r1.join("aligned.fasta");
    if !aligned_path.is_file() {
        write_fasta(&r1.join("raw.fasta"), &kept)?;
        let mut meta = String::from(header); meta.push('\n');
        for l in meta_text.lines().skip(1) { if !removed.contains(l.split('\t').nth(id_col).unwrap_or("")) { meta.push_str(l); meta.push('\n'); } }
        fs::write(r1.join("metadata.tsv"), meta)?;
        let mut tree = Tree::parse(&fs::read_to_string(fr.join("backbone.topology.newick"))?)?;
        for tip in &removed { tree.prune(tip)?; }
        let kept_ids: BTreeSet<String> = kept.iter().map(|(x, _)| x.clone()).collect();
        if tree.tips().into_iter().collect::<BTreeSet<_>>() != kept_ids { bail!("pruning produced the wrong taxon set"); }
        fs::write(r1.join("topology.newick"), tree.to_newick())?;
        let cmd = vec![s(&mafft), "--localpair".into(), "--maxiterate".into(), "1000".into(), "--thread".into(), o.threads.to_string(),
                       "--threadit".into(), "0".into(), s(&r1.join("raw.fasta"))];
        run(&cmd, &r1.join("mafft.log"), Some(&aligned_path))?;
        let aligned = read_fasta_ordered(&aligned_path)?;
        let am: BTreeMap<&String, &String> = aligned.iter().map(|(x, q)| (x, q)).collect();
        if am.len() != kept.len() || aligned.iter().map(|(_, q)| q.len()).collect::<BTreeSet<_>>().len() != 1 { bail!("invalid reference alignment"); }
        for (x, q) in &kept {
            let a = am.get(x).ok_or_else(|| anyhow!("{x} missing from the alignment"))?;
            if a.replace(['-', '.'], "") != q.replace(['-', '.'], "") { bail!("reference residues changed for {x}"); }
        }
        fs::write(r1.join("commands.json"), serde_json::to_string(&[cmd])?)?;
    }
    // 02_model
    let best_tree = r2.join("fit.raxml.bestTree");
    if !best_tree.is_file() {
        let cmd = vec![s(&tool(&o.raxml_ng)?), "--evaluate".into(), "--msa".into(), s(&aligned_path), "--tree".into(), s(&r1.join("topology.newick")),
                       "--model".into(), "GTR+G4".into(), "--threads".into(), o.threads.min(4).to_string(), "--seed".into(), "20260929".into(),
                       "--prefix".into(), s(&r2.join("fit")), "--force".into(), "perf_threads".into()];
        run(&cmd, &r2.join("raxml.log"), None)?;
        let fixed = Tree::parse(&fs::read_to_string(r1.join("topology.newick"))?)?;
        let fit = Tree::parse(&fs::read_to_string(&best_tree)?)?;
        if splits(&fixed) != splits(&fit) { bail!("raxml-ng changed the fixed topology"); }
        fs::write(r2.join("commands.json"), serde_json::to_string(&[cmd])?)?;
    }
    // 03_queries
    let qa = r3.join("queries.aligned.fasta");
    let queries = r3.join("queries.fasta");
    if !qa.is_file() {
        fs::copy(&o.queries, &queries)?;
        let cmd = vec![s(&mafft), "--auto".into(), "--addfragments".into(), s(&queries), "--keeplength".into(), "--mapout".into(),
                       "--thread".into(), o.threads.min(4).to_string(), s(&aligned_path)];
        run(&cmd, &r3.join("mafft.log"), Some(&r3.join("combined.fasta")))?;
        let qseq = read_fasta_ordered(&queries)?;
        let base = read_fasta_ordered(&aligned_path)?;
        let full: BTreeMap<String, String> = read_fasta_ordered(&r3.join("combined.fasta"))?.into_iter().collect();
        if base.iter().any(|(x, q)| full.get(x) != Some(q)) || full.len() != base.len() + qseq.len() { bail!("query addition altered the reference or the IDs"); }
        let out: Vec<(String, String)> = qseq.iter().map(|(x, _)| (x.clone(), full[x].clone())).collect();
        write_fasta(&qa, &out)?;
        fs::write(r3.join("commands.json"), serde_json::to_string(&[cmd])?)?;
    }
    // 04_epa
    let jplace = r4.join("epa_result.jplace");
    if !jplace.is_file() {
        let cmd = vec![s(&tool(&o.epa_ng)?), "--ref-msa".into(), s(&aligned_path), "--tree".into(), s(&best_tree), "--query".into(), s(&qa),
                       "--model".into(), s(&r2.join("fit.raxml.bestModel")), "--outdir".into(), s(&r4), "--threads".into(), o.threads.to_string(),
                       "--no-heur".into(), "--no-pre-mask".into(), "--filter-min-lwr".into(), "0".into(), "--filter-max".into(), "2000".into(), "--precision".into(), "12".into()];
        run(&cmd, &r4.join("epa.log"), None)?;
        fs::write(r4.join("commands.json"), serde_json::to_string(&[cmd])?)?;
    }
    // 05_summary
    let (rows, ties) = summarize_with_ties(&jplace, &r1.join("metadata.tsv"), &queries, &qa, 0.95)?;
    let table = r5.join("lsu_placements.copies.tsv");
    fs::write(&table, placement_tsv(&rows))?;
    let mut stages = Vec::new();
    let mut placement = None;
    if let Some(f) = &o.compare_with {
        for rel in ["01_reference/raw.fasta", "01_reference/metadata.tsv", "01_reference/topology.newick", "01_reference/aligned.fasta",
                    "02_model/fit.raxml.bestTree", "02_model/fit.raxml.bestModel", "03_queries/queries.aligned.fasta"] {
            let (a, b) = (o.out.join(rel), f.join(rel));
            let verdict = if !b.is_file() { "frozen file absent".to_owned() }
                else if fs::read(&a)? == fs::read(&b)? { "identical".to_owned() } else { "DIFFERS".to_owned() };
            stages.push((rel.to_owned(), verdict));
        }
        let frozen_table = f.join("05_summary/lsu_placements.copies.tsv");
        if frozen_table.is_file() { placement = Some(compare_up_to_ties(&rows, &ties, &frozen_table)?); }
    }
    // 06_evaluation_v2: the fold's queries named with its held cultures from this round's placements
    let named = match &o.naming {
        Some((fixture, plan, key)) => crate::name_withheld_round(fixture, plan, key, &table, &o.out.join("06_evaluation_v2"),
                                                                  o.compare_with.as_ref().map(|f| f.join("06_evaluation_v2")).as_deref())?,
        None => Vec::new(),
    };
    Ok(WithholdReport { stages, placement, table, named })
}

/// The V18 reference IDs a fold of the frozen plan withholds (holdout_plan.json).
pub fn fold_removed(plan: &Path, key: &str) -> Result<Vec<String>> {
    let v: serde_json::Value = serde_json::from_str(&fs::read_to_string(plan)?)?;
    let fold = v["folds"].as_array().and_then(|f| f.iter().find(|x| x["key"] == key))
        .ok_or_else(|| anyhow!("no fold {key} in {}", plan.display()))?;
    Ok(fold["removed_ref_ids"].as_array().map(|ids| ids.iter().filter_map(|x| x.as_str().map(str::to_owned)).collect()).unwrap_or_default())
}

/// A placement table's verdict against the frozen one.
pub fn verdict(m: &TableMatch) -> String {
    match m {
        TableMatch::Exact => "identical".into(),
        TableMatch::UpToTiedBestEdges(rows) => format!("identical up to tied best edges ({} rows)", rows.len()),
        TableMatch::Different(n, first) => format!("DIFFERS in {n} rows: {}", first.iter().take(3).cloned().collect::<Vec<_>>().join("; ")),
    }
}
