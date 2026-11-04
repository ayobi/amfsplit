//! LSU placement summary: the frozen `lsu_placement.py` `summarize` (M2b, release 0.3.1) in Rust.
//!
//! Reads an EPA-ng jplace, the reference metadata, and the raw and aligned query windows, and writes
//! one row per query (`lsu_placements.copies.tsv`): the credible edge set at 0.95 mass, the
//! species-clade (SC) and genus-clade (GC) units those edges carry, the best edge and the alignment
//! retention that decides bridge eligibility. This is the table every naming gate reads.
//!
//! Edges get their clade units from the tree in the jplace: each edge's descendant tips are oriented
//! away from the first outgroup tip; an edge whose tips are all AMF and carry one species (genus)
//! label belongs to the largest edge of that label that contains it (ties to the smaller edge ID),
//! so polyphyletic labels keep separate units.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};
use serde_json::Value;

use super::*;

pub const PLACEMENT_HEADERS: [&str; 29] = [
    "accession", "aligned_bases", "raw_bases", "retained_fraction", "deleted_insertions", "ambiguous_aligned_bases",
    "bridge_alignment_eligible", "best_edge", "best_edge_lwr", "best_pendant_length", "best_distal_length", "lwr_sum",
    "lwr_entropy_nats", "credible_mass_target", "credible_mass_retained", "credible_edge_count", "credible_edges",
    "credible_species_clade_units", "amf_edge_mass", "external_species_clade", "external_species_clade_unit",
    "external_species_clade_mass", "external_genus_clade", "external_genus_clade_unit", "external_genus_clade_mass",
    "external_rank", "external_species_in_credible_edges", "external_genera_in_credible_edges", "interpretation",
];
const INTERPRETATION: &str = "conditional placement support; external labels are not validated culture predictions";

/// Python's math.fsum (CPython mathmodule.c, Shewchuk partials with half-even correction).
fn fsum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut p: Vec<f64> = Vec::new();
    for mut x in values {
        let mut i = 0;
        for j in 0..p.len() {
            let mut y = p[j];
            if x.abs() < y.abs() { std::mem::swap(&mut x, &mut y); }
            let hi = x + y;
            let yr = hi - x;
            let lo = y - yr;
            if lo != 0.0 { p[i] = lo; i += 1; }
            x = hi;
        }
        p.truncate(i);
        if x != 0.0 { p.push(x); }
    }
    let mut n = p.len();
    if n == 0 { return 0.0; }
    n -= 1;
    let mut hi = p[n];
    let mut lo = 0.0;
    while n > 0 {
        let x = hi;
        n -= 1;
        let y = p[n];
        hi = x + y;
        let yr = hi - x;
        lo = y - yr;
        if lo != 0.0 { break; }
    }
    if n > 0 && ((lo < 0.0 && p[n - 1] < 0.0) || (lo > 0.0 && p[n - 1] > 0.0)) {
        let y = lo * 2.0;
        let x = hi + y;
        let yr = x - hi;
        if y == yr { hi = x; }
    }
    hi
}

/// Python's math.isclose.
fn isclose(a: f64, b: f64, rel_tol: f64, abs_tol: f64) -> bool {
    if a == b { return true; }
    if a.is_infinite() || b.is_infinite() { return false; }
    let diff = (b - a).abs();
    diff <= (rel_tol * b).abs() || diff <= (rel_tol * a).abs() || diff <= abs_tol
}

/// A JSON number as Python's json module returns it: int without a fraction or exponent, else float.
#[derive(Clone, Copy, Debug)]
enum PyNum { Int(i64), Float(f64) }

impl PyNum {
    fn from(v: &Value) -> Result<Self> {
        let n = v.as_number().ok_or_else(|| anyhow!("jplace field is not a number: {v}"))?;
        if n.is_f64() { Ok(Self::Float(n.as_f64().unwrap())) }
        else { Ok(Self::Int(n.as_i64().ok_or_else(|| anyhow!("jplace integer out of range: {n}"))?)) }
    }
    fn f(self) -> f64 { match self { Self::Int(i) => i as f64, Self::Float(x) => x } }
    fn repr(self) -> String { match self { Self::Int(i) => i.to_string(), Self::Float(x) => py_float(x) } }
}

struct Node { name: Option<String>, edge: Option<usize>, children: Vec<usize> }

/// Edge-tagged Newick as EPA-ng writes it: `(R0001:0.1{0},...):0.2{5}...;`. Nodes come out in
/// postorder, children in file order, as Bio.Phylo's `find_clades(order='postorder')` visits them.
fn parse_jplace_tree(s: &str) -> Result<Vec<Node>> {
    fn subtree(b: &[u8], i: &mut usize, nodes: &mut Vec<Node>) -> Result<usize> {
        let mut children = Vec::new();
        if b.get(*i) == Some(&b'(') {
            *i += 1;
            loop {
                children.push(subtree(b, i, nodes)?);
                match b.get(*i) {
                    Some(b',') => *i += 1,
                    Some(b')') => { *i += 1; break; }
                    _ => bail!("malformed jplace tree near byte {}", *i),
                }
            }
        }
        let start = *i;
        while *i < b.len() && !b":,(){};".contains(&b[*i]) { *i += 1; }
        let name = (*i > start).then(|| String::from_utf8_lossy(&b[start..*i]).trim().to_owned()).filter(|n| !n.is_empty());
        if b.get(*i) == Some(&b':') {
            *i += 1;
            while *i < b.len() && !b",(){};".contains(&b[*i]) { *i += 1; }
        }
        let mut edge = None;
        if b.get(*i) == Some(&b'{') {
            let start = *i + 1;
            while *i < b.len() && b[*i] != b'}' { *i += 1; }
            edge = Some(String::from_utf8_lossy(&b[start..*i]).parse::<usize>().context("edge tag is not a number")?);
            *i += 1;
        }
        nodes.push(Node { name, edge, children });
        Ok(nodes.len() - 1)
    }
    let b = s.trim().as_bytes();
    let mut nodes = Vec::new();
    let mut i = 0;
    subtree(b, &mut i, &mut nodes)?;
    if b.get(i) != Some(&b';') { bail!("jplace tree does not end with ';'"); }
    Ok(nodes)
}

struct EdgeInfo {
    tips: BTreeSet<String>,
    is_amf: bool,
    species: Vec<String>,
    genera: Vec<String>,
    unit: String,
    genus_unit: String,
}

/// lsu_placement.py `jplace_tree`: edges in postorder with their oriented tips, labels and units.
fn jplace_tree(tree: &str, meta: &BTreeMap<String, Row>) -> Result<(Vec<usize>, BTreeMap<usize, EdgeInfo>)> {
    let nodes = parse_jplace_tree(tree)?;
    let refs: BTreeSet<String> = meta.keys().cloned().collect();
    let tips: BTreeSet<String> = nodes.iter().filter(|n| n.children.is_empty()).filter_map(|n| n.name.clone()).collect();
    if tips != refs { bail!("jplace tree and metadata taxon sets differ"); }
    let sentinel = meta.iter().filter(|(_, r)| r.get("is_amf").map(String::as_str) == Some("0")).map(|(x, _)| x.clone())
        .min().ok_or_else(|| anyhow!("missing outgroup metadata"))?;
    let mut desc: Vec<BTreeSet<String>> = Vec::with_capacity(nodes.len());
    for n in &nodes {
        let set = if n.children.is_empty() { n.name.iter().cloned().collect() }
                  else { n.children.iter().flat_map(|c| desc[*c].iter().cloned()).collect() };
        desc.push(set);
    }
    let field = |t: &str, k: &str| meta[t].get(k).cloned().unwrap_or_default();
    let mut order = Vec::new();
    let mut edges = BTreeMap::new();
    for (k, n) in nodes.iter().enumerate() {
        let Some(e) = n.edge else { continue };
        let side: BTreeSet<String> = if desc[k].contains(&sentinel) { refs.difference(&desc[k]).cloned().collect() } else { desc[k].clone() };
        let amf = !side.is_empty() && side.iter().all(|t| field(t, "is_amf") == "1");
        let species: Vec<String> = if amf { side.iter().map(|t| field(t, "label")).collect::<BTreeSet<_>>().into_iter().collect() } else { vec![] };
        let genera: Vec<String> = if amf { side.iter().map(|t| field(t, "genus")).collect::<BTreeSet<_>>().into_iter().collect() } else { vec![] };
        if edges.insert(e, EdgeInfo { tips: side, is_amf: amf, species, genera, unit: format!("E{e}"), genus_unit: format!("E{e}") }).is_some() {
            bail!("edge {e} appears twice in the jplace tree");
        }
        order.push(e);
    }
    // collapse nested pure-label components onto their largest containing edge
    for genus in [false, true] {
        let mut pure: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for e in &order {
            let labels = if genus { &edges[e].genera } else { &edges[e].species };
            if labels.len() == 1 { pure.entry(labels[0].clone()).or_default().push(*e); }
        }
        for components in pure.values() {
            for e in components {
                let ds = &edges[e].tips;
                let p = components.iter().filter(|p| ds.is_subset(&edges[*p].tips))
                    .max_by(|a, b| edges[*a].tips.len().cmp(&edges[*b].tips.len()).then(b.cmp(a)))
                    .copied().unwrap_or(*e);
                let info = edges.get_mut(e).unwrap();
                if genus { info.genus_unit = format!("GC{p}"); } else { info.unit = format!("SC{p}"); }
            }
        }
    }
    Ok((order, edges))
}

struct Placement { edge: usize, lwr: PyNum, pendant: PyNum, distal: PyNum }

/// lsu_placement.py `credible_edges`: the smallest highest-weight set reaching `mass`, plus ties.
fn credible_edges(ps: &[Placement], mass: f64) -> Result<Vec<&Placement>> {
    if !(mass > 0.0 && mass <= 1.0) { bail!("credible mass must be in (0,1]"); }
    let mut sorted: Vec<&Placement> = ps.iter().collect();
    sorted.sort_by(|a, b| b.lwr.f().partial_cmp(&a.lwr.f()).unwrap().then(a.edge.cmp(&b.edge)));
    let total = fsum(sorted.iter().map(|p| p.lwr.f()));
    if !isclose(total, 1.0, 1e-9, 1e-6) { bail!("placement weights do not sum to 1: {total}"); }
    let (mut chosen, mut acc, mut cut) = (Vec::new(), 0.0f64, None::<f64>);
    for p in sorted {
        let w = p.lwr.f();
        if let Some(c) = cut { if !isclose(w, c, 1e-12, 0.0) { break; } }
        chosen.push(p);
        acc += w;
        if cut.is_none() && acc >= mass { cut = Some(w); }
    }
    Ok(chosen)
}

/// One row per query, sorted by accession, with `PLACEMENT_HEADERS` as columns.
pub fn summarize_placements(jplace: &Path, metadata: &Path, queries: &Path, aligned: &Path, mass: f64) -> Result<Vec<Row>> {
    Ok(summarize_with_ties(jplace, metadata, queries, aligned, mass)?.0)
}

/// Edges tied for the top like_weight_ratio, per query. The frozen `best_edge` is the first of them in
/// EPA-ng's list, and EPA-ng does not list tied placements in the same order on every run.
pub type TopTies = BTreeMap<String, BTreeSet<usize>>;

pub fn summarize_with_ties(jplace: &Path, metadata: &Path, queries: &Path, aligned: &Path, mass: f64) -> Result<(Vec<Row>, TopTies)> {
    let meta: BTreeMap<String, Row> = read_table(metadata)?.rows.into_iter()
        .map(|r| Ok((row_value(&r, "ref_id")?.to_owned(), r))).collect::<Result<_>>()?;
    let data: Value = serde_json::from_str(&fs::read_to_string(jplace).with_context(|| format!("reading {}", jplace.display()))?)?;
    let (order, edges) = jplace_tree(data.get("tree").and_then(Value::as_str).ok_or_else(|| anyhow!("jplace has no tree"))?, &meta)?;
    let raw = read_fasta(queries)?;
    let aligned = read_fasta(aligned)?;
    if raw.keys().ne(aligned.keys()) { bail!("aligned and raw query IDs differ"); }
    let fields: Vec<&str> = data.get("fields").and_then(Value::as_array).ok_or_else(|| anyhow!("jplace has no fields"))?
        .iter().map(|f| f.as_str().ok_or_else(|| anyhow!("jplace field name is not text"))).collect::<Result<_>>()?;
    let col = |name: &str| fields.iter().position(|f| *f == name).ok_or_else(|| anyhow!("jplace has no {name} field"));
    let (ce, cw, cp, cd) = (col("edge_num")?, col("like_weight_ratio")?, col("pendant_length")?, col("distal_length")?);

    let mut rows: Vec<Row> = Vec::new();
    let mut ties: TopTies = BTreeMap::new();
    let mut seen: BTreeSet<String> = BTreeSet::new();
    for rec in data.get("placements").and_then(Value::as_array).ok_or_else(|| anyhow!("jplace has no placements"))? {
        let mut ps = Vec::new();
        for p in rec.get("p").and_then(Value::as_array).ok_or_else(|| anyhow!("placement without p"))? {
            let p = p.as_array().ok_or_else(|| anyhow!("placement entry is not a list"))?;
            let edge = match PyNum::from(&p[ce])? { PyNum::Int(i) => usize::try_from(i)?, PyNum::Float(x) => x as usize };
            let lwr = PyNum::from(&p[cw])?;
            if !edges.contains_key(&edge) { bail!("unknown edge ID {edge}"); }
            if !lwr.f().is_finite() || !(0.0..=1.0).contains(&lwr.f()) { bail!("invalid like_weight_ratio"); }
            ps.push(Placement { edge, lwr, pendant: PyNum::from(&p[cp])?, distal: PyNum::from(&p[cd])? });
        }
        if ps.iter().map(|p| p.edge).collect::<BTreeSet<_>>().len() != ps.len() { bail!("duplicate edge for a query"); }
        let names: Vec<String> = match rec.get("n").and_then(Value::as_array).filter(|n| !n.is_empty()) {
            Some(n) => n.iter().map(|v| v.as_str().map(str::to_owned).ok_or_else(|| anyhow!("query name is not text"))).collect::<Result<_>>()?,
            None => rec.get("nm").and_then(Value::as_array).ok_or_else(|| anyhow!("placement without names"))?
                .iter().map(|v| v.get(0).and_then(Value::as_str).map(str::to_owned).ok_or_else(|| anyhow!("bad nm entry"))).collect::<Result<_>>()?,
        };
        let cs = credible_edges(&ps, mass)?;
        let mut best = &ps[0];
        for p in &ps[1..] { if p.lwr.f() > best.lwr.f() { best = p; } }
        let (mut species_mass, mut genus_mass): (Vec<(String, f64)>, Vec<(String, f64)>) = (Vec::new(), Vec::new());
        let mut amf_mass = 0.0f64;
        let add = |v: &mut Vec<(String, f64)>, unit: &str, w: f64| match v.iter_mut().find(|(u, _)| u == unit) {
            Some((_, m)) => *m += w,
            None => v.push((unit.to_owned(), 0.0 + w)),
        };
        for p in &ps {
            let info = &edges[&p.edge];
            let w = p.lwr.f();
            if info.is_amf { amf_mass += w; }
            if info.species.len() == 1 { add(&mut species_mass, &info.unit, w); }
            if info.genera.len() == 1 { add(&mut genus_mass, &info.genus_unit, w); }
        }
        let first_max = |v: &[(String, f64)]| v.iter().fold(None::<&(String, f64)>, |b, x| match b { Some(b) if x.1 <= b.1 => Some(b), _ => Some(x) })
            .map(|(u, m)| (u.clone(), *m)).unwrap_or_default();
        let (sp_unit, maxsp) = first_max(&species_mass);
        let (gn_unit, maxgn) = first_max(&genus_mass);
        let label_of = |unit: &str, genus: bool| order.iter().map(|e| &edges[e])
            .find(|d| if genus { d.genus_unit == unit && d.genera.len() == 1 } else { d.unit == unit && d.species.len() == 1 })
            .map(|d| if genus { d.genera[0].clone() } else { d.species[0].clone() }).unwrap_or_default();
        let best_species = label_of(&sp_unit, false);
        let best_genus = label_of(&gn_unit, true);
        let rank = if maxsp >= mass { "species_clade" } else if maxgn >= mass { "genus_clade" } else if amf_mass >= mass { "AMF_broad" } else { "unresolved" };
        let spp: BTreeSet<&String> = cs.iter().flat_map(|p| edges[&p.edge].species.iter()).collect();
        let gns: BTreeSet<&String> = cs.iter().flat_map(|p| edges[&p.edge].genera.iter()).collect();
        let units: BTreeSet<&String> = cs.iter().map(|p| &edges[&p.edge].unit).collect();
        let join = |s: &BTreeSet<&String>, sep: &str| s.iter().map(|x| x.as_str()).collect::<Vec<_>>().join(sep);
        let lwr_sum = fsum(ps.iter().map(|p| p.lwr.f()));
        let entropy = -fsum(ps.iter().filter(|p| p.lwr.f() > 0.0).map(|p| p.lwr.f() * p.lwr.f().ln()));
        let retained = fsum(cs.iter().map(|p| p.lwr.f()));
        let top: BTreeSet<usize> = ps.iter().filter(|p| p.lwr.f() == best.lwr.f()).map(|p| p.edge).collect();
        for x in names {
            if !seen.insert(x.clone()) || !raw.contains_key(&x) { bail!("duplicate or unknown placement query {x}"); }
            ties.insert(x.clone(), top.clone());
            let count = |s: &str, skip: &str| s.chars().filter(|c| !skip.contains(*c)).count();
            let aligned_bases = count(&aligned[&x], "-.");
            let raw_bases = count(&raw[&x], "-.");
            if aligned_bases > raw_bases { bail!("alignment adds query nucleotides for {x}"); }
            let retention = if raw_bases > 0 { aligned_bases as f64 / raw_bases as f64 } else { 0.0 };
            let usable = retention >= 0.9 && aligned_bases > 0;
            let values: [String; 29] = [
                x.clone(), aligned_bases.to_string(), raw_bases.to_string(), py_float(retention), (raw_bases - aligned_bases).to_string(),
                count(&aligned[&x], "ACGT-.").to_string(), if usable { "1" } else { "0" }.into(), best.edge.to_string(), best.lwr.repr(),
                best.pendant.repr(), best.distal.repr(), py_float(lwr_sum), py_float(entropy), py_float(mass), py_float(retained),
                cs.len().to_string(), cs.iter().map(|p| p.edge.to_string()).collect::<Vec<_>>().join(","), join(&units, ","),
                py_float(amf_mass), best_species.clone(), sp_unit.clone(), py_float(maxsp), best_genus.clone(), gn_unit.clone(),
                py_float(maxgn), rank.into(), join(&spp, ";"), join(&gns, ";"), INTERPRETATION.into(),
            ];
            rows.push(PLACEMENT_HEADERS.iter().map(|h| (*h).to_owned()).zip(values).collect());
        }
    }
    if seen.len() != raw.len() { bail!("not every query has a placement ({} of {})", seen.len(), raw.len()); }
    rows.sort_by(|a, b| a.get("accession").cmp(&b.get("accession")));
    Ok((rows, ties))
}

/// How a placement table compares with an expected one, allowing for tied best edges.
#[derive(Clone, Debug, PartialEq)]
pub enum TableMatch {
    /// byte for byte
    Exact,
    /// every cell equal except `best_edge` in these rows, where the expected edge is one of the edges
    /// tied for the top weight in this run
    UpToTiedBestEdges(Vec<String>),
    /// any other difference: (rows differing, first differences)
    Different(usize, Vec<String>),
}

pub fn compare_up_to_ties(rows: &[Row], ties: &TopTies, expected: &Path) -> Result<TableMatch> {
    if placement_tsv(rows) == fs::read(expected)? { return Ok(TableMatch::Exact); }
    let exp = read_table(expected)?;
    let mine: BTreeMap<&str, &Row> = rows.iter().map(|r| (r.get("accession").map(String::as_str).unwrap_or(""), r)).collect();
    let (mut tied, mut other, mut first) = (Vec::new(), BTreeSet::new(), Vec::new());
    if exp.headers.iter().map(String::as_str).ne(PLACEMENT_HEADERS) { bail!("{} does not have the placement table's columns", expected.display()); }
    if exp.rows.len() != rows.len() { other.insert("(row count)".to_owned()); }
    for e in &exp.rows {
        let x = row_value(e, "accession")?;
        let Some(m) = mine.get(x) else { other.insert(x.to_owned()); continue };
        for h in PLACEMENT_HEADERS {
            let (a, b) = (m.get(h).map(String::as_str).unwrap_or(""), e.get(h).map(String::as_str).unwrap_or(""));
            if a == b { continue; }
            let tie_ok = h == "best_edge" && b.parse::<usize>().ok().is_some_and(|edge| ties.get(x).is_some_and(|t| t.contains(&edge) && t.len() > 1));
            if tie_ok { tied.push(x.to_owned()); continue; }
            other.insert(x.to_owned());
            if first.len() < 6 { first.push(format!("{x} [{h}]: got {a:?}; expected {b:?}")); }
        }
    }
    Ok(if other.is_empty() { TableMatch::UpToTiedBestEdges(tied) } else { TableMatch::Different(other.len(), first) })
}

/// The rows as the frozen code writes them (csv.DictWriter, tab-separated, '\n').
pub fn placement_tsv(rows: &[Row]) -> Vec<u8> {
    let headers: Vec<String> = PLACEMENT_HEADERS.iter().map(|h| (*h).to_owned()).collect();
    serialize_table(&headers, rows)
}

#[derive(Clone, Debug)]
pub struct PlacementReport {
    pub rounds: Vec<(String, TableReport)>,
    pub available: usize,
}

impl PlacementReport {
    pub fn passed(&self) -> bool {
        !self.rounds.is_empty() && self.rounds.iter().all(|(_, t)| t.byte_identical && t.differing_cells == 0)
    }
    pub fn summary(&self) -> String {
        let mut s = String::new();
        let ok = self.rounds.iter().filter(|(_, t)| t.byte_identical && t.differing_cells == 0).count();
        for (key, t) in &self.rounds {
            if t.byte_identical && t.differing_cells == 0 { continue; }
            s.push_str(&format!("{key}: {} rows; byte-identical={}; differing cells={}\n", t.rows, t.byte_identical, t.differing_cells));
            for d in &t.first_differences { s.push_str("  "); s.push_str(d); s.push('\n'); }
        }
        s.push_str(&format!("placement summaries: {ok}/{} rounds byte-identical to the published tables ({} of the 51 rounds with their own placement were in the runs folder)\n",
            self.rounds.len(), self.available));
        s.push_str(if self.passed() { "PLACEMENT SUMMARY (lsu_placement.py summarize): PASS\n" } else { "PLACEMENT SUMMARY: FAIL (see differences above)\n" });
        s
    }
}

/// Summarize every round of the LSU-removal test found in an M4 runs folder (the --out of
/// external_holdout.py: <key>/01_reference/metadata.tsv, <key>/03_queries/, <key>/04_epa/epa_result.jplace)
/// and compare it with the round's published table in the benchmark's rounds/ (checked against
/// rounds/MANIFEST.json). The jplace files are large and stay local.
pub fn run_placement_parity(benchmark_dir: impl AsRef<Path>, runs_dir: impl AsRef<Path>, only: Option<&str>) -> Result<PlacementReport> {
    let m4 = benchmark_dir.as_ref().join("benchmarks/m4_external_holdout");
    let manifest_path = m4.join("rounds/MANIFEST.json");
    let manifest: Value = serde_json::from_str(&fs::read_to_string(&manifest_path)?)?;
    let runs = runs_dir.as_ref();
    let mut rounds = Vec::new();
    let mut available = 0;
    for (key, r) in manifest.get("rounds").and_then(Value::as_object).ok_or_else(|| anyhow!("rounds/MANIFEST.json has no rounds"))? {
        if r.get("own_placement").and_then(Value::as_bool) != Some(true) { continue; }
        let fold = runs.join(key);
        let jplace = fold.join("04_epa/epa_result.jplace");
        if !jplace.exists() { continue; }
        available += 1;
        if only.is_some_and(|o| !key.starts_with(o)) { continue; }
        let rel = "05_summary/lsu_placements.copies.tsv";
        let expected = m4.join("rounds").join(key).join(rel);
        let want = r.pointer(&format!("/files/{}", rel.replace('/', "~1"))).and_then(Value::as_str).ok_or_else(|| anyhow!("{key}: no {rel} in the manifest"))?;
        let got = sha256_file(&expected)?;
        if got != want { bail!("{key}/{rel} hash mismatch: {got}"); }
        let rows = summarize_placements(&jplace, &fold.join("01_reference/metadata.tsv"), &fold.join("03_queries/queries.fasta"),
                                        &fold.join("03_queries/queries.aligned.fasta"), 0.95)?;
        rounds.push((key.clone(), compare_table(&expected, &rows, key)?));
    }
    if rounds.is_empty() { bail!("no round of the LSU-removal test with a jplace in {}", runs.display()); }
    Ok(PlacementReport { rounds, available })
}

#[cfg(test)]
pub(crate) mod tests_support {
    pub fn fsum(v: Vec<f64>) -> f64 { super::fsum(v) }
}

// ---------------------------------------------------------------------------------------------
// The front half: MAFFT and EPA-ng with the frozen arguments (lsu_placement.py `build`, query part).

/// The frozen V18 reference of the paper (frozen_m2b/reference): alignment, RAxML-NG tree and model,
/// metadata. `place` refuses another reference unless asked, so that a run is the paper's method.
pub const REFERENCE_PINS: [(&str, &str); 4] = [
    ("reference.aligned.fasta", "1cd6426557f842d9e0085b95f02f3156406d43fb71f840a2c6891bece769ce33"),
    ("fit.raxml.bestTree", "d3da113b83e3d185ebfe10fd1958aa121e401c0dd035c74350ef607a74abf11b"),
    ("fit.raxml.bestModel", "2e7a4dd93fd75c60586bad52a2fd4e4c5a73043289e89ad87f686df21e2a6194"),
    ("reference_metadata.tsv", "5c38fd7fddd3896492644503359917d5fbd1ac3d4c818c0b4c7fc8a49ff5591d"),
];
/// The 903 culture LSU windows (culture reference v4, LSU_LROR-FLR2.fasta) and the frozen table they give.
pub const FROZEN_QUERIES_SHA: &str = "18679043a974d1bff4fdf798c0f74655d7444e1b0916a5c0c8356c796a9cca2f";
pub const FROZEN_TABLE_SHA: &str = "1de103b6ae45f0372dda9ab0799230c4d1243cb98add943e5bab67785e5e910d";

pub struct PlaceOptions {
    pub reference_dir: std::path::PathBuf,
    pub queries: std::path::PathBuf,
    pub out: std::path::PathBuf,
    pub mafft: String,
    pub epa_ng: String,
    /// MAFFT --thread; the frozen build used 4
    pub mafft_threads: usize,
    /// EPA-ng --threads; the frozen build used 8. Other counts give the same calls but can break exact
    /// ties between edges differently, so only 8 reproduces the frozen table byte for byte.
    pub epa_threads: usize,
    pub mass: f64,
    pub allow_other_reference: bool,
    /// Table to compare the result with. When the queries are the frozen 903 culture windows on the
    /// pinned reference at the frozen settings, the caller passes the frozen M2b table here.
    pub compare_with: Option<std::path::PathBuf>,
}

#[derive(Clone, Debug)]
pub struct PlaceReport {
    pub queries: usize,
    pub table_sha256: String,
    /// true when the queries are the frozen 903 culture windows on the pinned reference at the frozen settings
    pub frozen_inputs: bool,
    /// the comparison with `compare_with`, when one was given
    pub comparison: Option<TableMatch>,
}

/// FASTA in file order, as Biopython reads it: ID = first word of the header, sequence uppercased.
pub fn read_fasta_ordered(path: &Path) -> Result<Vec<(String, String)>> {
    let text = fs::read_to_string(path).with_context(|| format!("reading {}", path.display()))?;
    let mut out: Vec<(String, String)> = Vec::new();
    let mut seen = BTreeSet::new();
    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        if let Some(head) = line.strip_prefix('>') {
            let id = head.split_whitespace().next().ok_or_else(|| anyhow!("empty FASTA header in {}", path.display()))?.to_owned();
            if !seen.insert(id.clone()) { bail!("duplicate FASTA ID {id} in {}", path.display()); }
            out.push((id, String::new()));
        } else if let Some(last) = out.last_mut() {
            last.1.extend(line.chars().filter(|c| !c.is_whitespace()).map(|c| c.to_ascii_uppercase()));
        }
    }
    if out.is_empty() { bail!("empty FASTA: {}", path.display()); }
    Ok(out)
}

fn resolve_tool(name: &str) -> Result<std::path::PathBuf> {
    let p = Path::new(name);
    if p.components().count() > 1 { return Ok(p.to_path_buf()); }
    std::env::var_os("PATH").and_then(|paths| std::env::split_paths(&paths).map(|d| d.join(name)).find(|c| c.is_file()))
        .ok_or_else(|| anyhow!("{name} is not on PATH"))
}

fn run_logged(cmd: &[String], log: &Path, stdout: Option<&Path>) -> Result<()> {
    let err = fs::File::create(log)?;
    let mut c = std::process::Command::new(&cmd[0]);
    c.args(&cmd[1..]);
    match stdout {
        Some(path) => { c.stdout(fs::File::create(path)?).stderr(err); }
        None => { c.stdout(err.try_clone()?).stderr(err); }
    }
    let status = c.status().with_context(|| format!("starting {}", cmd[0]))?;
    if !status.success() { bail!("{} failed ({status}); see {}", cmd[0], log.display()); }
    Ok(())
}

fn tool_version(tool: &Path) -> String {
    std::process::Command::new(tool).arg("--version").output().ok()
        .map(|o| format!("{}{}", String::from_utf8_lossy(&o.stdout), String::from_utf8_lossy(&o.stderr)))
        .and_then(|t| t.lines().map(str::trim).find(|l| !l.is_empty()).map(str::to_owned)).unwrap_or_default()
}

/// Place LSU windows on the pinned V18 reference: MAFFT --addfragments --keeplength onto the
/// reference alignment, EPA-ng on the RAxML-NG tree and model, then the placement summary.
/// Writes into `out` (which must be new or empty): queries.fasta, combined.aligned.fasta,
/// queries.aligned.fasta, epa/, logs/, lsu_placements.copies.tsv and run.json.
pub fn run_place(o: &PlaceOptions) -> Result<PlaceReport> {
    let rd = &o.reference_dir;
    let mut pinned_ok = true;
    for (name, sha) in REFERENCE_PINS {
        let got = sha256_file(&rd.join(name))?;
        if got != sha {
            if !o.allow_other_reference { bail!("{} is not the frozen V18 reference ({got}); pass --allow-other-reference to use it anyway", rd.join(name).display()); }
            pinned_ok = false;
        }
    }
    if o.out.exists() && fs::read_dir(&o.out)?.next().is_some() { bail!("{} must be new or empty; earlier results are preserved", o.out.display()); }
    let (logs, epa) = (o.out.join("logs"), o.out.join("epa"));
    fs::create_dir_all(&logs)?;
    fs::create_dir_all(&epa)?;
    let mafft = resolve_tool(&o.mafft)?;
    let epa_ng = resolve_tool(&o.epa_ng)?;
    let queries = o.out.join("queries.fasta");
    fs::copy(&o.queries, &queries)?;
    let qseq = read_fasta_ordered(&queries)?;
    let refs = read_fasta_ordered(&rd.join("reference.aligned.fasta"))?;
    let s = |p: &Path| p.display().to_string();
    let mut commands: Vec<Vec<String>> = Vec::new();

    let combined_path = o.out.join("combined.aligned.fasta");
    let cmd: Vec<String> = vec![s(&mafft), "--auto".into(), "--addfragments".into(), s(&queries), "--keeplength".into(), "--mapout".into(),
                                "--thread".into(), o.mafft_threads.to_string(), s(&rd.join("reference.aligned.fasta"))];
    run_logged(&cmd, &logs.join("query_mafft.log"), Some(&combined_path))?;
    commands.push(cmd);
    let combined: BTreeMap<String, String> = read_fasta_ordered(&combined_path)?.into_iter().collect();
    let expected: BTreeSet<&String> = refs.iter().map(|(x, _)| x).chain(qseq.iter().map(|(x, _)| x)).collect();
    if combined.len() != expected.len() || combined.keys().any(|x| !expected.contains(x))
        || refs.iter().any(|(x, seq)| combined.get(x) != Some(seq)) {
        bail!("MAFFT changed the reference alignment or the IDs differ");
    }
    let aligned_path = o.out.join("queries.aligned.fasta");
    let mut aligned = String::new();
    for (x, _) in &qseq { aligned.push('>'); aligned.push_str(x); aligned.push('\n'); aligned.push_str(&combined[x]); aligned.push('\n'); }
    fs::write(&aligned_path, aligned)?;

    let cmd: Vec<String> = vec![s(&epa_ng), "--ref-msa".into(), s(&rd.join("reference.aligned.fasta")), "--tree".into(), s(&rd.join("fit.raxml.bestTree")),
                                "--query".into(), s(&aligned_path), "--model".into(), s(&rd.join("fit.raxml.bestModel")), "--outdir".into(), s(&epa),
                                "--threads".into(), o.epa_threads.to_string(), "--no-heur".into(), "--no-pre-mask".into(),
                                "--filter-min-lwr".into(), "0".into(), "--filter-max".into(), "2000".into(), "--precision".into(), "12".into()];
    run_logged(&cmd, &logs.join("epa.log"), None)?;
    commands.push(cmd);

    let (rows, ties) = summarize_with_ties(&epa.join("epa_result.jplace"), &rd.join("reference_metadata.tsv"), &queries, &aligned_path, o.mass)?;
    let table = o.out.join("lsu_placements.copies.tsv");
    fs::write(&table, placement_tsv(&rows))?;
    let table_sha256 = sha256_file(&table)?;
    let frozen_inputs = pinned_ok && sha256_file(&queries)? == FROZEN_QUERIES_SHA
        && o.mafft_threads == 4 && o.epa_threads == 8 && o.mass == 0.95;
    let comparison = match &o.compare_with { Some(t) => Some(compare_up_to_ties(&rows, &ties, t)?), None => None };
    if frozen_inputs && o.compare_with.is_none() && table_sha256 != FROZEN_TABLE_SHA {
        bail!("the frozen 903 windows gave a table other than the frozen one; pass the frozen table to see how it differs");
    }
    let comparison_text = comparison.as_ref().map(|c| match c {
        TableMatch::Exact => "byte-identical".to_owned(),
        TableMatch::UpToTiedBestEdges(r) => format!("identical except a tied best edge in {} rows: {}", r.len(), r.join(";")),
        TableMatch::Different(n, _) => format!("{n} rows differ beyond tied best edges"),
    });
    let tools = serde_json::json!({
        "mafft": {"path": s(&mafft), "version": tool_version(&mafft), "sha256": sha256_file(&mafft).unwrap_or_default()},
        "epa_ng": {"path": s(&epa_ng), "version": tool_version(&epa_ng), "sha256": sha256_file(&epa_ng).unwrap_or_default()},
    });
    let run = serde_json::json!({
        "commands": commands, "tools": tools, "reference_pinned": pinned_ok,
        "inputs": {"queries.fasta": sha256_file(&queries)?, "reference": REFERENCE_PINS.iter().map(|(n, _)| (n.to_string(), sha256_file(&rd.join(n)).unwrap_or_default())).collect::<BTreeMap<_, _>>()},
        "outputs": {"lsu_placements.copies.tsv": table_sha256, "queries": rows.len()},
        "frozen_inputs": frozen_inputs,
        "comparison": {"with": o.compare_with.as_ref().map(|p| s(p)), "result": comparison_text},
    });
    fs::write(o.out.join("run.json"), serde_json::to_string_pretty(&run)? + "\n")?;
    Ok(PlaceReport { queries: rows.len(), table_sha256, frozen_inputs, comparison })
}
