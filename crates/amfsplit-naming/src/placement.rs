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
        for x in names {
            if !seen.insert(x.clone()) || !raw.contains_key(&x) { bail!("duplicate or unknown placement query {x}"); }
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
    Ok(rows)
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
