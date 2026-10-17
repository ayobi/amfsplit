//! The Stage 1 run: reads -> per-read records, region and window FASTAs, summary.
//!
//! A port of `main()` in amfsplit.py v0.3.1. Output files, columns, FASTA headers and the
//! summary JSON are byte-identical to the prototype's (checked by the parity tests in
//! `parity/stage1`); the console report on stderr is informational only.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};

use crate::anchors::{
    AnchorRow, ITSX_REGIONS, Regions, confidence_of, itsx_regions, load_anchors,
    orientation_from_anchors,
};
use crate::concat::{ConcatParams, concatemer_screen};
use crate::panel::{Window, load_panel};
use crate::pycompat::{
    Json, Num, OrderedCounter, csv_row, float_repr, median, py_round, py_slice, quantiles,
};
use crate::search::{find_site, max_edits, region_bounds};
use crate::seqio::{SeqReader, revcomp};
use crate::window::{USABLE, WindowResult, cut_window};

/// The version of the Stage 1 algorithm this port reproduces; written to `summary.json`
/// exactly as the prototype writes it.
pub const STAGE1_ALGORITHM_VERSION: &str = "0.3.1";

/// Options of the Stage 1 run, named and defaulted as in the prototype's command line.
#[derive(Clone, Debug)]
pub struct Params {
    pub reads: PathBuf,
    pub anchors: Option<PathBuf>,
    pub panel: PathBuf,
    pub out: String,
    pub max_err_frac: f64,
    pub keep_primers: bool,
    pub open_ends: bool,
    pub only_anchored: bool,
    pub region_margin: i64,
    pub no_region_limit: bool,
    pub no_concat_screen: bool,
    pub concat_min_primers: i64,
    pub concat_min_sep: i64,
    pub concat_sim: f64,
    pub concat_strong: i64,
    pub concat_tol: f64,
    pub concat_err_frac: Option<f64>,
    pub keep_concatemers: bool,
}

impl Params {
    /// The prototype's defaults for everything but the paths.
    pub fn new(reads: PathBuf, anchors: Option<PathBuf>, panel: PathBuf, out: String) -> Self {
        Params {
            reads,
            anchors,
            panel,
            out,
            max_err_frac: 0.15,
            keep_primers: false,
            open_ends: false,
            only_anchored: false,
            region_margin: 50,
            no_region_limit: false,
            no_concat_screen: false,
            concat_min_primers: 2,
            concat_min_sep: 150,
            concat_sim: 0.40,
            concat_strong: 3,
            concat_tol: 0.05,
            concat_err_frac: None,
            keep_concatemers: false,
        }
    }

    pub fn output_path(&self, suffix: &str) -> PathBuf {
        PathBuf::from(format!("{}.{}", self.out, suffix))
    }
}

/// What a run produced, for the caller's report and provenance.
#[derive(Clone, Debug)]
pub struct RunReport {
    pub summary: Json,
    pub outputs: Vec<PathBuf>,
    pub console: String,
}

/// One records.tsv row under construction.
struct Rec<'a> {
    vals: Vec<String>,
    idx: &'a HashMap<&'a str, usize>,
}

impl Rec<'_> {
    fn set(&mut self, col: &str, v: String) {
        self.vals[self.idx[col]] = v;
    }
}

/// `orient_by_primers`: orientation of an unanchored read from the majority of primer sites.
fn orient_by_primers(seq: &[u8], windows: &[Window], frac: f64) -> Result<char> {
    let rc = revcomp(seq);
    let (mut plus, mut minus) = (0, 0);
    for w in windows {
        for site in [w.fwd.clone(), w.rev_site()] {
            let k = max_edits(site.len(), frac);
            if find_site(&site, seq, k)?.is_some() {
                plus += 1;
            }
            if find_site(&site, &rc, k)?.is_some() {
                minus += 1;
            }
        }
    }
    Ok(if plus == 0 && minus == 0 {
        '?'
    } else if plus >= minus {
        '+'
    } else {
        '-'
    })
}

/// `dist(v)` of the prototype's summary.
fn length_dist(v: &[i64]) -> Json {
    if v.is_empty() {
        return Json::Null;
    }
    let (q1, q3) = if v.len() > 1 {
        let q = quantiles(v, 4);
        (Json::Float(q[0]), Json::Float(q[2]))
    } else {
        (Json::Int(v[0]), Json::Int(v[0]))
    };
    Json::obj(vec![
        ("n", Json::Int(v.len() as i64)),
        ("median", Json::from(median(v))),
        ("q1", q1),
        ("q3", q3),
        ("min", Json::Int(*v.iter().min().expect("non-empty"))),
        ("max", Json::Int(*v.iter().max().expect("non-empty"))),
    ])
}

fn num_text(n: Num) -> String {
    match n {
        Num::Int(i) => i.to_string(),
        Num::Float(f) => float_repr(f),
    }
}

fn opt_text(v: Option<i64>) -> String {
    v.map_or(String::new(), |x| x.to_string())
}

/// Run Stage 1 on one set of reads and write `{out}.records.tsv`, one FASTA per ITSxRust region
/// and per panel window, and `{out}.summary.json`.
pub fn run(p: &Params) -> Result<RunReport> {
    let windows = load_panel(&p.panel)?;
    let anchors: HashMap<String, AnchorRow> = match &p.anchors {
        Some(a) => load_anchors(a)?,
        None => HashMap::new(),
    };

    // output files, opened (and created) up front as the prototype does
    let region_names: Vec<String> = ITSX_REGIONS
        .iter()
        .map(|s| s.to_string())
        .chain(windows.iter().map(|w| w.name.clone()))
        .collect();
    let mut outputs = Vec::new();
    let mut fasta: Vec<BufWriter<File>> = Vec::new();
    for name in &region_names {
        let path = p.output_path(&format!("{name}.fasta"));
        let f = File::create(&path).with_context(|| format!("cannot write {}", path.display()))?;
        fasta.push(BufWriter::new(f));
        outputs.push(path);
    }
    let rec_path = p.output_path("records.tsv");
    let mut rec_fh = BufWriter::new(
        File::create(&rec_path).with_context(|| format!("cannot write {}", rec_path.display()))?,
    );
    outputs.push(rec_path);

    const WIN_FIELDS: [&str; 7] = ["status", "start", "end", "len", "fwd_ed", "rev_ed", "multi"];
    let mut cols: Vec<String> = [
        "read_id",
        "read_len",
        "orient",
        "orient_source",
        "itsx_confidence",
    ]
    .iter()
    .map(|s| s.to_string())
    .collect();
    cols.extend(ITSX_REGIONS.iter().map(|r| format!("{r}_coords")));
    for w in &windows {
        cols.extend(WIN_FIELDS.iter().map(|k| format!("{}_{k}", w.name)));
    }
    cols.extend(
        [
            "windows_ok",
            "flank_consistent",
            "region_limited",
            "concatemer",
            "concat_evidence",
            "concat_kind",
        ]
        .iter()
        .map(|s| s.to_string()),
    );
    let col_index: HashMap<&str, usize> = cols
        .iter()
        .enumerate()
        .map(|(i, c)| (c.as_str(), i))
        .collect();
    rec_fh.write_all(csv_row(&cols).as_bytes())?;

    // primers for the concatemer screen: one entry per distinct site, first name kept; a name
    // given to two sites keeps its first position and its last site (dict semantics)
    let mut by_site: Vec<(Vec<u8>, String)> = Vec::new();
    for w in &windows {
        for (site, name) in [(w.fwd.clone(), &w.fwd_name), (w.rev_site(), &w.rev_name)] {
            if !by_site.iter().any(|(s, _)| *s == site) {
                by_site.push((site, name.clone()));
            }
        }
    }
    let mut screen_primers: Vec<(String, Vec<u8>)> = Vec::new();
    for (site, name) in &by_site {
        match screen_primers.iter_mut().find(|(n, _)| n == name) {
            Some(entry) => entry.1 = site.clone(),
            None => screen_primers.push((name.clone(), site.clone())),
        }
    }
    let concat_frac = p.concat_err_frac.unwrap_or(p.max_err_frac);
    let cparams = ConcatParams {
        frac: concat_frac,
        min_sep: p.concat_min_sep,
        min_primers: p.concat_min_primers,
        tol_frac: p.concat_tol,
        max_dist: p.concat_sim,
        strong_primers: p.concat_strong,
    };

    let (mut n_reads, mut n_processed, mut n_concat, mut n_limited) = (0i64, 0i64, 0i64, 0i64);
    let mut concat_kinds = OrderedCounter::default();
    let mut status_counts: Vec<OrderedCounter> = vec![OrderedCounter::default(); windows.len()];
    let mut lengths: HashMap<String, Vec<i64>> = HashMap::new();
    let mut orient_src = OrderedCounter::default();
    let mut region_counts: HashMap<&str, i64> = HashMap::new();
    let mut inconsistent: Vec<i64> = vec![0; windows.len()];

    for item in SeqReader::open(&p.reads)? {
        let (rid, raw) = item?;
        n_reads += 1;
        let row = anchors.get(&rid);
        if row.is_none() && p.only_anchored {
            orient_src.add("skipped_unanchored", 1);
            continue;
        }
        n_processed += 1;
        let anchored_orient = row.and_then(orientation_from_anchors);
        let source = if anchored_orient.is_some() {
            "itsxrust"
        } else {
            "primers"
        };
        let orient = match anchored_orient {
            Some(o) => o,
            None => orient_by_primers(&raw, &windows, p.max_err_frac)?,
        };
        let source = if orient == '?' { "unresolved" } else { source };
        orient_src.add(source, 1);
        let seq = if orient == '-' { revcomp(&raw) } else { raw };
        let n = seq.len() as i64;

        let mut rec = Rec {
            vals: vec![String::new(); cols.len()],
            idx: &col_index,
        };
        rec.set("read_id", rid.clone());
        rec.set("read_len", n.to_string());
        rec.set("orient", orient.to_string());
        rec.set("orient_source", source.to_string());
        rec.set(
            "itsx_confidence",
            row.map_or("unanchored".to_string(), confidence_of),
        );

        let regions = match row {
            Some(r) => itsx_regions(r, n).with_context(|| format!("read {rid}: bad anchors"))?,
            None => Regions::default(),
        };
        let full = regions.full_its;
        let limited = full.is_some() && !p.no_region_limit;
        rec.set(
            "region_limited",
            if limited { "yes" } else { "no" }.to_string(),
        );
        n_limited += i64::from(limited);

        // concatemer screen: whole read, independent of the region limit
        let screen = if !p.no_concat_screen && orient != '?' {
            concatemer_screen(&seq, &screen_primers, &cparams)
                .with_context(|| format!("read {rid}"))?
        } else {
            Default::default()
        };
        rec.set(
            "concatemer",
            if screen.is_concatemer { "yes" } else { "no" }.to_string(),
        );
        rec.set("concat_evidence", screen.evidence.clone());
        rec.set("concat_kind", screen.kind.to_string());
        if screen.is_concatemer {
            n_concat += 1;
            concat_kinds.add(screen.kind, 1);
        }
        let skip_output = screen.is_concatemer && !p.keep_concatemers;

        for (i, name) in ITSX_REGIONS.iter().enumerate() {
            let coords = regions.get(name);
            rec.set(
                &format!("{name}_coords"),
                coords.map_or(String::new(), |c| format!("{}-{}", c.0, c.1)),
            );
            if let Some((a, b)) = coords
                && b >= a
                && !skip_output
            {
                *region_counts.entry(name).or_insert(0) += 1;
                lengths.entry(name.to_string()).or_default().push(b - a + 1);
                let fh = &mut fasta[i];
                writeln!(fh, ">{rid} region={name} coords={a}-{b} orient={orient}")?;
                fh.write_all(py_slice(&seq, a - 1, b))?;
                fh.write_all(b"\n")?;
            }
        }

        if skip_output {
            // flagged, not dropped: the read stays in records.tsv with every window marked
            for (wi, w) in windows.iter().enumerate() {
                status_counts[wi].add("concatemer", 1);
                rec.set(&format!("{}_status", w.name), "concatemer".into());
            }
            rec.set("windows_ok", "0".into());
            rec.set("flank_consistent", "n/a".into());
            rec_fh.write_all(csv_row(&rec.vals).as_bytes())?;
            continue;
        }

        let mut n_ok = 0;
        let mut flank_ok = true;
        for (wi, w) in windows.iter().enumerate() {
            let res = if orient == '?' {
                WindowResult::unoriented()
            } else {
                let fb = limited.then(|| region_bounds(w.fwd_region, full, n, p.region_margin));
                let rb = limited.then(|| region_bounds(w.rev_region, full, n, p.region_margin));
                cut_window(w, &seq, p.max_err_frac, p.keep_primers, p.open_ends, fb, rb)
                    .with_context(|| format!("read {rid}, window {}", w.name))?
            };
            status_counts[wi].add(res.status, 1);
            let (start, end, len) = match res.span {
                Some((s, e, l)) => (Some(s), Some(e), Some(l)),
                None => (None, None, None),
            };
            rec.set(&format!("{}_status", w.name), res.status.to_string());
            rec.set(&format!("{}_start", w.name), opt_text(start));
            rec.set(&format!("{}_end", w.name), opt_text(end));
            rec.set(&format!("{}_len", w.name), opt_text(len));
            rec.set(&format!("{}_fwd_ed", w.name), opt_text(res.fwd_ed));
            rec.set(&format!("{}_rev_ed", w.name), opt_text(res.rev_ed));
            rec.set(&format!("{}_multi", w.name), opt_text(res.multi));
            if res.usable() {
                let (s, e, l) = res.span.expect("usable windows have coordinates");
                n_ok += 1;
                lengths.entry(w.name.clone()).or_default().push(l);
                let fh = &mut fasta[ITSX_REGIONS.len() + wi];
                writeln!(
                    fh,
                    ">{rid} window={} coords={s}-{e} status={} orient={orient} fwd_ed={} rev_ed={}",
                    w.name,
                    res.status,
                    opt_text(res.fwd_ed),
                    opt_text(res.rev_ed)
                )?;
                fh.write_all(py_slice(&seq, s - 1, e))?;
                fh.write_all(b"\n")?;
                // cross-check primer-derived windows against the HMM anchors
                if let Some(f) = full {
                    if w.name.starts_with("SSU") && e >= f.0 {
                        flank_ok = false;
                        inconsistent[wi] += 1;
                    }
                    if w.name.starts_with("LSU") && s <= f.1 {
                        flank_ok = false;
                        inconsistent[wi] += 1;
                    }
                }
            }
        }
        rec.set("windows_ok", n_ok.to_string());
        rec.set(
            "flank_consistent",
            if flank_ok { "yes" } else { "no" }.to_string(),
        );
        rec_fh.write_all(csv_row(&rec.vals).as_bytes())?;
    }

    for mut fh in fasta {
        fh.flush()?;
    }
    rec_fh.flush()?;

    let denom = n_processed.max(1);
    let pct =
        |count: i64, nd: usize| Json::Float(py_round((100 * count) as f64 / denom as f64, nd));
    let empty: Vec<i64> = Vec::new();
    let len_of = |name: &str| length_dist(lengths.get(name).unwrap_or(&empty));
    let summary = Json::obj(vec![
        ("reads", Json::Int(n_reads)),
        ("reads_processed", Json::Int(n_processed)),
        ("orientation_source", orient_src.to_json()),
        (
            "itsxrust_regions",
            Json::Obj(
                ITSX_REGIONS
                    .iter()
                    .map(|r| {
                        let reads = *region_counts.get(r).unwrap_or(&0);
                        (
                            r.to_string(),
                            Json::obj(vec![("reads", Json::Int(reads)), ("length", len_of(r))]),
                        )
                    })
                    .collect(),
            ),
        ),
        (
            "windows",
            Json::Obj(
                windows
                    .iter()
                    .enumerate()
                    .map(|(wi, w)| {
                        let sc = &status_counts[wi];
                        let usable: i64 = USABLE.iter().map(|k| sc.get(k)).sum();
                        (
                            w.name.clone(),
                            Json::obj(vec![
                                ("target", Json::Str(w.target.clone())),
                                (
                                    "primers",
                                    Json::Str(format!("{}/{}", w.fwd_name, w.rev_name)),
                                ),
                                ("ok", Json::Int(sc.get("ok"))),
                                ("open", Json::Int(sc.get("open_5p") + sc.get("open_3p"))),
                                ("usable", Json::Int(usable)),
                                ("ok_pct", pct(sc.get("ok"), 1)),
                                ("usable_pct", pct(usable, 1)),
                                ("status", sc.to_json()),
                                ("length", len_of(&w.name)),
                                ("crosses_its_anchor", Json::Int(inconsistent[wi])),
                            ]),
                        )
                    })
                    .collect(),
            ),
        ),
        (
            "concatemers",
            Json::obj(vec![
                ("flagged", Json::Int(n_concat)),
                ("pct", pct(n_concat, 2)),
                ("full", Json::Int(concat_kinds.get("full"))),
                ("partial", Json::Int(concat_kinds.get("partial"))),
                ("kept_in_fasta", Json::Bool(p.keep_concatemers)),
            ]),
        ),
        ("region_limited_reads", Json::Int(n_limited)),
        (
            "params",
            Json::obj(vec![
                ("version", Json::Str(STAGE1_ALGORITHM_VERSION.into())),
                ("max_err_frac", Json::Float(p.max_err_frac)),
                ("keep_primers", Json::Bool(p.keep_primers)),
                ("open_ends", Json::Bool(p.open_ends)),
                ("only_anchored", Json::Bool(p.only_anchored)),
                ("region_margin", Json::Int(p.region_margin)),
                ("region_limit", Json::Bool(!p.no_region_limit)),
                ("concat_screen", Json::Bool(!p.no_concat_screen)),
                ("concat_min_primers", Json::Int(p.concat_min_primers)),
                ("concat_min_sep", Json::Int(p.concat_min_sep)),
                ("concat_tol", Json::Float(p.concat_tol)),
                ("concat_sim", Json::Float(p.concat_sim)),
                ("concat_strong", Json::Int(p.concat_strong)),
                ("concat_err_frac", Json::Float(concat_frac)),
            ]),
        ),
    ]);
    let sum_path = p.output_path("summary.json");
    std::fs::write(&sum_path, summary.dump_indent1())
        .with_context(|| format!("cannot write {}", sum_path.display()))?;
    outputs.push(sum_path);

    // console report, in the prototype's words
    let mut console = String::new();
    let pct_text =
        |count: i64, nd: usize| float_repr(py_round((100 * count) as f64 / denom as f64, nd));
    console += &format!(
        "{n_reads} reads ({n_processed} processed)  orientation: {}\n",
        orient_src.py_dict_repr()
    );
    console += &format!(
        "  concatemers flagged: {n_concat} ({}%; {} full, {} partial)  region-limited reads: {n_limited}\n",
        pct_text(n_concat, 2),
        concat_kinds.get("full"),
        concat_kinds.get("partial")
    );
    let med = |name: &str| {
        lengths
            .get(name)
            .filter(|v| !v.is_empty())
            .map_or("-".to_string(), |v| num_text(median(v)))
    };
    for r in ITSX_REGIONS {
        console += &format!(
            "  {r:<24} {:>6} reads  median len {}\n",
            region_counts.get(r).unwrap_or(&0),
            med(r)
        );
    }
    for (wi, w) in windows.iter().enumerate() {
        let sc = &status_counts[wi];
        let usable: i64 = USABLE.iter().map(|k| sc.get(k)).sum();
        console += &format!(
            "  {:<24} {usable:>6} usable ({}%; {} open)  median len {}  {}\n",
            w.name,
            pct_text(usable, 1),
            sc.get("open_5p") + sc.get("open_3p"),
            med(&w.name),
            sc.py_dict_repr()
        );
    }
    Ok(RunReport {
        summary,
        outputs,
        console,
    })
}

/// Fail early on an output prefix whose directory does not exist.
pub fn check_out_prefix(out: &str) -> Result<()> {
    let parent = Path::new(out)
        .parent()
        .filter(|p| !p.as_os_str().is_empty());
    if let Some(dir) = parent
        && !dir.is_dir()
    {
        bail!("output directory {} does not exist", dir.display());
    }
    Ok(())
}
