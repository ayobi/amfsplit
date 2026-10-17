//! Concatemer screen (`concatemer_screen` of the prototype, v0.3.1).
//!
//! Two copies in one read repeat every primer site with the same spacing:
//!   tandem    each primer occurs at two sites at least `min_sep` apart, and at least
//!             `min_primers` of them share one offset (within max(40 nt, tol_frac of it));
//!   inverted  each primer also occurs on the opposite strand, and at least `min_primers` of
//!             them share one mirror sum (sense position + antisense position + site length).
//! A cluster of fewer than `strong_primers` primers must also pass a primer-free check: the two
//! supposed copies (up to 800 nt) differ by at most `max_dist` edits per base.

use anyhow::Result;

use crate::pycompat::{median, median_trunc, py_slice};
use crate::search::{find_all_sites, find_site, max_edits, norm_ed};
use crate::seqio::revcomp;

#[derive(Clone, Debug)]
pub struct ConcatParams {
    pub frac: f64,
    pub min_sep: i64,
    pub min_primers: i64,
    pub tol_frac: f64,
    pub max_dist: f64,
    pub strong_primers: i64,
}

const MIN_UNIT: i64 = 100;
const TOL_MIN: f64 = 40.0;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ConcatResult {
    pub is_concatemer: bool,
    pub evidence: String,
    /// "full", "partial" or "".
    pub kind: &'static str,
}

/// `_best_cluster`: the largest set of primers whose values agree within
/// max(TOL_MIN, tol_frac * |value|), names sorted.
fn best_cluster(named: &[(String, i64)], tol_frac: f64) -> Vec<String> {
    let mut best: Vec<String> = Vec::new();
    for (_, v0) in named {
        let tol = TOL_MIN.max(tol_frac * v0.abs() as f64);
        let members: Vec<String> = named
            .iter()
            .filter(|(_, v)| ((v - v0).abs() as f64) <= tol)
            .map(|(n, _)| n.clone())
            .collect();
        if members.len() > best.len() {
            best = members;
        }
    }
    best.sort();
    best
}

/// Screen one oriented read. `primers` are (name, rRNA-sense site) in the prototype's order.
pub fn concatemer_screen(
    seq: &[u8],
    primers: &[(String, Vec<u8>)],
    p: &ConcatParams,
) -> Result<ConcatResult> {
    let n = seq.len() as i64;
    let rc = revcomp(seq);
    // name -> (spacing, first site) and name -> (mirror sum, first site), insertion-ordered
    let mut offsets: Vec<(String, (i64, i64))> = Vec::new();
    let mut mirrors: Vec<(String, (i64, i64))> = Vec::new();
    for (name, site) in primers {
        let k = max_edits(site.len(), p.frac);
        let mut sites: Vec<i64> = find_all_sites(site, seq, k, p.min_sep, 3)?
            .into_iter()
            .map(|s| s.0)
            .collect();
        sites.sort_unstable();
        if sites.len() >= 2 {
            offsets.push((name.clone(), (sites[1] - sites[0], sites[0])));
        }
        let anti = find_site(site, &rc, k)?;
        if let (Some(a), Some(&first)) = (anti, sites.first()) {
            mirrors.push((
                name.clone(),
                (first + (n - 1 - a.end) + site.len() as i64, first),
            ));
        }
    }
    let lookup = |v: &[(String, (i64, i64))], name: &str| {
        v.iter().find(|(k, _)| k == name).expect("cluster member").1
    };
    let mut evidence: Vec<String> = Vec::new();
    let mut kinds: Vec<&'static str> = Vec::new();

    if !offsets.is_empty() {
        let values: Vec<(String, i64)> = offsets.iter().map(|(k, v)| (k.clone(), v.0)).collect();
        let cluster = best_cluster(&values, p.tol_frac);
        if cluster.len() as i64 >= p.min_primers {
            let period = median_trunc(
                &cluster
                    .iter()
                    .map(|t| lookup(&offsets, t).0)
                    .collect::<Vec<_>>(),
            );
            let pos = median_trunc(
                &cluster
                    .iter()
                    .map(|t| lookup(&offsets, t).1)
                    .collect::<Vec<_>>(),
            );
            let w = period.min(n - (pos + period)).min(800);
            let strong = cluster.len() as i64 >= p.strong_primers;
            let d = if w >= MIN_UNIT && !strong {
                Some(norm_ed(
                    py_slice(seq, pos, pos + w),
                    py_slice(seq, pos + period, pos + period + w),
                )?)
            } else {
                None
            };
            if strong || d.is_some_and(|d| d <= p.max_dist) {
                evidence.push(format!(
                    "tandem:{}@{}{}",
                    cluster.join(","),
                    period,
                    d.map_or(String::new(), |d| format!("~{d:.2}"))
                ));
                kinds.push(if period as f64 >= 0.4 * n as f64 {
                    "full"
                } else {
                    "partial"
                });
            }
        }
    }
    if !mirrors.is_empty() {
        let values: Vec<(String, i64)> = mirrors.iter().map(|(k, v)| (k.clone(), v.0)).collect();
        let cluster = best_cluster(&values, p.tol_frac);
        if cluster.len() as i64 >= p.min_primers {
            let m = median(
                &cluster
                    .iter()
                    .map(|t| lookup(&mirrors, t).0)
                    .collect::<Vec<_>>(),
            );
            let centre = (m.as_f64() / 2.0).trunc() as i64;
            let w = centre.min(n - centre).min(800);
            let strong = cluster.len() as i64 >= p.strong_primers;
            let d = if w >= MIN_UNIT && !strong {
                let right = revcomp(py_slice(seq, centre, centre + w));
                Some(norm_ed(py_slice(seq, centre - w, centre), &right)?)
            } else {
                None
            };
            if strong || d.is_some_and(|d| d <= p.max_dist) {
                evidence.push(format!(
                    "inverted:{}@{}{}",
                    cluster.join(","),
                    centre,
                    d.map_or(String::new(), |d| format!("~{d:.2}"))
                ));
                let c = centre as f64;
                kinds.push(if 0.35 * n as f64 <= c && c <= 0.65 * n as f64 {
                    "full"
                } else {
                    "partial"
                });
            }
        }
    }
    let kind = if kinds.is_empty() {
        ""
    } else if kinds.contains(&"full") {
        "full"
    } else {
        "partial"
    };
    Ok(ConcatResult {
        is_concatemer: !evidence.is_empty(),
        evidence: evidence.join(";"),
        kind,
    })
}
