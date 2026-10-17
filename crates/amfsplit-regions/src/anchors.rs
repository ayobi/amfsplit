//! ITSxRust per-read anchors (`--anchors-tsv`): ITS sub-regions, flanks and orientation.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

use crate::seqio::open_maybe_gz;

/// One anchors row, keeping only the columns Stage 1 reads. A value is `None` when the column is
/// missing from the file or the row is short (csv.DictReader gives None there), and `Some("")`
/// when the cell is empty.
#[derive(Clone, Debug, Default)]
pub struct AnchorRow {
    pub its1_start: Option<String>,
    pub its1_end: Option<String>,
    pub its2_start: Option<String>,
    pub its2_end: Option<String>,
    pub full_start: Option<String>,
    pub full_end: Option<String>,
    /// None: the file has no confidence column (reported as "unanchored", as the prototype does).
    pub confidence: Option<Option<String>>,
    pub strands: [Option<String>; 4],
}

pub fn load_anchors(path: &Path) -> Result<HashMap<String, AnchorRow>> {
    let mut text = String::new();
    open_maybe_gz(path)?
        .read_to_string(&mut text)
        .with_context(|| format!("{} is not UTF-8 text", path.display()))?;
    let mut lines = text.split('\n').map(|l| l.trim_end_matches('\r'));
    let header: Vec<&str> = match lines.next() {
        Some(h) if !h.is_empty() => h.split('\t').collect(),
        _ => return Ok(HashMap::new()),
    };
    let col = |n: &str| header.iter().position(|h| *h == n);
    let id = col("read_id").with_context(|| format!("{} has no read_id column", path.display()))?;
    let cols = [
        col("its1_start"),
        col("its1_end"),
        col("its2_start"),
        col("its2_end"),
        col("full_start"),
        col("full_end"),
    ];
    let conf = col("confidence");
    let strand_cols = [
        col("ssu_strand"),
        col("s58s_strand"),
        col("s58e_strand"),
        col("lsu_strand"),
    ];
    let mut out = HashMap::new();
    for line in lines {
        if line.is_empty() {
            continue; // csv.DictReader skips blank rows
        }
        let f: Vec<&str> = line.split('\t').collect();
        let get = |c: Option<usize>| c.and_then(|i| f.get(i)).map(|s| s.to_string());
        let rid = match f.get(id) {
            Some(r) => r.to_string(),
            None => continue,
        };
        let row = AnchorRow {
            its1_start: get(cols[0]),
            its1_end: get(cols[1]),
            its2_start: get(cols[2]),
            its2_end: get(cols[3]),
            full_start: get(cols[4]),
            full_end: get(cols[5]),
            confidence: conf.map(|c| f.get(c).map(|s| s.to_string())),
            strands: [
                get(strand_cols[0]),
                get(strand_cols[1]),
                get(strand_cols[2]),
                get(strand_cols[3]),
            ],
        };
        out.insert(rid, row); // a repeated read_id keeps the last row
    }
    Ok(out)
}

/// `_int(v)`: None for a missing or empty cell, otherwise Python's `int(v)`.
fn py_int_opt(v: &Option<String>) -> Result<Option<i64>> {
    match v.as_deref() {
        None | Some("") => Ok(None),
        Some(s) => s
            .trim()
            .parse::<i64>()
            .map(Some)
            .with_context(|| format!("anchor value is not an integer: {s:?}")),
    }
}

/// 1-based inclusive regions in the read's rRNA-sense orientation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Regions {
    pub ssu_flank: Option<(i64, i64)>,
    pub its1: Option<(i64, i64)>,
    pub s58: Option<(i64, i64)>,
    pub its2: Option<(i64, i64)>,
    pub full_its: Option<(i64, i64)>,
    pub lsu_flank: Option<(i64, i64)>,
}

pub const ITSX_REGIONS: [&str; 6] = ["ssu_flank", "its1", "s58", "its2", "full_its", "lsu_flank"];

impl Regions {
    pub fn get(&self, name: &str) -> Option<(i64, i64)> {
        match name {
            "ssu_flank" => self.ssu_flank,
            "its1" => self.its1,
            "s58" => self.s58,
            "its2" => self.its2,
            "full_its" => self.full_its,
            "lsu_flank" => self.lsu_flank,
            _ => None,
        }
    }
}

/// `itsx_regions(row, read_len)` of the prototype.
pub fn itsx_regions(row: &AnchorRow, read_len: i64) -> Result<Regions> {
    let pair = |a: &Option<String>, b: &Option<String>| -> Result<(Option<i64>, Option<i64>)> {
        Ok((py_int_opt(a)?, py_int_opt(b)?))
    };
    let both = |p: (Option<i64>, Option<i64>)| p.0.zip(p.1);
    let i1 = pair(&row.its1_start, &row.its1_end)?;
    let i2 = pair(&row.its2_start, &row.its2_end)?;
    let fl = pair(&row.full_start, &row.full_end)?;
    let its1 = both(i1);
    let its2 = both(i2);
    let full = both(fl);
    Ok(Regions {
        its1,
        its2,
        full_its: full,
        s58: match (its1, its2) {
            (Some(a), Some(b)) => Some((a.1 + 1, b.0 - 1)),
            _ => None,
        },
        ssu_flank: full.filter(|f| f.0 > 1).map(|f| (1, f.0 - 1)),
        lsu_flank: full.filter(|f| f.1 < read_len).map(|f| (f.1 + 1, read_len)),
    })
}

/// First strand among ssu, s58s, s58e, lsu that is "+" or "-".
pub fn orientation_from_anchors(row: &AnchorRow) -> Option<char> {
    row.strands.iter().find_map(|s| match s.as_deref() {
        Some("+") => Some('+'),
        Some("-") => Some('-'),
        _ => None,
    })
}

/// The value written to `itsx_confidence` for an anchored read.
pub fn confidence_of(row: &AnchorRow) -> String {
    match &row.confidence {
        None => "unanchored".into(),
        Some(v) => v.clone().unwrap_or_default(),
    }
}

use std::io::Read as _;
