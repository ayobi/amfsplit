//! Primer-window panel: one row per legacy barcode window (tab-separated, `#` comments).
//!
//! Required columns: window, target, fwd_name, fwd_seq, rev_name, rev_seq, min_len, max_len.
//! Optional: fwd_region, rev_region (ssu, lsu or any); when absent or empty, the region follows
//! the window name (SSU_*: ssu, LSU_*: lsu, AMPLICON_*: ssu then lsu, anything else: any).

use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::pycompat::py_strip;
use crate::seqio::{open_maybe_gz, revcomp};

/// Where a primer is searched, relative to the ITSxRust full-ITS anchors.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Region {
    Ssu,
    Lsu,
    Any,
}

impl Region {
    fn parse(s: &str) -> Option<Region> {
        match s {
            "ssu" => Some(Region::Ssu),
            "lsu" => Some(Region::Lsu),
            "any" => Some(Region::Any),
            _ => None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Window {
    pub name: String,
    pub target: String,
    pub fwd_name: String,
    pub fwd: Vec<u8>,
    pub rev_name: String,
    pub rev: Vec<u8>,
    pub min_len: i64,
    pub max_len: i64,
    pub fwd_region: Region,
    pub rev_region: Region,
}

impl Window {
    /// The reverse primer as it appears in the rRNA-sense read.
    pub fn rev_site(&self) -> Vec<u8> {
        revcomp(&self.rev)
    }
}

fn default_regions(name: &str) -> (Region, Region) {
    let n = name.to_uppercase();
    if n.starts_with("SSU") {
        (Region::Ssu, Region::Ssu)
    } else if n.starts_with("LSU") {
        (Region::Lsu, Region::Lsu)
    } else if n.starts_with("AMPLICON") {
        (Region::Ssu, Region::Lsu)
    } else {
        (Region::Any, Region::Any)
    }
}

/// Python's `int(s)` for the integers a panel holds.
fn py_int(s: &str, what: &str) -> Result<i64> {
    s.trim()
        .parse::<i64>()
        .with_context(|| format!("panel: {what} is not an integer: {s:?}"))
}

pub fn load_panel(path: &Path) -> Result<Vec<Window>> {
    let mut text = Vec::new();
    open_maybe_gz(path)?.read_to_end(&mut text)?;
    let text =
        String::from_utf8(text).with_context(|| format!("{} is not UTF-8", path.display()))?;
    // rows = (line for line in fh if line.strip() and not line.startswith("#"))
    let mut lines = text
        .split_inclusive('\n')
        .filter(|l| !py_strip(l.as_bytes()).is_empty() && !l.starts_with('#'))
        .map(|l| l.trim_end_matches(['\n', '\r']));
    let header: Vec<&str> = match lines.next() {
        Some(h) => h.split('\t').collect(),
        None => return Ok(Vec::new()),
    };
    let col = |name: &str| header.iter().position(|h| *h == name);
    let need = |name: &str| {
        col(name).with_context(|| format!("panel {} has no column {name:?}", path.display()))
    };
    let (iw, it, ifn, ifs, irn, irs, imin, imax) = (
        need("window")?,
        need("target")?,
        need("fwd_name")?,
        need("fwd_seq")?,
        need("rev_name")?,
        need("rev_seq")?,
        need("min_len")?,
        need("max_len")?,
    );
    let (ifr, irr) = (col("fwd_region"), col("rev_region"));
    let mut windows = Vec::new();
    let mut seen = std::collections::HashSet::new();
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        let get = |i: usize| -> Result<&str> {
            f.get(i)
                .copied()
                .with_context(|| format!("panel row has too few fields: {line:?}"))
        };
        let name = get(iw)?.to_string();
        let (dfr, drr) = default_regions(&name);
        let region = |idx: Option<usize>, default: Region| -> Result<Region> {
            match idx.and_then(|i| f.get(i)).map(|s| s.to_lowercase()) {
                Some(s) if !s.is_empty() => Region::parse(&s).with_context(|| {
                    format!("panel window {name}: region must be one of ssu, lsu, any")
                }),
                _ => Ok(default),
            }
        };
        let w = Window {
            target: get(it)?.to_string(),
            fwd_name: get(ifn)?.to_string(),
            fwd: get(ifs)?.to_uppercase().into_bytes(),
            rev_name: get(irn)?.to_string(),
            rev: get(irs)?.to_uppercase().into_bytes(),
            min_len: py_int(get(imin)?, "min_len")?,
            max_len: py_int(get(imax)?, "max_len")?,
            fwd_region: region(ifr, dfr)?,
            rev_region: region(irr, drr)?,
            name,
        };
        if !seen.insert(w.name.clone()) {
            // the prototype opens one output file per name, so a repeated name would lose output
            bail!("panel {}: window {} appears twice", path.display(), w.name);
        }
        windows.push(w);
    }
    Ok(windows)
}

use std::io::Read as _;
