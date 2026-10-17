//! IUPAC-aware primer-site search (edlib infix alignment), as in the prototype.

use anyhow::{Result, bail};

use amfsplit_edlib::{Mode, Task, align};

use crate::panel::Region;
use crate::pycompat::py_slice;

/// IUPAC codes and the bases each one equals, in the prototype's order (`EQUALITIES`).
const IUPAC: [(u8, &[u8]); 11] = [
    (b'R', b"AG"),
    (b'Y', b"CT"),
    (b'S', b"CG"),
    (b'W', b"AT"),
    (b'K', b"GT"),
    (b'M', b"AC"),
    (b'B', b"CGT"),
    (b'D', b"AGT"),
    (b'H', b"ACT"),
    (b'V', b"ACG"),
    (b'N', b"ACGT"),
];

fn equalities() -> &'static [(u8, u8)] {
    use std::sync::OnceLock;
    static EQ: OnceLock<Vec<(u8, u8)>> = OnceLock::new();
    EQ.get_or_init(|| {
        IUPAC
            .iter()
            .flat_map(|(code, bases)| bases.iter().map(move |&b| (*code, b)))
            .collect()
    })
}

/// Best infix match of a primer site: 0-based inclusive start and end, edits, and the number of
/// equally good end positions. `start` is None only in edlib's empty-target special case.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Hit {
    pub start: Option<i64>,
    pub end: i64,
    pub edits: i64,
    pub n_best: i64,
}

impl Hit {
    /// The start position; the prototype would fail (TypeError) where it is absent.
    pub fn start_or_err(&self) -> Result<i64> {
        match self.start {
            Some(s) => Ok(s),
            None => {
                bail!("primer search on an empty sequence (the v0.3.1 prototype fails here too)")
            }
        }
    }
}

/// `find_site(site, seq, max_err)`.
pub fn find_site(site: &[u8], seq: &[u8], max_err: i64) -> Result<Option<Hit>> {
    let k = i32::try_from(max_err).unwrap_or(i32::MAX);
    let r = align(site, seq, Mode::Hw, Task::Locations, k, equalities())?;
    if r.edit_distance < 0 {
        return Ok(None);
    }
    let Some(&(start, end)) = r.locations.first() else {
        bail!("edlib returned a distance but no location");
    };
    Ok(Some(Hit {
        start: start.map(i64::from),
        end: i64::from(end),
        edits: i64::from(r.edit_distance),
        n_best: r.locations.len() as i64,
    }))
}

/// `max_edits(primer, frac)`: at least 1, otherwise the integer part of len * frac.
pub fn max_edits(site_len: usize, frac: f64) -> i64 {
    ((site_len as f64 * frac).trunc() as i64).max(1)
}

/// `find_site_in(site, seq, max_err, lo, hi)`: the search restricted to `seq[lo:hi]`,
/// coordinates on the whole read.
pub fn find_site_in(
    site: &[u8],
    seq: &[u8],
    max_err: i64,
    lo: i64,
    hi: i64,
) -> Result<Option<Hit>> {
    if hi - lo < site.len() as i64 - max_err {
        return Ok(None);
    }
    let Some(hit) = find_site(site, py_slice(seq, lo, hi), max_err)? else {
        return Ok(None);
    };
    let start = hit.start_or_err()?;
    Ok(Some(Hit {
        start: Some(start + lo),
        end: hit.end + lo,
        ..hit
    }))
}

/// `region_bounds(region, full, read_len, margin)`: 0-based half-open search bounds.
pub fn region_bounds(
    region: Region,
    full: Option<(i64, i64)>,
    read_len: i64,
    margin: i64,
) -> (i64, i64) {
    match (full, region) {
        (None, _) | (_, Region::Any) => (0, read_len),
        (Some(f), Region::Ssu) => (0, read_len.min(f.0 - 1 + margin)),
        (Some(f), Region::Lsu) => ((f.1 - margin).max(0), read_len),
    }
}

/// `find_all_sites(site, seq, max_err, min_sep, max_sites)`: the best match, then the best
/// match again with the first one and `min_sep` around it masked by 'X', and so on.
pub fn find_all_sites(
    site: &[u8],
    seq: &[u8],
    max_err: i64,
    min_sep: i64,
    max_sites: usize,
) -> Result<Vec<(i64, i64, i64)>> {
    let mut masked = seq.to_vec();
    let mut sites = Vec::new();
    for _ in 0..max_sites {
        let Some(hit) = find_site(site, &masked, max_err)? else {
            break;
        };
        let start = hit.start_or_err()?;
        sites.push((start, hit.end, hit.edits));
        let n = masked.len() as i64;
        let lo = (start - min_sep).max(0);
        let hi = n.min(hit.end + min_sep + 1);
        let mut next = py_slice(&masked, 0, lo).to_vec();
        next.extend(std::iter::repeat_n(b'X', (hi - lo).max(0) as usize));
        next.extend_from_slice(py_slice(&masked, hi, n));
        masked = next;
    }
    Ok(sites)
}

/// `_norm_ed(a, b)`: global edit distance per base of the longer sequence; 1.0 if either is empty.
pub fn norm_ed(a: &[u8], b: &[u8]) -> Result<f64> {
    if a.is_empty() || b.is_empty() {
        return Ok(1.0);
    }
    let r = align(a, b, Mode::Nw, Task::Distance, -1, &[])?;
    Ok(f64::from(r.edit_distance) / a.len().max(b.len()) as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iupac_primer_matches_any_listed_base() {
        let hit = find_site(b"ACGRT", b"TTTACGGTTT", 0).unwrap().unwrap();
        assert_eq!((hit.start, hit.end, hit.edits), (Some(3), 7, 0));
    }

    #[test]
    fn read_n_matches_primer_bases_but_not_codes() {
        assert!(find_site(b"ACGT", b"ACNT", 0).unwrap().is_some());
        assert!(find_site(b"ACRT", b"ACNT", 0).unwrap().is_none());
    }

    #[test]
    fn restricted_search_reports_whole_read_coordinates() {
        let seq = b"ACGTTTTTTTTTACGT";
        let hit = find_site_in(b"ACGT", seq, 0, 5, 16).unwrap().unwrap();
        assert_eq!((hit.start, hit.end), (Some(12), 15));
        assert!(find_site_in(b"ACGT", seq, 0, 5, 7).unwrap().is_none());
    }

    #[test]
    fn edit_budget_is_at_least_one() {
        assert_eq!(max_edits(5, 0.15), 1);
        assert_eq!(max_edits(20, 0.15), 3);
        assert_eq!(max_edits(22, 0.20), 4);
    }
}
