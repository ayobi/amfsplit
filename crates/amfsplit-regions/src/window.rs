//! Cutting one legacy window from an oriented read (`cut_window` of the prototype).

use anyhow::Result;

use crate::panel::Window;
use crate::search::{find_site_in, max_edits};

/// Window statuses that yield a usable window.
pub const USABLE: [&str; 3] = ["ok", "open_5p", "open_3p"];

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WindowResult {
    pub status: &'static str,
    /// 1-based inclusive coordinates and length, when a window was cut.
    pub span: Option<(i64, i64, i64)>,
    pub fwd_ed: Option<i64>,
    pub rev_ed: Option<i64>,
    /// None only for reads that could not be oriented (no search at all).
    pub multi: Option<i64>,
}

impl WindowResult {
    pub fn unoriented() -> Self {
        WindowResult {
            status: "unoriented",
            span: None,
            fwd_ed: None,
            rev_ed: None,
            multi: None,
        }
    }

    pub fn usable(&self) -> bool {
        USABLE.contains(&self.status)
    }
}

/// Cut window `w` from `seq`. Bounds are 0-based half-open search regions per primer (None: the
/// whole read).
pub fn cut_window(
    w: &Window,
    seq: &[u8],
    frac: f64,
    keep_primers: bool,
    open_ends: bool,
    fwd_bounds: Option<(i64, i64)>,
    rev_bounds: Option<(i64, i64)>,
) -> Result<WindowResult> {
    let n = seq.len() as i64;
    let fb = fwd_bounds.unwrap_or((0, n));
    let rb = rev_bounds.unwrap_or((0, n));
    let rev_site = w.rev_site();
    let f = find_site_in(&w.fwd, seq, max_edits(w.fwd.len(), frac), fb.0, fb.1)?;
    let r = find_site_in(&rev_site, seq, max_edits(rev_site.len(), frac), rb.0, rb.1)?;
    let multi = i64::from(f.is_some_and(|h| h.n_best > 1) || r.is_some_and(|h| h.n_best > 1));
    let base = WindowResult {
        status: "",
        span: None,
        fwd_ed: f.map(|h| h.edits),
        rev_ed: r.map(|h| h.edits),
        multi: Some(multi),
    };
    let with = |status: &'static str| WindowResult {
        status,
        ..base.clone()
    };
    // starts are always present here: find_site_in has already rejected empty search regions
    let start = |h: crate::search::Hit| h.start.expect("start present after find_site_in");
    let (s0, e0, open_tag): (i64, i64, Option<&'static str>) = match (f, r) {
        (None, None) => return Ok(with("both_missing")),
        (None, Some(r)) => {
            // fwd site not found: open only if the read cuts its search region and the window
            // would still be no longer than max_len; otherwise the read covers it: missing
            let (s0, e0) = (0, if keep_primers { r.end } else { start(r) - 1 });
            let beyond_read = fb.0 == 0 && (e0 - s0 + 1) <= w.max_len;
            if !(open_ends && beyond_read) {
                return Ok(with("missing_5p"));
            }
            (s0, e0, Some("open_5p"))
        }
        (Some(f), None) => {
            let (s0, e0) = (if keep_primers { start(f) } else { f.end + 1 }, n - 1);
            let beyond_read = rb.1 == n && (e0 - s0 + 1) <= w.max_len;
            if !(open_ends && beyond_read) {
                return Ok(with("missing_3p"));
            }
            (s0, e0, Some("open_3p"))
        }
        (Some(f), Some(r)) => {
            if keep_primers {
                (start(f), r.end, None)
            } else {
                (f.end + 1, start(r) - 1, None)
            }
        }
    };
    if e0 < s0 {
        return Ok(with("out_of_order"));
    }
    let length = e0 - s0 + 1;
    let min_len = if open_tag.is_some() {
        w.min_len.div_euclid(2)
    } else {
        w.min_len
    };
    let status = if length < min_len {
        "too_short"
    } else if length > w.max_len {
        "too_long"
    } else {
        open_tag.unwrap_or("ok")
    };
    Ok(WindowResult {
        span: Some((s0 + 1, e0 + 1, length)),
        ..with(status)
    })
}
