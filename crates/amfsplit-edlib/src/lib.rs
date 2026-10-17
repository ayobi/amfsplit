//! Minimal safe wrapper around a vendored copy of edlib (Martin Šošić, MIT licence).
//!
//! The C++ source is the one shipped with the Python `edlib` package 1.3.9.post1, so results
//! equal those of `edlib.align(...)` in Python, which the Stage 1 prototype used. Only what the
//! amfsplit port needs is exposed: one call returning the edit distance and the locations, in
//! the shape the Python binding returns them.

use std::os::raw::{c_char, c_int, c_uchar};

#[repr(C)]
#[derive(Clone, Copy)]
struct EdlibEqualityPair {
    first: c_char,
    second: c_char,
}

#[repr(C)]
struct EdlibAlignConfig {
    k: c_int,
    mode: c_int,
    task: c_int,
    additional_equalities: *const EdlibEqualityPair,
    additional_equalities_length: c_int,
}

#[repr(C)]
struct EdlibAlignResult {
    status: c_int,
    edit_distance: c_int,
    end_locations: *mut c_int,
    start_locations: *mut c_int,
    num_locations: c_int,
    alignment: *mut c_uchar,
    alignment_length: c_int,
    alphabet_length: c_int,
}

unsafe extern "C" {
    fn edlibAlign(
        query: *const c_char,
        query_length: c_int,
        target: *const c_char,
        target_length: c_int,
        config: EdlibAlignConfig,
    ) -> EdlibAlignResult;
    fn edlibFreeAlignResult(result: EdlibAlignResult);
}

/// Alignment method, as in edlib (`NW` global, `SHW` prefix, `HW` infix).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Nw,
    Shw,
    Hw,
}

/// What edlib computes. `Path` is not exposed: nothing in the port needs it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Task {
    /// Edit distance and end locations.
    Distance,
    /// Edit distance, end and start locations.
    Locations,
}

/// Result in the shape of the Python binding: `edit_distance` is -1 when larger than k, and
/// `locations` holds `(start, end)` pairs, 0-based inclusive, `start` absent when not computed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Alignment {
    pub edit_distance: i32,
    pub locations: Vec<(Option<i32>, i32)>,
    pub alphabet_length: i32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EdlibError;

impl std::fmt::Display for EdlibError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "edlib reported an error")
    }
}

impl std::error::Error for EdlibError {}

/// Align `query` to `target`. `k` < 0 means no limit. `equalities` are extra pairs of bytes
/// that count as equal (edlib treats each pair symmetrically; equality is not transitive).
pub fn align(
    query: &[u8],
    target: &[u8],
    mode: Mode,
    task: Task,
    k: i32,
    equalities: &[(u8, u8)],
) -> Result<Alignment, EdlibError> {
    let pairs: Vec<EdlibEqualityPair> = equalities
        .iter()
        .map(|&(a, b)| EdlibEqualityPair {
            first: a as c_char,
            second: b as c_char,
        })
        .collect();
    let config = EdlibAlignConfig {
        k,
        mode: match mode {
            Mode::Nw => 0,
            Mode::Shw => 1,
            Mode::Hw => 2,
        },
        task: match task {
            Task::Distance => 0,
            Task::Locations => 1,
        },
        additional_equalities: if pairs.is_empty() {
            std::ptr::null()
        } else {
            pairs.as_ptr()
        },
        additional_equalities_length: pairs.len() as c_int,
    };
    // SAFETY: the pointers are valid for the given lengths for the duration of the call, edlib
    // copies the sequences internally, and the result is freed exactly once below.
    let res = unsafe {
        edlibAlign(
            query.as_ptr() as *const c_char,
            query.len() as c_int,
            target.as_ptr() as *const c_char,
            target.len() as c_int,
            config,
        )
    };
    let out = if res.status != 0 {
        Err(EdlibError)
    } else {
        let mut locations = Vec::new();
        if res.num_locations > 0 && !res.end_locations.is_null() {
            let n = res.num_locations as usize;
            // SAFETY: edlib allocated num_locations ints for each non-null location array.
            let ends = unsafe { std::slice::from_raw_parts(res.end_locations, n) };
            let starts = if res.start_locations.is_null() {
                None
            } else {
                Some(unsafe { std::slice::from_raw_parts(res.start_locations, n) })
            };
            for i in 0..n {
                locations.push((starts.map(|s| s[i]), ends[i]));
            }
        }
        Ok(Alignment {
            edit_distance: res.edit_distance,
            locations,
            alphabet_length: res.alphabet_length,
        })
    };
    // SAFETY: res came from edlibAlign and is not used afterwards.
    unsafe { edlibFreeAlignResult(res) };
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn infix_locations() {
        let a = align(b"ACGT", b"TTACGTTT", Mode::Hw, Task::Locations, -1, &[]).unwrap();
        assert_eq!(a.edit_distance, 0);
        assert_eq!(a.locations, vec![(Some(2), 5)]);
    }

    #[test]
    fn k_limit_gives_minus_one() {
        let a = align(b"AAAAAAAA", b"CCCCCCCC", Mode::Hw, Task::Locations, 2, &[]).unwrap();
        assert_eq!(a.edit_distance, -1);
        assert!(a.locations.is_empty());
    }

    #[test]
    fn iupac_equalities() {
        let eq = [(b'R', b'A'), (b'R', b'G')];
        let a = align(b"ARGT", b"TTAGGTTT", Mode::Hw, Task::Locations, 0, &eq).unwrap();
        assert_eq!(a.edit_distance, 0);
        // R does not equal C
        let b = align(b"ARGT", b"TTACGTTT", Mode::Hw, Task::Locations, 0, &eq).unwrap();
        assert_eq!(b.edit_distance, -1);
    }

    #[test]
    fn global_distance() {
        let a = align(b"ACGTACGT", b"ACGTTCGT", Mode::Nw, Task::Distance, -1, &[]).unwrap();
        assert_eq!(a.edit_distance, 1);
    }

    #[test]
    fn empty_target_infix() {
        // edlib's special case: distance = query length even above k, end -1, no start.
        let a = align(b"ACGT", b"", Mode::Hw, Task::Locations, 1, &[]).unwrap();
        assert_eq!(a.edit_distance, 4);
        assert_eq!(a.locations, vec![(None, -1)]);
    }
}
