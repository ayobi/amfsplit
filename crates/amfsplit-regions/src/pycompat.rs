//! Small pieces of Python semantics that the Stage 1 outputs depend on byte for byte.
//!
//! The prototype (`amfsplit.py` v0.3.1) writes integers, floats, slices, medians, quartiles,
//! rounded percentages, CSV fields and JSON exactly as CPython does; this module reproduces the
//! parts it uses, and nothing more.

/// Python's `str.isspace()` restricted to ASCII: the characters `str.split()` and `str.strip()`
/// treat as whitespace in ASCII text.
pub fn is_py_space(b: u8) -> bool {
    matches!(
        b,
        b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c | 0x1c | 0x1d | 0x1e | 0x1f
    )
}

/// `s.strip()`.
pub fn py_strip(s: &[u8]) -> &[u8] {
    let start = s.iter().position(|&b| !is_py_space(b)).unwrap_or(s.len());
    let end = s
        .iter()
        .rposition(|&b| !is_py_space(b))
        .map_or(start, |i| i + 1);
    &s[start..end]
}

/// `s.split()[0]`, or None when there is no token (Python raises IndexError there).
pub fn py_first_token(s: &[u8]) -> Option<&[u8]> {
    let start = s.iter().position(|&b| !is_py_space(b))?;
    let len = s[start..]
        .iter()
        .position(|&b| is_py_space(b))
        .unwrap_or(s.len() - start);
    Some(&s[start..start + len])
}

/// Resolve one slice bound the way CPython does for a sequence of length `n`.
fn bound(i: i64, n: i64) -> i64 {
    if i < 0 { (i + n).max(0) } else { i.min(n) }
}

/// `s[start:end]` with Python's handling of negative and out-of-range bounds.
pub fn py_slice(s: &[u8], start: i64, end: i64) -> &[u8] {
    let n = s.len() as i64;
    let (a, b) = (bound(start, n), bound(end, n));
    if a >= b {
        &s[0..0]
    } else {
        &s[a as usize..b as usize]
    }
}

/// A Python number as the prototype produces it: `int` or `float`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Num {
    Int(i64),
    Float(f64),
}

impl Num {
    pub fn as_f64(self) -> f64 {
        match self {
            Num::Int(i) => i as f64,
            Num::Float(f) => f,
        }
    }
}

/// `statistics.median` of integers: the middle value (int) for odd n, the mean of the two
/// middle values (float) for even n. Panics on empty input, as Python raises.
pub fn median(values: &[i64]) -> Num {
    let mut v = values.to_vec();
    v.sort_unstable();
    let n = v.len();
    assert!(n > 0, "median of no data");
    if n % 2 == 1 {
        Num::Int(v[n / 2])
    } else {
        // (a + b) / 2 is Python true division: exact operands, one correctly rounded division
        Num::Float((v[n / 2 - 1] + v[n / 2]) as f64 / 2.0)
    }
}

/// `int(statistics.median(values))`: truncation toward zero.
pub fn median_trunc(values: &[i64]) -> i64 {
    match median(values) {
        Num::Int(i) => i,
        Num::Float(f) => f.trunc() as i64,
    }
}

/// `statistics.quantiles(values, n=n)` with the default 'exclusive' method, for integer data
/// with at least two values.
pub fn quantiles(values: &[i64], n: i64) -> Vec<f64> {
    let mut data = values.to_vec();
    data.sort_unstable();
    let ld = data.len() as i64;
    assert!(ld >= 2, "quantiles need at least two values here");
    let m = ld + 1;
    (1..n)
        .map(|i| {
            let mut j = i * m / n; // floor division, all operands positive
            j = j.clamp(1, ld - 1);
            let delta = i * m - j * n;
            let num = data[(j - 1) as usize] * (n - delta) + data[j as usize] * delta;
            num as f64 / n as f64
        })
        .collect()
}

/// `statistics.median` of floats: the middle value for odd n, the mean of the two middle
/// values for even n. Panics on empty input, as Python raises.
pub fn median_f64(values: &[f64]) -> f64 {
    let mut v = values.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    let n = v.len();
    assert!(n > 0, "median of no data");
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

/// `Counter(items).most_common(n)`: counts in first-seen order, then a stable sort by count,
/// highest first (ties keep first-seen order, as `sorted(..., reverse=True)` and
/// `heapq.nlargest` do).
pub fn most_common<'a, I: IntoIterator<Item = &'a str>>(
    items: I,
    n: Option<usize>,
) -> Vec<(String, i64)> {
    let mut c = OrderedCounter::default();
    for it in items {
        c.add(it, 1);
    }
    let mut v: Vec<(String, i64)> = c.items().to_vec();
    v.sort_by_key(|x| std::cmp::Reverse(x.1)); // stable: ties keep first-seen order
    if let Some(n) = n {
        v.truncate(n);
    }
    v
}

/// Rows of a tab-separated file as Python's `csv.reader(fh, delimiter="\t")` parses them: a field
/// that starts with a quote is quoted (doubled quotes inside it are one quote, and it may span
/// lines); blank rows are dropped, as `csv.DictReader` drops them.
pub fn read_tsv_rows(text: &str) -> Vec<Vec<String>> {
    let mut rows = Vec::new();
    let mut row: Vec<String> = Vec::new();
    let mut field = String::new();
    let mut chars = text.chars().peekable();
    let (mut in_quotes, mut at_field_start, mut row_has_content) = (false, true, false);
    while let Some(c) = chars.next() {
        if in_quotes {
            if c == '"' {
                if chars.peek() == Some(&'"') {
                    chars.next();
                    field.push('"');
                } else {
                    in_quotes = false;
                }
            } else {
                field.push(c);
            }
            continue;
        }
        match c {
            '"' if at_field_start => {
                in_quotes = true;
                at_field_start = false;
                row_has_content = true;
            }
            '\t' => {
                row.push(std::mem::take(&mut field));
                at_field_start = true;
                row_has_content = true;
            }
            '\r' if chars.peek() == Some(&'\n') => {}
            '\n' | '\r' => {
                if row_has_content || !field.is_empty() {
                    row.push(std::mem::take(&mut field));
                    rows.push(std::mem::take(&mut row));
                }
                at_field_start = true;
                row_has_content = false;
            }
            c => {
                field.push(c);
                at_field_start = false;
                row_has_content = true;
            }
        }
    }
    if row_has_content || !field.is_empty() {
        row.push(field);
        rows.push(row);
    }
    rows
}

/// `round(x, ndigits)` for a float: Python and Rust both round the exact binary value, ties to
/// even, so formatting with the precision and parsing back gives Python's result.
pub fn py_round(x: f64, ndigits: usize) -> f64 {
    format!("{:.*}", ndigits, x)
        .parse()
        .expect("formatted float parses")
}

/// `repr(float)`: shortest round-trip digits, fixed notation for 1e-4 <= |x| < 1e16 (with a
/// trailing ".0" when integral), otherwise scientific with a signed two-digit exponent.
pub fn float_repr(x: f64) -> String {
    if x.is_nan() {
        return "nan".into();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    let a = x.abs();
    if x == 0.0 || (1e-4..1e16).contains(&a) {
        let s = format!("{x:?}"); // Rust's Debug is shortest round-trip with ".0" when integral
        if s.contains('e') {
            // not expected in this range; fall through to the general form below
        } else {
            return s;
        }
    }
    let s = format!("{x:e}"); // e.g. 1e-5, 1.5e16
    let (mant, exp) = s.split_once('e').expect("scientific form");
    let exp: i32 = exp.parse().expect("exponent");
    format!("{mant}e{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs())
}

/// One field as `csv.writer(..., delimiter='\t')` writes it with QUOTE_MINIMAL: quoted (inner
/// quotes doubled) when it holds a tab, a quote, CR or LF.
pub fn csv_field(s: &str) -> std::borrow::Cow<'_, str> {
    if s.bytes().any(|b| matches!(b, b'\t' | b'"' | b'\n' | b'\r')) {
        format!("\"{}\"", s.replace('"', "\"\"")).into()
    } else {
        s.into()
    }
}

/// One row of a tab-separated `csv.writer` with LF line endings. A row of a single empty field
/// is written as `""`, as Python does.
pub fn csv_row(fields: &[String]) -> String {
    if fields.len() == 1 && fields[0].is_empty() {
        return "\"\"\n".into();
    }
    let mut out = String::new();
    for (i, f) in fields.iter().enumerate() {
        if i > 0 {
            out.push('\t');
        }
        out.push_str(&csv_field(f));
    }
    out.push('\n');
    out
}

/// A JSON value with the key order kept, written as `json.dump(value, fh, indent=1)` writes it.
#[derive(Clone, Debug, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Int(i64),
    Float(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl From<Num> for Json {
    fn from(n: Num) -> Self {
        match n {
            Num::Int(i) => Json::Int(i),
            Num::Float(f) => Json::Float(f),
        }
    }
}

impl Json {
    pub fn obj<K: Into<String>>(pairs: Vec<(K, Json)>) -> Json {
        Json::Obj(pairs.into_iter().map(|(k, v)| (k.into(), v)).collect())
    }

    /// `json.dumps(self, indent=1)` (ensure_ascii=True, no trailing newline).
    pub fn dump_indent1(&self) -> String {
        let mut out = String::new();
        self.write(&mut out, 0);
        out
    }

    fn write(&self, out: &mut String, level: usize) {
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            Json::Int(i) => out.push_str(&i.to_string()),
            Json::Float(f) => out.push_str(&json_float(*f)),
            Json::Str(s) => write_json_str(out, s),
            Json::Arr(items) => {
                if items.is_empty() {
                    out.push_str("[]");
                    return;
                }
                out.push('[');
                for (i, v) in items.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    out.push('\n');
                    out.push_str(&" ".repeat(level + 1));
                    v.write(out, level + 1);
                }
                out.push('\n');
                out.push_str(&" ".repeat(level));
                out.push(']');
            }
            Json::Obj(pairs) => {
                if pairs.is_empty() {
                    out.push_str("{}");
                    return;
                }
                out.push('{');
                for (i, (k, v)) in pairs.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    out.push('\n');
                    out.push_str(&" ".repeat(level + 1));
                    write_json_str(out, k);
                    out.push_str(": ");
                    v.write(out, level + 1);
                }
                out.push('\n');
                out.push_str(&" ".repeat(level));
                out.push('}');
            }
        }
    }
}

/// Python's json float output: `float.__repr__`, with NaN and infinities as JSON extensions.
fn json_float(f: f64) -> String {
    if f.is_nan() {
        "NaN".into()
    } else if f.is_infinite() {
        if f > 0.0 {
            "Infinity".into()
        } else {
            "-Infinity".into()
        }
    } else {
        float_repr(f)
    }
}

fn write_json_str(out: &mut String, s: &str) {
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            c if (c as u32) < 0x20 || (c as u32) > 0x7e => {
                let mut buf = [0u16; 2];
                for unit in c.encode_utf16(&mut buf) {
                    out.push_str(&format!("\\u{:04x}", unit));
                }
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

/// A dict with insertion order, as the prototype's `Counter`s and dicts keep it.
#[derive(Clone, Debug, Default)]
pub struct OrderedCounter {
    items: Vec<(String, i64)>,
}

impl OrderedCounter {
    pub fn add(&mut self, key: &str, by: i64) {
        match self.items.iter_mut().find(|(k, _)| k == key) {
            Some((_, v)) => *v += by,
            None => self.items.push((key.to_string(), by)),
        }
    }

    /// `counter[key]`: 0 when absent, without inserting it.
    pub fn get(&self, key: &str) -> i64 {
        self.items
            .iter()
            .find(|(k, _)| k == key)
            .map_or(0, |(_, v)| *v)
    }

    pub fn items(&self) -> &[(String, i64)] {
        &self.items
    }

    pub fn to_json(&self) -> Json {
        Json::Obj(
            self.items
                .iter()
                .map(|(k, v)| (k.clone(), Json::Int(*v)))
                .collect(),
        )
    }

    /// `str(dict(counter))`, for the console report.
    pub fn py_dict_repr(&self) -> String {
        let inner: Vec<String> = self
            .items
            .iter()
            .map(|(k, v)| format!("'{k}': {v}"))
            .collect();
        format!("{{{}}}", inner.join(", "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slices_follow_python() {
        let s = b"ABCDEFG";
        assert_eq!(py_slice(s, 2, 4), b"CD");
        assert_eq!(py_slice(s, -3, 100), b"EFG");
        assert_eq!(py_slice(s, 5, 2), b"");
        assert_eq!(py_slice(s, 0, -1), b"ABCDEF");
        assert_eq!(py_slice(s, -100, 2), b"AB");
    }

    #[test]
    fn median_and_quantiles_follow_python() {
        assert_eq!(median(&[3, 1, 2]), Num::Int(2));
        assert_eq!(median(&[4, 1, 2, 3]), Num::Float(2.5));
        assert_eq!(median_trunc(&[1, 2]), 1);
        // statistics.quantiles([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], n=4) == [2.75, 5.5, 8.25]
        assert_eq!(
            quantiles(&[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 4),
            vec![2.75, 5.5, 8.25]
        );
        // statistics.quantiles([5, 1], n=4) == [0.0, 3.0, 6.0] (the exclusive method extrapolates)
        assert_eq!(quantiles(&[5, 1], 4), vec![0.0, 3.0, 6.0]);
    }

    #[test]
    fn rounding_and_repr_follow_python() {
        assert_eq!(py_round(0.25, 1), 0.2);
        assert_eq!(py_round(2.675, 2), 2.67);
        assert_eq!(float_repr(466.0), "466.0");
        assert_eq!(float_repr(0.15), "0.15");
        assert_eq!(float_repr(1e-5), "1e-05");
        assert_eq!(float_repr(1e16), "1e+16");
        assert_eq!(float_repr(96.9), "96.9");
        assert_eq!(float_repr(0.0), "0.0");
    }

    #[test]
    fn csv_quoting_follows_python() {
        // csv.writer(delimiter='\t', lineterminator='\n').writerow(['a\rb', 'c"d', 'e f', '', '3', 'x\ty'])
        let row: Vec<String> = ["a\rb", "c\"d", "e f", "", "3", "x\ty"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert_eq!(csv_row(&row), "\"a\rb\"\t\"c\"\"d\"\te f\t\t3\t\"x\ty\"\n");
        assert_eq!(csv_row(&[String::new()]), "\"\"\n");
    }

    #[test]
    fn json_follows_python_indent1() {
        let j = Json::obj(vec![
            ("a", Json::Str("é\u{1}\"\\".into())),
            ("c", Json::Obj(vec![])),
            (
                "d",
                Json::obj(vec![("x", Json::Float(1.0)), ("y", Json::Null)]),
            ),
        ]);
        assert_eq!(
            j.dump_indent1(),
            "{\n \"a\": \"\\u00e9\\u0001\\\"\\\\\",\n \"c\": {},\n \"d\": {\n  \"x\": 1.0,\n  \"y\": null\n }\n}"
        );
    }

    #[test]
    fn arrays_and_counters_follow_python() {
        let j = Json::obj(vec![
            (
                "t",
                Json::Arr(vec![Json::Arr(vec![Json::Str("a".into()), Json::Int(2)])]),
            ),
            ("e", Json::Arr(vec![])),
        ]);
        assert_eq!(
            j.dump_indent1(),
            "{\n \"t\": [\n  [\n   \"a\",\n   2\n  ]\n ],\n \"e\": []\n}"
        );
        // Counter("abbcc").most_common(2) == [('b', 2), ('c', 2)]
        assert_eq!(
            most_common(["a", "b", "b", "c", "c"], Some(2)),
            vec![("b".to_string(), 2), ("c".to_string(), 2)]
        );
        assert_eq!(median_f64(&[3.0, 1.0, 2.0, 10.0]), 2.5);
    }

    #[test]
    fn tsv_rows_follow_python_csv() {
        let rows = read_tsv_rows("a\tb\n\n\"x\ty\"\t\"q\"\"q\"\r\nlast");
        assert_eq!(
            rows,
            vec![
                vec!["a".to_string(), "b".to_string()],
                vec!["x\ty".to_string(), "q\"q".to_string()],
                vec!["last".to_string()],
            ]
        );
    }

    #[test]
    fn strip_and_tokens() {
        assert_eq!(py_strip(b"  ACGT\r\n"), b"ACGT");
        assert_eq!(py_first_token(b"read1 desc\n"), Some(&b"read1"[..]));
        assert_eq!(py_first_token(b"   \n"), None);
    }
}
