//! Triage parity: `amfsplit triage` and `amfsplit validate` must reproduce triage.py and
//! validate_windows.py (amfsplit v0.3.1) byte for byte. Golden files: parity/triage/make_golden.py.
//!
//! The triage cases run on stored inputs (Stage 1 windows plus test windows, and the vsearch hits
//! of validate_windows.py). The validate case needs vsearch 2.27.0 and the Delavaux V18 database
//! with outgroups, which are not stored here: set AMFSPLIT_TEST_VSEARCH and AMFSPLIT_TEST_V18 to
//! run it; without them it is skipped, and says so.

use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

use flate2::read::GzDecoder;
use sha2::{Digest, Sha256};

fn dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../parity/triage")
        .canonicalize()
        .unwrap()
}

fn sha256(path: &Path) -> String {
    let buf = fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    format!("{:x}", Sha256::digest(&buf))
}

fn tsv(path: &Path) -> Vec<Vec<String>> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect()
}

fn unpack_inputs(to: &Path) {
    for entry in fs::read_dir(dir().join("inputs")).unwrap() {
        let p = entry.unwrap().path();
        let name = p.file_name().unwrap().to_string_lossy().into_owned();
        if let Some(plain) = name.strip_suffix(".gz") {
            let mut s = Vec::new();
            GzDecoder::new(fs::File::open(&p).unwrap())
                .read_to_end(&mut s)
                .unwrap();
            fs::write(to.join(plain), s).unwrap();
        }
    }
}

/// Compare every file in golden/<case>/MANIFEST.tsv with `files` (label -> path).
fn check(case: &str, files: &dyn Fn(&str) -> PathBuf, failures: &mut Vec<String>) {
    let manifest = tsv(&dir().join("golden").join(case).join("MANIFEST.tsv"));
    for row in &manifest[1..] {
        let path = files(&row[0]);
        if !path.exists() {
            failures.push(format!("{case}: {} missing", row[0]));
        } else if sha256(&path) != row[1] {
            failures.push(format!("{case}: {} differs", row[0]));
        }
    }
}

#[test]
fn triage_outputs_match_the_prototype() {
    let cases = tsv(&dir().join("cases.tsv"));
    let work = std::env::temp_dir().join(format!("amfsplit-triage-test-{}", std::process::id()));
    let mut failures = Vec::new();
    for c in &cases[1..] {
        let (case, opts) = (&c[0], &c[1]);
        let d = work.join(case);
        fs::create_dir_all(&d).unwrap();
        unpack_inputs(&d);
        let opts = opts
            .replace("{hits}", &dir().join("inputs/hits.tsv").to_string_lossy())
            .replace("{out}", &d.join("out").to_string_lossy())
            .replace(
                "{records}",
                &d.join("windows.records.tsv").to_string_lossy(),
            )
            .replace("{prefix}", &d.join("windows").to_string_lossy());
        let res = Command::new(env!("CARGO_BIN_EXE_amfsplit"))
            .arg("triage")
            .args(opts.split_whitespace())
            .output()
            .unwrap();
        if !res.status.success() {
            failures.push(format!(
                "{case}: exit {:?}: {}",
                res.status.code(),
                String::from_utf8_lossy(&res.stderr)
            ));
            continue;
        }
        fs::write(d.join("stdout.txt"), &res.stdout).unwrap();
        let files = |label: &str| -> PathBuf {
            if label == "stdout.txt" || label.ends_with(".amf.fasta") {
                d.join(label)
            } else {
                d.join(format!("out.{label}"))
            }
        };
        check(case, &files, &mut failures);
    }
    let _ = fs::remove_dir_all(&work);
    assert!(
        failures.is_empty(),
        "{} mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}

#[test]
fn validate_outputs_match_the_prototype() {
    let (Ok(vsearch), Ok(v18)) = (
        std::env::var("AMFSPLIT_TEST_VSEARCH"),
        std::env::var("AMFSPLIT_TEST_V18"),
    ) else {
        eprintln!(
            "SKIPPED validate parity: set AMFSPLIT_TEST_VSEARCH (vsearch 2.27.0) and AMFSPLIT_TEST_V18 (V18_LSUDB_052025.fasta)"
        );
        return;
    };
    let d = std::env::temp_dir().join(format!("amfsplit-validate-test-{}", std::process::id()));
    fs::create_dir_all(&d).unwrap();
    unpack_inputs(&d);
    // restore the V18 sequences of the OUT_* test windows (stored masked)
    let mut v18_seqs = std::collections::HashMap::new();
    let mut cur = String::new();
    for line in fs::read_to_string(&v18).unwrap().lines() {
        if let Some(h) = line.strip_prefix('>') {
            cur = h.split_whitespace().next().unwrap_or("").to_string();
            v18_seqs.insert(cur.clone(), String::new());
        } else if let Some(s) = v18_seqs.get_mut(&cur) {
            s.push_str(line.trim());
        }
    }
    let masked = fs::read_to_string(d.join("windows.LSU_LROR-FLR2.fasta")).unwrap();
    let mut out = fs::File::create(d.join("lsu.fasta")).unwrap();
    let mut restore: Option<String> = None;
    for line in masked.lines() {
        if let Some(id) = line.strip_prefix(">OUT_") {
            restore = Some(v18_seqs[id.split_whitespace().next().unwrap()].to_uppercase());
            writeln!(out, "{line}").unwrap();
        } else if line.starts_with('>') {
            restore = None;
            writeln!(out, "{line}").unwrap();
        } else {
            writeln!(
                out,
                "{}",
                restore.take().unwrap_or_else(|| line.to_string())
            )
            .unwrap();
        }
    }
    drop(out);
    let res = Command::new(env!("CARGO_BIN_EXE_amfsplit"))
        .args(["validate", "-w"])
        .arg(d.join("lsu.fasta"))
        .arg("-r")
        .arg(&v18)
        .arg("-o")
        .arg(d.join("val"))
        .args(["--id", "0.70", "--threads", "1", "--vsearch", &vsearch])
        .output()
        .unwrap();
    assert!(
        res.status.success(),
        "validate failed: {}",
        String::from_utf8_lossy(&res.stderr)
    );
    fs::write(d.join("stdout.txt"), &res.stdout).unwrap();
    let mut failures = Vec::new();
    let files = |label: &str| -> PathBuf {
        if label == "stdout.txt" {
            d.join(label)
        } else {
            d.join(format!("val.{label}"))
        }
    };
    check("validate", &files, &mut failures);
    let _ = fs::remove_dir_all(&d);
    assert!(
        failures.is_empty(),
        "{} mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}
