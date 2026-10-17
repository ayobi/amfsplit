//! `amfsplit triage --amf-reference`: a target is AMF if and only if its ID is a record of the AMF
//! reference. Not in the v0.3.1 prototype, so there is no golden file of its own. On V18 the
//! AMF-only file and the regex `_[A-Z][a-z]+aceae_` select the same 460 of 495 records (checked
//! when the module was written), so --amf-reference must reproduce the golden case
//! `triage_amf_pattern` call for call: triage.tsv, amf.ids, stdout and every AMF-only FASTA,
//! and the summary except its params.
//!
//! The self-contained test builds the AMF reference from the hit targets that the pattern
//! selects. With AMFSPLIT_TEST_V18_AMF set to V18_LSUDB_052025_AMFONLY.fasta, the same check
//! runs on the real file.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use flate2::read::GzDecoder;
use sha2::{Digest, Sha256};

const V18_AMF_SHA256: &str = "e6794cbbe166735a94e44b146aae2ecffb35ebfee61cf6ae8ed2e48ba9227198";

fn dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../parity/triage")
        .canonicalize()
        .unwrap()
}

fn sha256(path: &Path) -> String {
    format!("{:x}", Sha256::digest(fs::read(path).unwrap()))
}

fn workdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("amfsplit-amfref-{tag}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    for entry in fs::read_dir(dir().join("inputs")).unwrap() {
        let p = entry.unwrap().path();
        let name = p.file_name().unwrap().to_string_lossy().into_owned();
        if let Some(plain) = name.strip_suffix(".gz") {
            let mut s = Vec::new();
            GzDecoder::new(fs::File::open(&p).unwrap())
                .read_to_end(&mut s)
                .unwrap();
            fs::write(d.join(plain), s).unwrap();
        }
    }
    d
}

fn triage(d: &Path, extra: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_amfsplit"))
        .arg("triage")
        .arg("--hits")
        .arg(dir().join("inputs/hits.tsv"))
        .arg("--out")
        .arg(d.join("out"))
        .arg("--records")
        .arg(d.join("windows.records.tsv"))
        .arg("--fasta-prefix")
        .arg(d.join("windows"))
        .args(extra)
        .output()
        .unwrap()
}

/// Distinct hit targets, in first-seen order.
fn hit_targets() -> Vec<String> {
    let mut seen = Vec::new();
    for line in fs::read_to_string(dir().join("inputs/hits.tsv"))
        .unwrap()
        .lines()
    {
        let t = line.split('\t').nth(1).unwrap().to_string();
        if !seen.contains(&t) {
            seen.push(t);
        }
    }
    seen
}

fn is_family_named(t: &str) -> bool {
    // the regex _[A-Z][a-z]+aceae_, searched anywhere in the name
    t.split('_').collect::<Vec<_>>().windows(3).any(|w| {
        let f = w[1];
        f.len() > 5
            && f.ends_with("aceae")
            && f.chars().next().is_some_and(|c| c.is_ascii_uppercase())
            && f[1..f.len() - 5].chars().all(|c| c.is_ascii_lowercase())
    })
}

/// Everything but the summary must equal golden triage_amf_pattern; the summary must equal it
/// apart from params, which must name the AMF reference.
fn assert_matches_pattern_case(d: &Path, out: &Output, reference: &Path, n_records: i64) {
    assert!(
        out.status.success(),
        "triage failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    fs::write(d.join("stdout.txt"), &out.stdout).unwrap();
    let gold = dir().join("golden/triage_amf_pattern");
    let manifest = fs::read_to_string(gold.join("MANIFEST.tsv")).unwrap();
    let mut failures = Vec::new();
    for row in manifest.lines().skip(1) {
        let f: Vec<&str> = row.split('\t').collect();
        if f[0] == "triage.summary.json" {
            continue;
        }
        let mine = if f[0] == "stdout.txt" || f[0].ends_with(".amf.fasta") {
            d.join(f[0])
        } else {
            d.join(format!("out.{}", f[0]))
        };
        if sha256(&mine) != f[1] {
            failures.push(f[0].to_string());
        }
    }
    assert!(
        failures.is_empty(),
        "differ from triage_amf_pattern: {failures:?}"
    );

    let mut mine: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(d.join("out.triage.summary.json")).unwrap())
            .unwrap();
    let mut golden: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(gold.join("triage.summary.json")).unwrap())
            .unwrap();
    let params = mine.as_object_mut().unwrap().remove("params").unwrap();
    golden.as_object_mut().unwrap().remove("params");
    assert_eq!(mine, golden, "summary differs outside params");
    assert_eq!(params["amf_pattern"], serde_json::Value::Null);
    assert_eq!(params["outgroups"], serde_json::Value::Null);
    assert_eq!(params["amf_reference"]["records"], n_records);
    assert_eq!(params["amf_reference"]["sha256"], sha256(reference));
    assert_eq!(params["min_id"], 90.0);
}

#[test]
fn amf_reference_reproduces_the_pattern_case() {
    let d = workdir("derived");
    let amf: Vec<String> = hit_targets()
        .into_iter()
        .filter(|t| is_family_named(t))
        .collect();
    let reference = d.join("amf_reference.fasta");
    fs::write(
        &reference,
        amf.iter()
            .map(|t| format!(">{t} amf\nN\n"))
            .collect::<String>(),
    )
    .unwrap();
    let out = triage(&d, &["--amf-reference", reference.to_str().unwrap()]);
    assert_matches_pattern_case(&d, &out, &reference, amf.len() as i64);
    // the V18 hint is for --outgroups only
    assert!(
        !String::from_utf8_lossy(&out.stderr).contains("note:"),
        "unexpected note"
    );
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn amf_reference_v18_file() {
    let Ok(v18_amf) = std::env::var("AMFSPLIT_TEST_V18_AMF") else {
        eprintln!(
            "SKIPPED: set AMFSPLIT_TEST_V18_AMF to V18_LSUDB_052025_AMFONLY.fasta to run this test"
        );
        return;
    };
    let reference = PathBuf::from(&v18_amf);
    assert_eq!(
        sha256(&reference),
        V18_AMF_SHA256,
        "not V18_LSUDB_052025_AMFONLY.fasta: {v18_amf}"
    );
    let d = workdir("v18");
    let out = triage(&d, &["--amf-reference", &v18_amf]);
    assert_matches_pattern_case(&d, &out, &reference, 460);
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn no_amf_rule_is_an_error() {
    let d = workdir("norule");
    let out = triage(&d, &[]);
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(!out.status.success());
    assert!(
        err.contains("say which hit targets are AMF") && err.contains("--amf-reference"),
        "{err}"
    );
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn reference_from_another_release_is_an_error() {
    let d = workdir("mismatch");
    let reference = d.join("other.fasta");
    fs::write(&reference, ">XX000001.1_Glomeraceae_Glomus_sp\nACGT\n").unwrap();
    let out = triage(&d, &["--amf-reference", reference.to_str().unwrap()]);
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(!out.status.success());
    assert!(err.contains("same database release"), "{err}");
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn every_target_amf_gives_a_note() {
    let d = workdir("allamf");
    let reference = d.join("all.fasta");
    fs::write(
        &reference,
        hit_targets()
            .iter()
            .map(|t| format!(">{t}\nN\n"))
            .collect::<String>(),
    )
    .unwrap();
    let out = triage(&d, &["--amf-reference", reference.to_str().unwrap()]);
    assert!(out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("no read can be called non_AMF"), "{err}");
    let _ = fs::remove_dir_all(&d);
}

#[test]
fn amf_reference_conflicts_with_other_rules() {
    let d = workdir("conflict");
    for other in [["--amf-pattern", "x"], ["--outgroups", "Mortierella"]] {
        let out = triage(&d, &["--amf-reference", "/dev/null", other[0], other[1]]);
        assert!(
            !out.status.success(),
            "{other:?} accepted with --amf-reference"
        );
    }
    let _ = fs::remove_dir_all(&d);
}
