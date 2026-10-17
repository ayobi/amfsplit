//! Stage 1 parity: `amfsplit regions` must reproduce the outputs of the prototype
//! (amfsplit.py v0.3.1) byte for byte on every case in parity/stage1/cases.tsv.
//!
//! The golden files were written by parity/stage1/make_golden.py; this test needs no Python.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

use sha2::{Digest, Sha256};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .unwrap()
}

fn sha256(path: &Path) -> String {
    let mut f = fs::File::open(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut buf = Vec::new();
    f.read_to_end(&mut buf).unwrap();
    format!("{:x}", Sha256::digest(&buf))
}

fn tsv(path: &Path) -> Vec<Vec<String>> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect()
}

#[test]
fn stage1_outputs_match_the_prototype() {
    let dir = root().join("parity/stage1");
    let cases = tsv(&dir.join("cases.tsv"));
    let header = &cases[0];
    let col = |name: &str| header.iter().position(|h| h == name).unwrap();
    let work = std::env::temp_dir().join(format!("amfsplit-parity-test-{}", std::process::id()));
    let mut failures = Vec::new();
    for c in &cases[1..] {
        let case = &c[col("case")];
        let out_dir = work.join(case);
        fs::create_dir_all(&out_dir).unwrap();
        let out = out_dir.join("out");
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_amfsplit"));
        cmd.arg("regions")
            .arg("-i")
            .arg(dir.join(&c[col("reads")]))
            .arg("-p")
            .arg(dir.join(&c[col("panel")]))
            .arg("-o")
            .arg(&out);
        if !c[col("anchors")].is_empty() {
            cmd.arg("-a").arg(dir.join(&c[col("anchors")]));
        }
        cmd.args(c[col("options")].split_whitespace());
        let res = cmd.output().unwrap();
        if !res.status.success() {
            failures.push(format!(
                "{case}: exit {:?}: {}",
                res.status.code(),
                String::from_utf8_lossy(&res.stderr)
            ));
            continue;
        }
        let manifest = tsv(&dir.join("golden").join(case).join("MANIFEST.tsv"));
        for row in &manifest[1..] {
            let (file, digest) = (&row[0], &row[1]);
            let mine = PathBuf::from(format!("{}.{file}", out.display()));
            if !mine.exists() {
                failures.push(format!("{case}: {file} missing"));
            } else if &sha256(&mine) != digest {
                failures.push(format!(
                    "{case}: {file} differs (run parity/stage1/check_parity.py for a field diff)"
                ));
            }
        }
        assert!(
            PathBuf::from(format!("{}.run.json", out.display())).exists(),
            "{case}: run.json missing"
        );
    }
    let _ = fs::remove_dir_all(&work);
    assert!(
        failures.is_empty(),
        "{} mismatch(es):\n{}",
        failures.len(),
        failures.join("\n")
    );
}
