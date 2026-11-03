//! `amfsplit name --reads` on 12 culture sequences (one per species, renamed q_*), with their ITSxRust
//! anchors and cached VT/SH evidence and LSU placements, so no external tool is needed. Checks the
//! whole path, reads to windows (`regions` in-process) to calls, against a committed golden.

use std::fs;
use std::path::Path;
use std::process::Command;

#[test]
fn name_from_reads_matches_its_golden() {
    let smoke = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../parity/naming/name_smoke");
    let out = std::env::temp_dir().join(format!("amfsplit-name-smoke-{}", std::process::id()));
    let _ = fs::remove_dir_all(&out);
    fs::create_dir_all(out.join("work/place")).unwrap();
    fs::copy(smoke.join("work/evidence.json"), out.join("work/evidence.json")).unwrap();
    fs::copy(smoke.join("work/place/lsu_placements.copies.tsv"), out.join("work/place/lsu_placements.copies.tsv")).unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_amfsplit"))
        .args(["name", "-i"]).arg(smoke.join("input/reads.fasta")).arg("-a").arg(smoke.join("input/reads.anchors.tsv"))
        .args(["--db", "unused-everything-cached", "-o"]).arg(&out)
        .status().expect("running amfsplit");
    assert!(status.success(), "amfsplit name failed");
    for f in ["calls.tsv", "rows.tsv"] {
        assert!(fs::read(out.join(f)).unwrap() == fs::read(smoke.join("golden").join(f)).unwrap(), "{f} differs from the golden");
    }
    let _ = fs::remove_dir_all(&out);
}
