//! amfsplit: long-read SSU–ITS–LSU metabarcoding of arbuscular mycorrhizal fungi.
//!
//! v0.1 under construction. `regions` is the Rust port of Stage 1 (amfsplit.py v0.3.1); naming
//! and the decision steps follow.

use std::io::Read;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};
use sha2::{Digest, Sha256};

use amfsplit_regions::pycompat::Json;
use amfsplit_regions::stage1::{Params, check_out_prefix, run};

#[derive(Parser)]
#[command(
    name = "amfsplit",
    version,
    about = "Long-read AMF metabarcoding: regions, names and honest uncertainty"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Cut ITS regions (from ITSxRust anchors) and legacy primer windows from long reads
    /// (Stage 1; same options and outputs as amfsplit.py v0.3.1).
    Regions(RegionsArgs),
}

#[derive(Args)]
struct RegionsArgs {
    /// Reads (FASTA/FASTQ, .gz ok)
    #[arg(short = 'i', long)]
    reads: PathBuf,
    /// ITSxRust --anchors-tsv for the same reads
    #[arg(short = 'a', long)]
    anchors: Option<PathBuf>,
    /// Primer-window panel TSV
    #[arg(short = 'p', long)]
    panel: PathBuf,
    /// Output prefix
    #[arg(short = 'o', long)]
    out: String,
    /// Max edits per primer as a fraction of its length
    #[arg(long, default_value_t = 0.15)]
    max_err_frac: f64,
    /// Include primer sites in windows (default: primer-trimmed)
    #[arg(long)]
    keep_primers: bool,
    /// When one primer site lies outside the amplicon, cut the window to the read end and flag
    /// it open_5p/open_3p (e.g. NS31 windows on nested amplicons that start downstream of NS31)
    #[arg(long)]
    open_ends: bool,
    /// Skip reads ITSxRust could not anchor (e.g. off-target host rDNA)
    #[arg(long)]
    only_anchored: bool,
    /// Search each primer inside its ITSxRust region widened by this many nt
    #[arg(long, default_value_t = 50, allow_negative_numbers = true)]
    region_margin: i64,
    /// Search every primer over the whole read (the v0.2 behaviour)
    #[arg(long)]
    no_region_limit: bool,
    /// Skip the concatemer screen
    #[arg(long)]
    no_concat_screen: bool,
    /// Different primers that must show a repeat (or an inverted site) to flag a read
    #[arg(long, default_value_t = 2)]
    concat_min_primers: i64,
    /// Min distance between two sites of one primer to count as two copies
    #[arg(long, default_value_t = 150)]
    concat_min_sep: i64,
    /// Max edits per base between the two supposed copies (unrelated DNA is ~0.5)
    #[arg(long, default_value_t = 0.40)]
    concat_sim: f64,
    /// Primers agreeing on one spacing that need no copy-similarity check
    #[arg(long, default_value_t = 3)]
    concat_strong: i64,
    /// Primers must agree on the copy spacing within this fraction (min 40 nt)
    #[arg(long, default_value_t = 0.05)]
    concat_tol: f64,
    /// Max edits per primer for the concatemer screen (default: --max-err-frac)
    #[arg(long)]
    concat_err_frac: Option<f64>,
    /// Write windows and regions of flagged concatemers too (they stay flagged in records.tsv)
    #[arg(long)]
    keep_concatemers: bool,
}

fn sha256_of(path: &Path) -> Result<(String, u64)> {
    let mut f =
        std::fs::File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let mut h = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    let mut total = 0u64;
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        total += n as u64;
        h.update(&buf[..n]);
    }
    Ok((format!("{:x}", h.finalize()), total))
}

/// `{out}.run.json`: what produced the outputs. Kept apart from summary.json, which stays
/// byte-identical to the prototype's.
fn write_run_json(p: &Params, argv: &[String], outputs: &[PathBuf]) -> Result<PathBuf> {
    let file_entry = |path: &Path| -> Result<Json> {
        let (digest, bytes) = sha256_of(path)?;
        Ok(Json::obj(vec![
            ("path", Json::Str(path.display().to_string())),
            ("sha256", Json::Str(digest)),
            ("bytes", Json::Int(bytes as i64)),
        ]))
    };
    let mut inputs = vec![
        ("reads".to_string(), file_entry(&p.reads)?),
        ("panel".to_string(), file_entry(&p.panel)?),
    ];
    if let Some(a) = &p.anchors {
        inputs.push(("anchors".to_string(), file_entry(a)?));
    }
    let mut outs = Vec::new();
    for o in outputs {
        let (digest, bytes) = sha256_of(o)?;
        let name = o
            .file_name()
            .map_or(String::new(), |n| n.to_string_lossy().into_owned());
        outs.push((
            name,
            Json::obj(vec![
                ("sha256", Json::Str(digest)),
                ("bytes", Json::Int(bytes as i64)),
            ]),
        ));
    }
    let run = Json::obj(vec![
        ("tool", Json::Str("amfsplit".into())),
        ("tool_version", Json::Str(env!("CARGO_PKG_VERSION").into())),
        ("command", Json::Str("regions".into())),
        (
            "stage1_algorithm_version",
            Json::Str(amfsplit_regions::STAGE1_ALGORITHM_VERSION.into()),
        ),
        (
            "edlib",
            Json::Str("vendored C++ core of Python edlib 1.3.9.post1".into()),
        ),
        ("argv", Json::Str(argv.join(" "))),
        ("inputs", Json::Obj(inputs)),
        ("outputs", Json::Obj(outs)),
    ]);
    let path = p.output_path("run.json");
    std::fs::write(&path, run.dump_indent1() + "\n")
        .with_context(|| format!("cannot write {}", path.display()))?;
    Ok(path)
}

fn regions(a: RegionsArgs, argv: &[String]) -> Result<()> {
    check_out_prefix(&a.out)?;
    let mut p = Params::new(a.reads, a.anchors, a.panel, a.out);
    p.max_err_frac = a.max_err_frac;
    p.keep_primers = a.keep_primers;
    p.open_ends = a.open_ends;
    p.only_anchored = a.only_anchored;
    p.region_margin = a.region_margin;
    p.no_region_limit = a.no_region_limit;
    p.no_concat_screen = a.no_concat_screen;
    p.concat_min_primers = a.concat_min_primers;
    p.concat_min_sep = a.concat_min_sep;
    p.concat_sim = a.concat_sim;
    p.concat_strong = a.concat_strong;
    p.concat_tol = a.concat_tol;
    p.concat_err_frac = a.concat_err_frac;
    p.keep_concatemers = a.keep_concatemers;
    let report = run(&p)?;
    write_run_json(&p, argv, &report.outputs)?;
    eprint!("{}", report.console);
    Ok(())
}

fn main() -> Result<()> {
    let argv: Vec<String> = std::env::args().collect();
    let cli = Cli::parse();
    match cli.command {
        Command::Regions(a) => regions(a, &argv),
    }
}
