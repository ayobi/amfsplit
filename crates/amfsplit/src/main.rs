//! amfsplit: long-read SSU–ITS–LSU metabarcoding of arbuscular mycorrhizal fungi.
//!
//! v0.1 under construction. `regions`, `validate` and `triage` are the Rust port of Stage 1
//! (amfsplit.py, validate_windows.py and triage.py v0.3.1); naming and the decision steps follow.

use std::io::Read;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};
use sha2::{Digest, Sha256};

use amfsplit_regions::pycompat::Json;
use amfsplit_regions::stage1::{Params, check_out_prefix, run};
use amfsplit_regions::triage::{TriageParams, ValidateParams, triage, validate};

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
    /// Best vsearch hit of each window against a reference (validate_windows.py v0.3.1); writes
    /// {out}.hits.tsv, the input of `triage`, and {out}.summary.json.
    Validate(ValidateArgs),
    /// Call each read AMF, unresolved, non_AMF, no_hit or no_window from its best hit, and write
    /// AMF-only copies of the Stage 1 FASTAs (triage.py v0.3.1).
    Triage(TriageArgs),
    /// Run the frozen M2b naming reproducibility gate (both holdouts, per-copy golden parity).
    Gate0(Gate0Args),
}

#[derive(Args)]
struct ValidateArgs {
    /// Window FASTA from `amfsplit regions`
    #[arg(short = 'w', long)]
    windows: PathBuf,
    /// Reference FASTA (for triage: an LSU reference that also holds non-AMF lineages)
    #[arg(short = 'r', long = "ref")]
    reference: PathBuf,
    /// Output prefix
    #[arg(short = 'o', long)]
    out: String,
    /// Minimum identity to report
    #[arg(long, default_value_t = 0.70)]
    id: f64,
    #[arg(long, default_value_t = 4)]
    threads: i64,
    /// vsearch executable (the hits were validated with vsearch 2.27.0)
    #[arg(long, default_value = "vsearch", env = "AMFSPLIT_VSEARCH")]
    vsearch: PathBuf,
}

#[derive(Args)]
struct TriageArgs {
    /// Hits TSV from `amfsplit validate`
    #[arg(long)]
    hits: PathBuf,
    /// Output prefix
    #[arg(long)]
    out: String,
    /// Stage 1 records TSV (adds no_window calls)
    #[arg(long)]
    records: Option<PathBuf>,
    /// Stage 1 output prefix; writes <file>.amf.fasta beside every <prefix>.*.fasta
    #[arg(long)]
    fasta_prefix: Option<String>,
    /// FASTA whose record IDs are the AMF targets; every other target is non-AMF. With the
    /// Delavaux V18 LSU database, search V18_LSUDB_052025.fasta and give
    /// V18_LSUDB_052025_AMFONLY.fasta here (recommended)
    #[arg(long, conflicts_with_all = ["amf_pattern", "outgroups"])]
    amf_reference: Option<PathBuf>,
    /// Regex marking AMF targets; overrides --outgroups (e.g. 'Glomeromycota' for a reference
    /// whose headers carry the lineage)
    #[arg(long)]
    amf_pattern: Option<String>,
    /// Comma-separated non-AMF genera, read from the second underscore-separated field of the
    /// target name. amfsplit v0.3.1 defaulted to Mortierella,Rhodotorula,Exophiala,Citrus, the
    /// outgroups of an older Delavaux release
    #[arg(long)]
    outgroups: Option<String>,
    /// Identity (%) for a confident AMF call (~85 for ONT R9.4.1)
    #[arg(long, default_value_t = 90.0)]
    min_id: f64,
    /// Include unresolved reads in the AMF subsets
    #[arg(long)]
    keep_unresolved: bool,
}

#[derive(Args)]
struct Gate0Args {
    /// Frozen naming parity fixture directory (defaults to this checkout's parity/naming).
    #[arg(long)]
    fixtures: Option<PathBuf>,
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

fn validate_cmd(a: ValidateArgs) -> Result<()> {
    check_out_prefix(&a.out)?;
    let p = ValidateParams {
        windows: a.windows,
        reference: a.reference,
        out: a.out,
        id: a.id,
        threads: a.threads,
        vsearch: a.vsearch,
    };
    print!("{}", validate(&p)?);
    Ok(())
}

fn triage_cmd(a: TriageArgs) -> Result<()> {
    check_out_prefix(&a.out)?;
    let mut p = TriageParams::new(a.hits, a.out);
    p.records = a.records;
    p.fasta_prefix = a.fasta_prefix;
    p.outgroups = a.outgroups;
    p.amf_pattern = a.amf_pattern;
    p.amf_reference = a.amf_reference;
    p.min_id = a.min_id;
    p.keep_unresolved = a.keep_unresolved;
    let report = triage(&p)?;
    print!("{}", report.stdout);
    for note in report.notes {
        eprintln!("{note}");
    }
    Ok(())
}

fn gate0_cmd(a: Gate0Args) -> Result<()> {
    let fixtures = a.fixtures.unwrap_or_else(|| {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../parity/naming")
    });
    let report = amfsplit_naming::run_gate0(&fixtures)?;
    print!("{}", report.summary());
    if !report.passed() {
        anyhow::bail!("frozen naming reproducibility gate failed");
    }
    Ok(())
}

fn main() -> std::process::ExitCode {
    let argv: Vec<String> = std::env::args().collect();
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Regions(a) => regions(a, &argv),
        Command::Validate(a) => validate_cmd(a),
        Command::Triage(a) => triage_cmd(a),
        Command::Gate0(a) => gate0_cmd(a),
    };
    match result {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(e) => {
            // the message and its causes, without anyhow's backtrace (RUST_BACKTRACE is often set)
            eprintln!("error: {e:#}");
            std::process::ExitCode::FAILURE
        }
    }
}
