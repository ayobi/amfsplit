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
    /// Run the frozen M2b naming reproducibility gate (both holdouts, per-copy golden parity):
    /// Gate 0a, the nearest-reference LSU arm, and Gate 0b step 1, the LSU species-clade arm.
    Gate0(Gate0Args),
    /// Summarize an EPA-ng jplace into the placement table the naming arms read
    /// (lsu_placements.copies.tsv; frozen lsu_placement.py summarize).
    PlaceSummary(PlaceSummaryArgs),
    /// Place LSU windows on the frozen V18 reference: MAFFT --addfragments, EPA-ng, then the
    /// placement summary, with the frozen arguments and thread counts.
    Place(PlaceArgs),
    /// Gate 3: the paper's calls (both culture-only holdouts, species clades plus N6) from the tool's
    /// own BLAST evidence and EPA-ng placements instead of the frozen files.
    Gate3(Gate3Args),
    /// Check MaarjAM, UNITE and the frozen V18 reference against their pins and build the VT and SH
    /// BLAST databases once, with a manifest of the files and tools
    Db(DbArgs),
    /// Name query sequences with the paper's method (species clades plus N6), from the four windows
    /// `amfsplit regions` writes
    Name(NameArgs),
    /// Apply an opt-in Protocol v3 rule to a named-row table (development estimates, dress rehearsal)
    Arm(ArmArgs),
    /// A withholding round: V18 with a taxon's sequences removed, re-placed as the frozen LSU-removal rounds did
    /// (prune, MAFFT L-INS-i, raxml-ng --evaluate, MAFFT --addfragments, EPA-ng, summary)
    Withhold(WithholdArgs),
}

#[derive(Args)]
struct WithholdArgs {
    /// The frozen V18 reference folder (reference.unaligned.fasta, reference_metadata.tsv, backbone.topology.newick)
    #[arg(long)]
    frozen_reference: PathBuf,
    /// A fold of the frozen plan (holdout_plan.json key, e.g. species_diversispora_celata)
    #[arg(long, requires = "plan")]
    fold: Option<String>,
    /// benchmarks/m4_external_holdout/holdout_plan.json
    #[arg(long)]
    plan: Option<PathBuf>,
    /// Or the V18 reference IDs to withhold, comma-separated
    #[arg(long, value_delimiter = ',')]
    removed: Vec<String>,
    /// LSU windows to place (FASTA)
    #[arg(long)]
    queries: PathBuf,
    #[arg(short, long)]
    out: PathBuf,
    /// The rounds' thread count (default 8): MAFFT L-INS-i and EPA-ng use it, raxml-ng and the query MAFFT at most 4
    #[arg(long, default_value_t = 8)]
    threads: usize,
    #[arg(long, default_value = "mafft")]
    mafft: String,
    #[arg(long, default_value = "raxml-ng")]
    raxml_ng: String,
    #[arg(long, default_value = "epa-ng")]
    epa_ng: String,
    /// A frozen round folder (runs/<key>) to compare every stage with
    #[arg(long)]
    compare_with: Option<PathBuf>,
}

fn withhold_cmd(a: WithholdArgs) -> Result<()> {
    let removed = match (&a.fold, &a.plan) {
        (Some(key), Some(plan)) => amfsplit_naming::withhold::fold_removed(plan, key)?,
        _ => a.removed.clone(),
    };
    if removed.is_empty() { anyhow::bail!("nothing to withhold: give --fold with --plan, or --removed"); }
    let o = amfsplit_naming::withhold::WithholdOptions { frozen_reference: a.frozen_reference, removed, queries: a.queries, out: a.out,
        threads: a.threads, mafft: a.mafft, raxml_ng: a.raxml_ng, epa_ng: a.epa_ng, compare_with: a.compare_with };
    let r = amfsplit_naming::withhold::run_withhold(&o)?;
    println!("placement table: {}", r.table.display());
    for (stage, verdict) in &r.stages { println!("  {stage:40} {verdict}"); }
    if let Some(p) = &r.placement { println!("  {:40} {}", "05_summary/lsu_placements.copies.tsv", amfsplit_naming::withhold::verdict(p)); }
    Ok(())
}

#[derive(Args)]
struct ArmArgs {
    /// The rules to apply, comma-separated, run in protocol order: r5 or r5s (the Rhizoglomus motif
    /// step, by VT or with the SSU fallback), r1a or r1b (unnamed cultures as taxa: inside or block),
    /// r2 (the inside gate), r7 (genus-clade fallback, only with an R3 arm), r8 or r8d (novelty on
    /// one marker: flag, or demote), then one R3 flag arm:
    /// r3a, r3a_all (all markers) or r3b (placement)
    #[arg(long, value_delimiter = ',', required = true)]
    rule: Vec<String>,
    /// Named-row table (.tsv or .tsv.gz) with held_cultures, e.g. a v2 golden or a published run's
    /// table; repeat for several (with --out-dir), sharing one distance cache
    #[arg(long, required = true)]
    table: Vec<PathBuf>,
    /// Windows prefix for sequences outside the culture reference (<prefix>.full_its.fasta etc.)
    #[arg(short = 'w', long)]
    windows: Option<PathBuf>,
    /// Output table (one --table)
    #[arg(short = 'o', long, conflicts_with = "out_dir")]
    out: Option<PathBuf>,
    /// Output folder (several --table): each written as <its round folder>__<its file name>
    #[arg(long)]
    out_dir: Option<PathBuf>,
    /// Score frozen against arm under a design: culture, species, rounds, genus, longglodb or eukaryome
    #[arg(long, value_parser = ["culture", "species", "rounds", "genus", "longglodb", "eukaryome"])]
    score: Option<String>,
    /// The EUKARYOME test's locked roles (benchmarks/eukaryome/predictions.predictions.tsv), for --score eukaryome
    #[arg(long)]
    euk_predictions: Option<PathBuf>,
    /// The queries' own placement table (e.g. runs_euk/A/05_summary/lsu_placements.copies.tsv.gz), for R7
    /// on sequences outside the culture reference
    #[arg(long)]
    query_placements: Option<PathBuf>,
    /// R4: the LongGloDB windows prefix (full_its, s58 and the inside markers)
    #[arg(long)]
    r4_longglodb: Option<PathBuf>,
    /// R4: a table with accession, culture and organism for those LongGloDB ASVs
    #[arg(long)]
    r4_longglodb_table: Option<PathBuf>,
    /// R4: the name checks (benchmarks/protocol_v3/r4_name_checks.tsv), verified before fitting
    #[arg(long)]
    r4_name_checks: Option<PathBuf>,
    /// R4: pooled (the rule), or v4 or longglodb (cross-fitted development estimates)
    #[arg(long, default_value = "pooled")]
    r4_fit: String,
    /// R4: the percentile of the refitted limits (decision 6, amended 3 November 2026: 97.5, the lowest
    /// cross-fitted false known among percentiles whose pooled limits demote no more than F's)
    #[arg(long, default_value_t = 97.5)]
    r4_percentile: f64,
    /// R6: the triage table (`amfsplit triage` --out prefix + .triage.tsv)
    #[arg(long)]
    triage: Option<PathBuf>,
    /// S (--rule s, the candidate default): its R5 variant
    #[arg(long, default_value = "r5s")]
    s_r5: String,
    /// S: its R3 variant
    #[arg(long, default_value = "r3b")]
    s_r3: String,
}

fn arm_cmd(a: ArmArgs) -> Result<()> {
    let extra = match &a.windows { Some(p) => amfsplit_naming::read_windows(p)?, None => Default::default() };
    let pairs: Vec<(PathBuf, PathBuf)> = match (&a.out, &a.out_dir) {
        (Some(o), None) if a.table.len() == 1 => vec![(a.table[0].clone(), o.clone())],
        (None, Some(d)) => {
            std::fs::create_dir_all(d)?;
            a.table.iter().map(|t| {
                let round = t.parent().and_then(|p| p.parent()).and_then(|p| p.file_name()).map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                let file = t.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default();
                (t.clone(), d.join(format!("{round}__{file}")))
            }).collect()
        }
        _ => anyhow::bail!("give -o for one table, or --out-dir for several"),
    };
    let rules: Vec<&str> = if a.rule.len() == 1 && a.rule[0] == "s" { amfsplit_naming::s_chain(&a.s_r5, &a.s_r3)? } else { a.rule.iter().map(String::as_str).collect() };
    let r4 = match (&a.r4_longglodb, &a.r4_longglodb_table, &a.r4_name_checks) {
        (Some(w), Some(t), Some(c)) => Some(amfsplit_naming::R4Inputs::prepare(&naming_data(), &a.r4_fit, a.r4_percentile, &amfsplit_naming::read_windows(w)?, t, c)?),
        _ => None,
    };
    let r = amfsplit_naming::run_rules_triage(&naming_data(), &rules, &pairs, &extra, a.query_placements.as_deref(), r4.as_ref(), a.triage.as_deref())?;
    if let Some(design) = &a.score {
        let scores = amfsplit_naming::score_results(&naming_data(), design, &r.results, a.euk_predictions.as_deref())?;
        println!("{:<20} {:<28} {:>14} {:>14}   arm - frozen, 95% CI (points; {} species-cluster replicates, seed {})", "design", "outcome", "frozen", "arm",
                 amfsplit_naming::BOOTSTRAP_REPLICATES, amfsplit_naming::BOOTSTRAP_SEED);
        for s in &scores {
            println!("{:<20} {:<28} {:>14} {:>14}   {:+.2} to {:+.2}", s.design, s.outcome, format!("{}/{}", s.frozen, s.n), format!("{}/{}", s.arm, s.n), s.diff_low, s.diff_high);
        }
        let path = match (&a.out, &a.out_dir) { (Some(o), _) => o.with_extension("scores.tsv"), (None, Some(d)) => d.join("scores.tsv"), _ => unreachable!() };
        std::fs::write(&path, amfsplit_naming::scores_tsv(&scores))?;
        println!("scores written to {}", path.display());
    }
    println!("{} on {} rows in {} rounds of {} table(s): {} calls changed or flagged (R5 D4, R1a U1, R1b U2, R2 N7, R7 G2, R8 N8, R6 T1/T2, R3 G1; {} known calls before)",
             rules.join("+"), r.rows, r.rounds, pairs.len(), r.demoted, r.known_before);
    if r.limits.len() <= 3 { for l in &r.limits { println!("  limits: {l}"); } } else { println!("  limits: {} distinct sets (refitted per round)", r.limits.len()); }
    Ok(())
}

#[derive(Args)]
struct DbArgs {
    /// Output folder (new or empty)
    #[arg(short = 'o', long)]
    out: PathBuf,
    /// MaarjAM FASTA
    #[arg(long)]
    maarjam: PathBuf,
    /// UNITE general FASTA release 10.0, 19 February 2025 (.gz as downloaded)
    #[arg(long)]
    unite: PathBuf,
    /// V18 full FASTA as released (V18_LSUDB_052025.fasta, with its outgroups): the reference alignment is
    /// rebuilt from it exactly as the frozen build did and checked against the pinned one
    #[arg(long, conflicts_with = "reference")]
    v18_full: Option<PathBuf>,
    /// Or an existing frozen V18 reference folder (reference.aligned.fasta, fit.raxml.bestTree,
    /// fit.raxml.bestModel, reference_metadata.tsv)
    #[arg(long)]
    reference: Option<PathBuf>,
    /// MAFFT threads for rebuilding the reference alignment (the result does not depend on it)
    #[arg(long, default_value_t = 8)]
    mafft_threads: usize,
    #[arg(long, default_value = "blastn")]
    blastn: String,
    #[arg(long, default_value = "makeblastdb")]
    makeblastdb: String,
    #[arg(long, default_value = "mafft")]
    mafft: String,
    #[arg(long, default_value = "epa-ng")]
    epa_ng: String,
    /// Accept files that differ from the pinned ones (results are then not the paper's method)
    #[arg(long)]
    allow_unpinned: bool,
}

fn db_cmd(a: DbArgs) -> Result<()> {
    let m = amfsplit_naming::build_db(&amfsplit_naming::DbOptions {
        out: a.out.clone(), maarjam: a.maarjam, unite: a.unite, reference_dir: a.reference, v18_full: a.v18_full,
        shipped_v18: naming_data().join("v18"), mafft_threads: a.mafft_threads,
        blastn: a.blastn, makeblastdb: a.makeblastdb, mafft: a.mafft, epa_ng: a.epa_ng, allow_unpinned: a.allow_unpinned,
    })?;
    for arm in ["VT", "SH"] {
        println!("{arm}: {} references from {} records", m["arms"][arm]["references"], m["arms"][arm]["records_read"]);
    }
    println!("db written to {} ({})", a.out.display(), if m["pinned"].as_bool() == Some(true) { "all files pinned" } else { "UNPINNED files: not the paper's method" });
    Ok(())
}

#[derive(Args)]
struct NameArgs {
    /// Prefix of the `amfsplit regions` output: <prefix>.ssu_flank.fasta, .full_its, .s58, .LSU_LROR-FLR2
    #[arg(short = 'w', long, conflicts_with = "reads")]
    windows: Option<PathBuf>,
    /// Or the reads themselves (FASTA/FASTQ, .gz OK): ITSxRust anchors, then `regions` with the paper's
    /// Stage 1 options and the shipped panel, then naming
    #[arg(short = 'i', long)]
    reads: Option<PathBuf>,
    /// ITSxRust --anchors-tsv for the reads, if already made (skips ITSxRust)
    #[arg(short = 'a', long, requires = "reads")]
    anchors: Option<PathBuf>,
    /// ITSx fungal HMM profile for ITSxRust (F.hmm; the paper's has SHA-256 a8c38c5c...)
    #[arg(long)]
    hmm: Option<PathBuf>,
    /// Read type: hifi (the paper's) or ont (ITSxRust ONT thresholds, 0.20 primer edit budget)
    #[arg(long, default_value = "hifi", value_parser = ["hifi", "ont"])]
    preset: String,
    /// Primer panel for `regions` (default: the shipped amf_legacy_windows.tsv)
    #[arg(long)]
    panel: Option<PathBuf>,
    #[arg(long, default_value = "itsxrust")]
    itsxrust: String,
    /// nhmmer threads for ITSxRust
    #[arg(long, default_value_t = 8)]
    hmmer_cpu: usize,
    /// Folder made by `amfsplit db`
    #[arg(long)]
    db: PathBuf,
    /// Output folder: calls.tsv, rows.tsv, skipped.tsv, run.json, work/ (BLAST and EPA-ng results, reused on rerun)
    #[arg(short = 'o', long)]
    out: PathBuf,
    #[arg(long, default_value_t = 4)]
    blast_threads: usize,
    /// EPA-ng --threads (8 as in the frozen runs)
    #[arg(long, default_value_t = 8)]
    epa_threads: usize,
    /// Override the db manifest's blastn
    #[arg(long)]
    blastn: Option<String>,
    /// Override the db manifest's mafft
    #[arg(long)]
    mafft: Option<String>,
    /// Override the db manifest's epa-ng
    #[arg(long)]
    epa_ng: Option<String>,
    /// Compare the named rows with this table (.tsv or .tsv.gz) on the shared columns
    #[arg(long)]
    compare_with: Option<PathBuf>,
    /// V (Protocol v3, opt-in): VT coverage over the part of each MaarjAM reference the sequence spans
    #[arg(long)]
    vt_spanned: bool,
}

/// The shipped primer panel: beside an installed binary, else in this checkout.
fn default_panel() -> PathBuf {
    if let Some(prefix) = std::env::current_exe().ok().and_then(|e| e.parent().and_then(|b| b.parent()).map(Path::to_path_buf)) {
        let p = prefix.join("share/amfsplit/panels/amf_legacy_windows.tsv");
        if p.is_file() { return p; }
    }
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../panels/amf_legacy_windows.tsv")
}

/// Reads to windows: ITSxRust anchors (unless given), then `regions` with the Stage 1 options of the
/// paper (HiFi: --max-err-frac 0.15 --only-anchored --open-ends; ONT: 0.20 and ITSxRust's ONT
/// thresholds as in run_stage1.sh). Returns the windows prefix.
fn reads_to_windows(a: &NameArgs, reads: &Path) -> Result<PathBuf> {
    let work = a.out.join("work/regions");
    std::fs::create_dir_all(&work)?;
    let anchors = match &a.anchors {
        Some(t) => t.clone(),
        None => {
            let hmm = a.hmm.as_ref().ok_or_else(|| anyhow::anyhow!("--reads needs --hmm (the ITSx fungal profile F.hmm) or --anchors (ITSxRust --anchors-tsv)"))?;
            let itsx = a.out.join("work/itsx");
            std::fs::create_dir_all(&itsx)?;
            let anchors = itsx.join("reads.anchors.tsv");
            let mut cmd = std::process::Command::new(&a.itsxrust);
            cmd.arg("extract").arg("-i").arg(reads).arg("--hmm").arg(hmm).args(["--preset", &a.preset]);
            if a.preset == "ont" { cmd.args(["--inc-e", "1e-3", "--min-anchor-score", "15", "--max-anchor-evalue", "1e-3"]); }
            cmd.args(["--hmmer-cpu", &a.hmmer_cpu.to_string(), "--region", "all", "-o"]).arg(itsx.join("reads"))
                .args(["--output-format", "fasta", "--plain-ids", "--anchors-tsv"]).arg(&anchors)
                .arg("--qc-json").arg(itsx.join("reads.qc.json")).arg("--write-skipped").arg(itsx.join("reads.skipped.fasta"));
            let log = std::fs::File::create(itsx.join("itsxrust.log"))?;
            let status = cmd.stdout(log.try_clone()?).stderr(log).status()
                .with_context(|| format!("starting {} (install ITSxRust 0.3.0, e.g. from Bioconda)", a.itsxrust))?;
            if !status.success() { anyhow::bail!("ITSxRust failed ({status}); see {}", itsx.join("itsxrust.log").display()); }
            use sha2::{Digest, Sha256};
            let hmm_sha = format!("{:x}", Sha256::digest(std::fs::read(hmm)?));
            std::fs::write(itsx.join("hmm.sha256"), format!("{hmm_sha}  {}\n", hmm.display()))?;
            if !hmm_sha.starts_with("a8c38c5c") { eprintln!("note: {} is not the paper's F.hmm (SHA-256 a8c38c5c...)", hmm.display()); }
            anchors
        }
    };
    let prefix = work.join("reads");
    let panel = a.panel.clone().unwrap_or_else(default_panel);
    let err = if a.preset == "ont" { "0.20" } else { "0.15" };
    let argv: Vec<String> = ["amfsplit", "regions", "-i"].iter().map(|s| s.to_string())
        .chain([reads.display().to_string(), "-a".into(), anchors.display().to_string(), "-p".into(), panel.display().to_string(),
                "-o".into(), prefix.display().to_string(), "--max-err-frac".into(), err.into(), "--only-anchored".into(), "--open-ends".into()])
        .collect();
    match Cli::try_parse_from(&argv)?.command {
        Command::Regions(r) => regions(r, &argv)?,
        _ => unreachable!(),
    }
    Ok(prefix)
}

fn name_cmd(a: NameArgs) -> Result<()> {
    let windows = match (&a.windows, &a.reads) {
        (Some(w), None) => w.clone(),
        (None, Some(reads)) => reads_to_windows(&a, reads)?,
        _ => anyhow::bail!("give --windows (the prefix `amfsplit regions` wrote) or --reads"),
    };
    let r = amfsplit_naming::run_name(&amfsplit_naming::NameOptions {
        windows, db: a.db, out: a.out.clone(), fixture_dir: naming_data(),
        blast_threads: a.blast_threads, epa_threads: a.epa_threads, blastn: a.blastn, mafft: a.mafft, epa_ng: a.epa_ng, compare_with: a.compare_with,
        vt_spanned: a.vt_spanned,
    })?;
    println!("{} sequences named, {} skipped (missing windows); {}", r.named, r.skipped, a.out.join("calls.tsv").display());
    for (status, n) in &r.status_counts { println!("  {n:5}  {status}"); }
    if let Some((n, diff)) = &r.comparison {
        if diff.is_empty() { println!("compared with the given table: {n} rows, every shared column identical"); }
        else {
            println!("compared with the given table: {n} rows; differing cells by column:");
            for (c, k) in diff { println!("  {c}: {k}"); }
            std::process::exit(1);
        }
    }
    Ok(())
}

#[derive(Args)]
struct Gate3Args {
    /// MaarjAM FASTA (VT); defaults to $AMFSPLIT_MAARJAM
    #[arg(long)]
    maarjam: Option<PathBuf>,
    /// UNITE general FASTA release 10.0, 19 February 2025 (SH); defaults to $AMFSPLIT_UNITE
    #[arg(long)]
    unite: Option<PathBuf>,
    /// V18 AMF-only LSU FASTA (LSU evidence); defaults to $AMFSPLIT_V18
    #[arg(long)]
    v18: Option<PathBuf>,
    /// Frozen V18 reference folder for the placements (reference.aligned.fasta, fit.raxml.bestTree,
    /// fit.raxml.bestModel, reference_metadata.tsv); defaults to $AMFSPLIT_V18_REFERENCE
    #[arg(long)]
    reference: Option<PathBuf>,
    /// Work folder; BLAST evidence and placements made here are reused by later runs
    #[arg(long, default_value = "gate3_work")]
    work: PathBuf,
    #[arg(long, default_value_t = 4)]
    blast_threads: usize,
    #[arg(long, default_value = "blastn")]
    blastn: String,
    #[arg(long, default_value = "makeblastdb")]
    makeblastdb: String,
    #[arg(long, default_value = "mafft")]
    mafft: String,
    #[arg(long, default_value = "epa-ng")]
    epa_ng: String,
}

fn gate3_cmd(a: Gate3Args) -> Result<()> {
    let env = |k: &str| std::env::var_os(k).filter(|v| !v.is_empty()).map(PathBuf::from);
    let need = |p: Option<PathBuf>, flag: &str, var: &str| p.or_else(|| env(var)).ok_or_else(|| anyhow::anyhow!("gate3 needs --{flag} (or {var})"));
    let databases: std::collections::BTreeMap<String, PathBuf> = [
        ("VT".to_owned(), need(a.maarjam, "maarjam", "AMFSPLIT_MAARJAM")?),
        ("SH".to_owned(), need(a.unite, "unite", "AMFSPLIT_UNITE")?),
        ("LSU".to_owned(), need(a.v18, "v18", "AMFSPLIT_V18")?),
    ].into_iter().collect();
    let report = amfsplit_naming::run_gate3(&amfsplit_naming::Gate3Options {
        fixture_dir: naming_data(),
        databases, reference_dir: need(a.reference, "reference", "AMFSPLIT_V18_REFERENCE")?, work: a.work,
        blast_threads: a.blast_threads, blastn: a.blastn, makeblastdb: a.makeblastdb, mafft: a.mafft, epa_ng: a.epa_ng,
    })?;
    print!("{}", report.summary());
    if !report.passed() { std::process::exit(1); }
    Ok(())
}

#[derive(Args)]
struct PlaceArgs {
    /// LSU windows to place (e.g. LSU_LROR-FLR2.fasta from `amfsplit regions`)
    #[arg(short = 'q', long)]
    queries: PathBuf,
    /// Reference folder: reference.aligned.fasta, fit.raxml.bestTree, fit.raxml.bestModel,
    /// reference_metadata.tsv (the frozen V18 reference; checked by hash)
    #[arg(short = 'r', long)]
    reference: PathBuf,
    /// Output folder (new or empty)
    #[arg(short = 'o', long)]
    out: PathBuf,
    #[arg(long, default_value = "mafft")]
    mafft: String,
    #[arg(long, default_value = "epa-ng")]
    epa_ng: String,
    /// MAFFT --thread (the frozen build used 4)
    #[arg(long, default_value_t = 4)]
    mafft_threads: usize,
    /// EPA-ng --threads (the frozen build used 8; other counts give the same calls but can break
    /// exact ties between edges differently)
    #[arg(long, default_value_t = 8)]
    epa_threads: usize,
    #[arg(long, default_value_t = 0.95)]
    mass: f64,
    /// Use a reference whose files differ from the frozen V18 reference
    #[arg(long)]
    allow_other_reference: bool,
    /// Compare the result with this placement table (tied best edges allowed). For the frozen 903
    /// culture windows the frozen table in this checkout's parity/naming/placement is used by default.
    #[arg(long)]
    compare_with: Option<PathBuf>,
}

fn place_cmd(a: PlaceArgs) -> Result<()> {
    use amfsplit_naming::TableMatch;
    let frozen_table = naming_data().join("placement/lsu_placements.copies.tsv");
    // the frozen table is the default comparison for the frozen query set
    let frozen_queries = std::fs::read(&a.queries).ok()
        .map(|b| { use sha2::{Digest, Sha256}; format!("{:x}", Sha256::digest(&b)) == amfsplit_naming::FROZEN_QUERIES_SHA }).unwrap_or(false);
    let compare_with = a.compare_with.clone().or_else(|| (frozen_queries && frozen_table.is_file()).then(|| frozen_table.clone()));
    let report = amfsplit_naming::run_place(&amfsplit_naming::PlaceOptions {
        reference_dir: a.reference, queries: a.queries, out: a.out.clone(), mafft: a.mafft, epa_ng: a.epa_ng,
        mafft_threads: a.mafft_threads, epa_threads: a.epa_threads, mass: a.mass, allow_other_reference: a.allow_other_reference,
        compare_with: compare_with.clone(),
    })?;
    println!("{} queries placed; {} (SHA-256 {})", report.queries, a.out.join("lsu_placements.copies.tsv").display(), report.table_sha256);
    let what = if report.frozen_inputs { "FROZEN M2b PLACEMENT" } else { "COMPARISON TABLE" };
    match &report.comparison {
        Some(TableMatch::Exact) => println!("{what} REPRODUCED byte for byte"),
        Some(TableMatch::UpToTiedBestEdges(rows)) => {
            println!("{what} REPRODUCED up to tied best edges: in {} row(s) ({}) best_edge names another of the edges tied for the top weight;",
                     rows.len(), rows.join(", "));
            println!("  EPA-ng lists tied placements in run-dependent order. Every other cell is identical, and no naming step reads best_edge.");
        }
        Some(TableMatch::Different(n, first)) => {
            println!("{what} NOT REPRODUCED: {n} row(s) differ beyond tied best edges");
            for d in first { println!("  {d}"); }
            std::process::exit(1);
        }
        None => {}
    }
    Ok(())
}

#[derive(Args)]
struct PlaceSummaryArgs {
    /// EPA-ng output (epa_result.jplace)
    #[arg(long)]
    jplace: PathBuf,
    /// Reference metadata (ref_id, is_amf, label, genus) of the tree in the jplace
    #[arg(long)]
    metadata: PathBuf,
    /// Raw query windows (LSU)
    #[arg(long)]
    queries: PathBuf,
    /// The same queries as aligned to the reference (MAFFT --addfragments --keeplength)
    #[arg(long)]
    aligned: PathBuf,
    /// Output TSV
    #[arg(short = 'o', long)]
    out: PathBuf,
    /// Credible mass
    #[arg(long, default_value_t = 0.95)]
    mass: f64,
}

/// The culture reference and frozen files the naming commands read (the layout of parity/naming):
/// $AMFSPLIT_DATA if set, else share/amfsplit/naming beside an installed binary, else this checkout.
fn naming_data() -> PathBuf {
    if let Some(d) = std::env::var_os("AMFSPLIT_DATA").filter(|v| !v.is_empty()) { return PathBuf::from(d); }
    if let Some(prefix) = std::env::current_exe().ok().and_then(|e| e.parent().and_then(|b| b.parent()).map(Path::to_path_buf)) {
        let d = prefix.join("share/amfsplit/naming");
        if d.join("inputs/benchmark").is_dir() { return d; }
    }
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../parity/naming")
}

fn place_summary_cmd(a: PlaceSummaryArgs) -> Result<()> {
    let rows = amfsplit_naming::summarize_placements(&a.jplace, &a.metadata, &a.queries, &a.aligned, a.mass)?;
    std::fs::write(&a.out, amfsplit_naming::placement_tsv(&rows))?;
    println!("{} queries summarized into {}", rows.len(), a.out.display());
    Ok(())
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
    /// Which method to check: baseline (Gate 0a), species-clades (Gate 0b step 1), n6 (species
    /// clades plus the ITS distance rule, the paper's method; Gate 0b step 2), rounds (the 52
    /// LSU-removal and genus rounds; needs --benchmark), lgdb (the LongGloDB test; needs
    /// --benchmark and --lgdb), euk (the EUKARYOME test; needs --benchmark and --euk) or all.
    #[arg(long, default_value = "all", value_parser = ["baseline", "species-clades", "n6", "rounds", "lgdb", "euk", "placement", "evidence", "all"])]
    policy: String,
    /// amf-novelty-benchmark checkout holding the archived rounds (release 0.3.2 or later);
    /// defaults to $AMFSPLIT_BENCHMARK.
    #[arg(long)]
    benchmark: Option<PathBuf>,
    /// Only the rounds whose key starts with this (e.g. species_glomus).
    #[arg(long)]
    round: Option<String>,
    /// Folder with the LongGloDB test windows and test evidence (never redistributed);
    /// defaults to $AMFSPLIT_LGDB. Used with --benchmark.
    #[arg(long)]
    lgdb: Option<PathBuf>,
    /// Folder with the four EUKARYOME test windows (overlap.win.*.fasta); defaults to $AMFSPLIT_EUK.
    /// Used with --benchmark.
    #[arg(long)]
    euk: Option<PathBuf>,
    /// Runs folder of the LSU-removal test (the --out of external_holdout.py, with each round's
    /// jplace); defaults to $AMFSPLIT_M4_RUNS. Used with --benchmark for --policy placement.
    #[arg(long)]
    m4_runs: Option<PathBuf>,
    /// MaarjAM FASTA (VT evidence, for --policy evidence); defaults to $AMFSPLIT_MAARJAM
    #[arg(long)]
    maarjam: Option<PathBuf>,
    /// UNITE general FASTA release 10.0, 19 February 2025 (SH evidence); defaults to $AMFSPLIT_UNITE
    #[arg(long)]
    unite: Option<PathBuf>,
    /// V18 AMF-only LSU FASTA (LSU evidence); defaults to $AMFSPLIT_V18
    #[arg(long)]
    v18: Option<PathBuf>,
    /// BLAST+ threads for --policy evidence
    #[arg(long, default_value_t = 4)]
    blast_threads: usize,
    /// Only these evidence arms (VT, SH, LSU; repeatable); default all three
    #[arg(long)]
    evidence_arm: Vec<String>,
    /// Check an even sample of about this many culture copies instead of all 903
    #[arg(long)]
    evidence_sample: Option<usize>,
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
        naming_data()
    });
    let mut passed = true;
    if matches!(a.policy.as_str(), "baseline" | "all") {
        let report = amfsplit_naming::run_gate0(&fixtures)?;
        print!("{}", report.summary());
        passed &= report.passed();
    }
    if matches!(a.policy.as_str(), "species-clades" | "all") {
        let report = amfsplit_naming::run_gate0b_species_clades(&fixtures)?;
        print!("{}", report.summary());
        passed &= report.passed();
    }
    if matches!(a.policy.as_str(), "n6" | "all") {
        let report = amfsplit_naming::run_gate0b_n6(&fixtures)?;
        print!("{}", report.summary());
        passed &= report.passed();
    }
    if matches!(a.policy.as_str(), "rounds" | "all") {
        let bench = a.benchmark.clone().or_else(|| std::env::var_os("AMFSPLIT_BENCHMARK").filter(|v| !v.is_empty()).map(PathBuf::from));
        match bench {
            Some(b) => {
                let report = amfsplit_naming::run_gate0b_rounds(&fixtures, &b, a.round.as_deref())?;
                print!("{}", report.summary());
                passed &= report.passed();
            }
            None if a.policy == "rounds" => anyhow::bail!("--policy rounds needs --benchmark or AMFSPLIT_BENCHMARK"),
            None => println!("LSU-removal rounds: SKIPPED (no --benchmark or AMFSPLIT_BENCHMARK)"),
        }
    }
    if matches!(a.policy.as_str(), "lgdb" | "all") {
        let env = |k: &str| std::env::var_os(k).filter(|v| !v.is_empty()).map(PathBuf::from);
        let bench = a.benchmark.clone().or_else(|| env("AMFSPLIT_BENCHMARK"));
        let lg = a.lgdb.clone().or_else(|| env("AMFSPLIT_LGDB"));
        match (bench, lg) {
            (Some(b), Some(l)) => {
                let report = amfsplit_naming::run_gate0b_lgdb(&fixtures, &b, &l)?;
                print!("{}", report.summary());
                passed &= report.passed();
            }
            _ if a.policy == "lgdb" => anyhow::bail!("--policy lgdb needs --benchmark and --lgdb (or AMFSPLIT_BENCHMARK and AMFSPLIT_LGDB)"),
            _ => println!("LongGloDB: SKIPPED (needs --benchmark and --lgdb, or AMFSPLIT_BENCHMARK and AMFSPLIT_LGDB)"),
        }
    }
    if matches!(a.policy.as_str(), "euk" | "all") {
        let env = |k: &str| std::env::var_os(k).filter(|v| !v.is_empty()).map(PathBuf::from);
        let bench = a.benchmark.clone().or_else(|| env("AMFSPLIT_BENCHMARK"));
        let euk = a.euk.clone().or_else(|| env("AMFSPLIT_EUK"));
        match (bench, euk) {
            (Some(b), Some(e)) => {
                let report = amfsplit_naming::run_gate0b_euk(&fixtures, &b, &e)?;
                print!("{}", report.summary());
                passed &= report.passed();
            }
            _ if a.policy == "euk" => anyhow::bail!("--policy euk needs --benchmark and --euk (or AMFSPLIT_BENCHMARK and AMFSPLIT_EUK)"),
            _ => println!("EUKARYOME: SKIPPED (needs --benchmark and --euk, or AMFSPLIT_BENCHMARK and AMFSPLIT_EUK)"),
        }
    }
    if matches!(a.policy.as_str(), "placement" | "all") {
        let env = |k: &str| std::env::var_os(k).filter(|v| !v.is_empty()).map(PathBuf::from);
        let bench = a.benchmark.clone().or_else(|| env("AMFSPLIT_BENCHMARK"));
        let runs = a.m4_runs.clone().or_else(|| env("AMFSPLIT_M4_RUNS"));
        match (bench, runs) {
            (Some(b), Some(r)) => {
                let report = amfsplit_naming::run_placement_parity(&b, &r, a.round.as_deref())?;
                print!("{}", report.summary());
                passed &= report.passed();
            }
            _ if a.policy == "placement" => anyhow::bail!("--policy placement needs --benchmark and --m4-runs (or AMFSPLIT_BENCHMARK and AMFSPLIT_M4_RUNS)"),
            _ => println!("Placement summaries: SKIPPED (needs --benchmark and --m4-runs, or AMFSPLIT_BENCHMARK and AMFSPLIT_M4_RUNS)"),
        }
    }
    if matches!(a.policy.as_str(), "evidence" | "all") {
        let env = |k: &str| std::env::var_os(k).filter(|v| !v.is_empty()).map(PathBuf::from);
        let dbs = [("VT", a.maarjam.clone().or_else(|| env("AMFSPLIT_MAARJAM"))),
                   ("SH", a.unite.clone().or_else(|| env("AMFSPLIT_UNITE"))),
                   ("LSU", a.v18.clone().or_else(|| env("AMFSPLIT_V18")))];
        if dbs.iter().all(|(_, p)| p.is_some()) {
            let dbs: std::collections::BTreeMap<String, PathBuf> = dbs.into_iter().map(|(k, p)| (k.to_owned(), p.unwrap())).collect();
            let work = std::env::temp_dir().join(format!("amfsplit-evidence-{}", std::process::id()));
            let report = amfsplit_naming::run_evidence_parity(&fixtures, &dbs, &work, a.blast_threads, "makeblastdb", "blastn",
                                                              &a.evidence_arm, a.evidence_sample);
            let _ = std::fs::remove_dir_all(&work);
            let report = report?;
            print!("{}", report.summary());
            passed &= report.passed();
        } else if a.policy == "evidence" {
            anyhow::bail!("--policy evidence needs --maarjam, --unite and --v18 (or AMFSPLIT_MAARJAM, AMFSPLIT_UNITE and AMFSPLIT_V18), and makeblastdb and blastn on PATH");
        } else {
            println!("Evidence: SKIPPED (needs --maarjam, --unite and --v18, or AMFSPLIT_MAARJAM, AMFSPLIT_UNITE and AMFSPLIT_V18)");
        }
    }
    if !passed {
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
        Command::PlaceSummary(a) => place_summary_cmd(a),
        Command::Place(a) => place_cmd(a),
        Command::Gate3(a) => gate3_cmd(a),
        Command::Db(a) => db_cmd(a),
        Command::Name(a) => name_cmd(a),
        Command::Arm(a) => arm_cmd(a),
        Command::Withhold(a) => withhold_cmd(a),
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
