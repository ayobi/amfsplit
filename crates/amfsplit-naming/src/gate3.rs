//! Gate 3: the paper's calls from the tool's own evidence and placements.
//!
//! Gate 0 named the 903 culture copies from the frozen evidence files and the frozen placement table.
//! Gate 3 rebuilds both from the culture windows instead: BLAST+ against MaarjAM, UNITE and V18 for
//! the VT, SH and LSU evidence (best units by the frozen first-in-BLAST-order rule), MAFFT and EPA-ng
//! on the frozen V18 reference for the placements. Then it names both culture-only holdouts with
//! species clades plus the ITS distance rule (N6) and compares them with the goldens. What it keeps
//! from the culture reference is the reference itself: the windows, which copy belongs to which
//! culture, and the names. Built evidence and placements are cached in the work folder, so a
//! rerun repeats only the naming.

use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

use anyhow::{anyhow, Result};

use super::evidence::{culture_evidence, load_evidence, save_evidence};
use super::placement::{compare_up_to_ties, summarize_with_ties, TableMatch};
use super::*;

pub struct Gate3Options {
    pub fixture_dir: PathBuf,
    /// VT MaarjAM, SH UNITE (full release), LSU V18 AMF-only
    pub databases: BTreeMap<String, PathBuf>,
    /// frozen V18 reference folder for `place`
    pub reference_dir: PathBuf,
    pub work: PathBuf,
    pub blast_threads: usize,
    pub blastn: String,
    pub makeblastdb: String,
    pub mafft: String,
    pub epa_ng: String,
}

#[derive(Clone, Debug)]
pub struct ArmCheck {
    pub arm: String,
    pub copies: usize,
    pub evidence_equal: usize,
    pub best_equal: usize,
}

#[derive(Clone, Debug)]
pub struct Gate3Report {
    pub evidence_reused: bool,
    pub placement_reused: bool,
    pub arms: Vec<ArmCheck>,
    pub placement: TableMatch,
    pub naming: Gate0Report,
}

impl Gate3Report {
    /// The calls match the goldens, the evidence and best units match the frozen ones, and the
    /// placement table matches the frozen one up to tied best edges.
    pub fn passed(&self) -> bool {
        self.naming.passed()
            && self.arms.iter().all(|a| a.evidence_equal == a.copies && a.best_equal == a.copies)
            && !matches!(self.placement, TableMatch::Different(..))
    }
    pub fn summary(&self) -> String {
        let mut s = String::new();
        for a in &self.arms {
            s.push_str(&format!("evidence {}: {}/{} copies equal to the frozen evidence; best unit {}/{}{}\n", a.arm, a.evidence_equal, a.copies,
                                a.best_equal, a.copies, if self.evidence_reused { " (cached)" } else { "" }));
        }
        let cached = if self.placement_reused { " (cached)" } else { "" };
        s.push_str(&match &self.placement {
            TableMatch::Exact => format!("placement: identical to the frozen table{cached}\n"),
            TableMatch::UpToTiedBestEdges(r) => format!("placement: identical to the frozen table except a tied best edge in {} row(s): {}{cached}\n", r.len(), r.join(", ")),
            TableMatch::Different(n, first) => format!("placement: {n} row(s) differ beyond tied best edges{cached}\n  {}\n", first.join("\n  ")),
        });
        s.push_str(&self.naming.summary());
        s.push_str(if self.passed() {
            "GATE 3, the paper's calls from the tool's own BLAST evidence and EPA-ng placements: PASS\n"
        } else {
            "GATE 3: FAIL (see above)\n"
        });
        s
    }
}

pub fn run_gate3(o: &Gate3Options) -> Result<Gate3Report> {
    let fx = &o.fixture_dir;
    let bench = fx.join("inputs/benchmark");
    fs::create_dir_all(&o.work)?;

    // 1. evidence (cached)
    let cache = o.work.join("evidence.json");
    let evidence_reused = cache.is_file();
    let ev = if evidence_reused { load_evidence(&cache)? } else {
        let blast_dir = o.work.join("blast");
        let ev = culture_evidence(fx, &o.databases, &blast_dir, o.blast_threads, &o.makeblastdb, &o.blastn)?;
        let _ = fs::remove_dir_all(&blast_dir);
        save_evidence(&cache, &ev)?;
        ev
    };

    // 2. placements (cached); the tie-aware comparison re-reads EPA-ng's output either way
    let place_dir = o.work.join("place");
    let table_path = place_dir.join("lsu_placements.copies.tsv");
    let frozen_table = fx.join("placement/lsu_placements.copies.tsv");
    let placement_reused = table_path.is_file();
    if !placement_reused {
        run_place(&PlaceOptions {
            reference_dir: o.reference_dir.clone(), queries: bench.join("ref/fasta/LSU_LROR-FLR2.fasta"), out: place_dir.clone(),
            mafft: o.mafft.clone(), epa_ng: o.epa_ng.clone(), mafft_threads: 4, epa_threads: 8, mass: 0.95,
            allow_other_reference: false, compare_with: Some(frozen_table.clone()),
        })?;
    }
    let (rows, ties) = summarize_with_ties(&place_dir.join("epa/epa_result.jplace"), &o.reference_dir.join("reference_metadata.tsv"),
                                           &place_dir.join("queries.fasta"), &place_dir.join("queries.aligned.fasta"), 0.95)?;
    if placement_tsv(&rows) != fs::read(&table_path)? { return Err(anyhow!("cached placement table does not match its EPA-ng output")); }
    let placement = compare_up_to_ties(&rows, &ties, &frozen_table)?;

    // 3. the frozen data with the rebuilt evidence and best units in place of the frozen ones
    let mut data = FrozenData::load(&bench)?;
    let mut arms = Vec::new();
    for (arm, (queries, best)) in ev {
        let frozen = data.evidence.get_mut(&arm).ok_or_else(|| anyhow!("no {arm} arm in the frozen data"))?;
        let frozen_best = data.best.get_mut(&arm).ok_or_else(|| anyhow!("no {arm} best units in the frozen data"))?;
        if queries.len() != frozen.queries.len() || queries.keys().any(|x| !frozen.queries.contains_key(x)) {
            return Err(anyhow!("{arm}: rebuilt evidence covers other copies than the frozen evidence"));
        }
        let evidence_equal = queries.iter().filter(|(x, q)| {
            let f = &frozen.queries[*x];
            f.anchor_pident == q.anchor_pident && f.unit_pident == q.unit_pident
        }).count();
        let best_equal = best.iter().filter(|(x, b)| frozen_best.get(*x) == Some(b)).count();
        arms.push(ArmCheck { arm: arm.clone(), copies: queries.len(), evidence_equal, best_equal });
        frozen.queries = queries;
        *frozen_best = best;
    }
    let placements = read_row_map(&read_table(&table_path)?);
    let naming = run_policy_on(fx, "epa_species_clades_n6", data, placements)?;
    Ok(Gate3Report { evidence_reused, placement_reused, arms, placement, naming })
}
