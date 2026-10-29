//! Standalone implementation of the Python N6 full-ITS novelty rule.
//!
//! This module deliberately does not modify the existing evaluator.
//! Its purpose is to establish parity with the frozen Python reference
//! before N6 is integrated into the naming workflow.

use std::collections::BTreeMap;

use anyhow::{anyhow, bail, Result};

use super::{
    percent_distance, py_float, py_round_float, quantiles_exclusive, row_value, FrozenData, Row,
    StringSet,
};

type CultureDistances = BTreeMap<String, f64>;
type CopyDistances = BTreeMap<String, CultureDistances>;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Pool {
    All,
    Main,
}

impl Pool {
    fn name(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::Main => "main",
        }
    }
}

#[derive(Debug)]
pub(super) struct ItsRule {
    /// Minimum distance from each copy to each culture, for each pool.
    dmin_all: CopyDistances,
    dmin_main: CopyDistances,

    /// Threshold cache, keyed by the exact training-culture set.
    threshold_cache: BTreeMap<StringSet, (Option<f64>, usize)>,
}

impl ItsRule {
    /// Construct the pairwise ITS distance cache.
    ///
    /// Accession ordering and pair orientation intentionally match Python:
    /// sort accessions, calculate each unordered pair once, and store the
    /// resulting distance symmetrically.
    pub(super) fn new(data: &FrozenData) -> Result<Self> {
        let its = data
            .fa
            .get("full_its")
            .ok_or_else(|| anyhow!("missing full_its reference sequences"))?;

        let its_ok: Vec<&String> = data
            .copies
            .keys()
            .filter(|x| its.contains_key(*x))
            .collect();

        let mains: StringSet = its_ok
            .iter()
            .filter(|x| {
                data.copies
                    .get(x.as_str())
                    .is_some_and(|meta| meta.copy_class == "main")
            })
            .map(|x| (*x).clone())
            .collect();

        let mut dmin_all: CopyDistances = BTreeMap::new();
        let mut dmin_main: CopyDistances = BTreeMap::new();

        for i in 0..its_ok.len() {
            let x = its_ok[i];
            let seq_x = &its[x.as_str()];

            for y in its_ok.iter().skip(i + 1) {
                let seq_y = &its[y.as_str()];

                // Preserve Python's argument order: dist(its[x], its[y]).
                let d = percent_distance(seq_x, seq_y, -1)?
                    .ok_or_else(|| anyhow!("unbounded ITS distance failed for {x} and {y}"))?;

                let culture_x = &data.copies[x.as_str()].culture;
                let culture_y = &data.copies[y.as_str()].culture;

                // Store the same distance in both directions.
                Self::update_min(&mut dmin_all, x.as_str(), culture_y, d);
                Self::update_min(&mut dmin_all, y.as_str(), culture_x, d);

                // The main pool includes a comparison only when the
                // reference target copy is main-type.
                if mains.contains(*y) {
                    Self::update_min(&mut dmin_main, x.as_str(), culture_y, d);
                }

                if mains.contains(x) {
                    Self::update_min(&mut dmin_main, y.as_str(), culture_x, d);
                }
            }
        }

        Ok(Self {
            dmin_all,
            dmin_main,
            threshold_cache: BTreeMap::new(),
        })
    }

    fn update_min(distances: &mut CopyDistances, query: &str, culture: &str, distance: f64) {
        let by_culture = distances.entry(query.to_owned()).or_default();

        let current = by_culture
            .entry(culture.to_owned())
            .or_insert(f64::INFINITY);

        if distance < *current {
            *current = distance;
        }
    }

    /// Calculate the Python-compatible within-species training threshold.
    ///
    /// Returns (threshold, number of distances contributing to it).
    /// A threshold is unavailable when fewer than two distances exist.
    pub(super) fn threshold(
        &mut self,
        data: &FrozenData,
        training: &StringSet,
        quantile: usize,
    ) -> Result<(Option<f64>, usize)> {
        if !(1..=99).contains(&quantile) {
            bail!("N6 quantile must be between 1 and 99; got {quantile}");
        }

        if let Some(cached) = self.threshold_cache.get(training) {
            // The frozen Python reference uses a fixed quantile per rule
            // instance. This implementation accepts a quantile argument,
            // so only cache the default q99 threshold.
            if quantile == 99 {
                return Ok(*cached);
            }
        }

        let mut values: Vec<f64> = Vec::new();
        let mut used: StringSet = StringSet::new();

        for (species, cultures) in &data.by_species {
            let tc: StringSet = cultures.intersection(training).cloned().collect();

            if tc.len() < 2 {
                continue;
            }

            for culture in &tc {
                let Some(copies) = data.by_culture.get(culture) else {
                    continue;
                };

                for y in copies {
                    if data.copies[y].copy_class != "main" {
                        continue;
                    }

                    let Some(distances) = self.dmin_main.get(y) else {
                        continue;
                    };

                    let mut nearest: Option<f64> = None;

                    for other_culture in tc.iter().filter(|c| *c != culture) {
                        if let Some(&d) = distances.get(other_culture) {
                            if nearest.is_none_or(|old| d < old) {
                                nearest = Some(d);
                            }
                        }
                    }

                    if let Some(d) = nearest {
                        values.push(d);
                        used.insert(culture.clone());
                    }
                }
            }

            // Each value is based only on cultures belonging to this species.
            let _ = species;
        }

        if !used.is_subset(training) {
            bail!("held-out culture entered the ITS threshold");
        }

        let n = values.len();

        let threshold = if n >= 2 {
            Some(quantiles_exclusive(&values, 100)[quantile - 1])
        } else {
            None
        };

        if quantile == 99 {
            self.threshold_cache
                .insert(training.clone(), (threshold, n));
        }

        Ok((threshold, n))
    }

    /// Find the nearest training-culture distance, excluding the query's own
    /// culture. The cache already contains the minimum per culture.
    pub(super) fn nearest(
        &self,
        data: &FrozenData,
        accession: &str,
        training: &StringSet,
        pool: Pool,
    ) -> Result<Option<f64>> {
        let own_culture = &data
            .copies
            .get(accession)
            .ok_or_else(|| anyhow!("unknown copy {accession}"))?
            .culture;

        let distances = match pool {
            Pool::All => &self.dmin_all,
            Pool::Main => &self.dmin_main,
        };

        let Some(by_culture) = distances.get(accession) else {
            return Ok(None);
        };

        Ok(by_culture
            .iter()
            .filter(|(culture, _)| training.contains(*culture) && *culture != own_culture)
            .map(|(_, distance)| *distance)
            .min_by(f64::total_cmp))
    }

    /// Apply N6 to cloned or original output rows.
    ///
    /// `rec` must contain the original reconciliation rows, keyed by
    /// accession, so placement comes from the original genus call.
    pub(super) fn apply(
        &mut self,
        data: &FrozenData,
        rows: &mut [Row],
        rec: &BTreeMap<String, Row>,
        training: &StringSet,
        pool: Pool,
        quantile: usize,
    ) -> Result<usize> {
        let (threshold, n) = self.threshold(data, training, quantile)?;
        let mut fired = 0;

        for row in rows {
            let accession = row_value(row, "accession")?.to_owned();

            let distance = self.nearest(data, &accession, training, pool)?;

            row.insert(
                "its_nearest_training_pct".into(),
                distance
                    .map(|d| py_float(py_round_float(d, 3)))
                    .unwrap_or_default(),
            );

            row.insert(
                "its_threshold_pct".into(),
                threshold
                    .map(|t| py_float(py_round_float(t, 3)))
                    .unwrap_or_default(),
            );

            row.insert("its_threshold_n".into(), n.to_string());
            row.insert("its_pool".into(), pool.name().to_owned());

            let status = row_value(row, "status")?;
            let hit = should_fire(distance, threshold, status);

            row.insert("its_rule_fired".into(), if hit { "1" } else { "0" }.into());

            if !hit {
                continue;
            }

            let d = distance.expect("a firing rule has a distance");
            let t = threshold.expect("a firing rule has a threshold");

            let original = rec
                .get(&accession)
                .ok_or_else(|| anyhow!("no reconciliation row for {accession}"))?;

            let genus = row_value(original, "call_genera")?.to_owned();

            row.insert("status".into(), "novel species, genus known".into());
            row.insert("placement".into(), genus);
            row.insert("novelty_claim".into(), "candidate_only".into());

            let reason = format!(
                "N6 full ITS {:.1}% from the nearest training culture, above {:.1}% (the {}th percentile of within-species distance between training cultures)",
                d, t, quantile
            );

            let existing = row_value(row, "reasons")?.to_owned();

            let combined = if existing.is_empty() {
                reason
            } else {
                format!("{existing}; {reason}")
            };

            row.insert("reasons".into(), combined);
            fired += 1;
        }

        Ok(fired)
    }
}

fn should_fire(distance: Option<f64>, threshold: Option<f64>, status: &str) -> bool {
    distance.is_some_and(|d| threshold.is_some_and(|t| d > t)) && status.starts_with("known")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    fn fixture() -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../parity/naming")
            .canonicalize()
            .expect("naming parity fixture should exist")
    }

    #[test]
    fn firing_condition_matches_python_strict_threshold_rule() {
        assert!(should_fire(Some(5.1), Some(5.0), "known species"));
        assert!(!should_fire(Some(5.0), Some(5.0), "known species"));
        assert!(!should_fire(Some(4.9), Some(5.0), "known species"));
        assert!(!should_fire(None, Some(5.0), "known species"));
        assert!(!should_fire(Some(5.1), None, "known species"));
        assert!(!should_fire(Some(5.1), Some(5.0), "novel lineage"));
        assert!(!should_fire(Some(5.1), Some(5.0), "unplaced"));
    }

    fn json_string_array(value: &serde_json::Value, context: &str) -> Vec<String> {
        value
            .as_array()
            .unwrap_or_else(|| panic!("{context}: expected an array"))
            .iter()
            .map(|item| {
                item.as_str()
                    .unwrap_or_else(|| panic!("{context}: expected string values"))
                    .to_owned()
            })
            .collect()
    }

    fn assert_optional_float_close(expected: Option<f64>, actual: Option<f64>, context: &str) {
        match (expected, actual) {
            (Some(expected), Some(actual)) => assert!(
                (expected - actual).abs() < 1e-9,
                "{context}: Python expected {expected:.15}, Rust got {actual:.15}"
            ),
            (None, None) => {}
            (expected, actual) => {
                panic!("{context}: Python expected {expected:?}, Rust got {actual:?}")
            }
        }
    }

    #[test]
    fn all_37_species_folds_match_python_reference() {
        use serde_json::Value;
        use std::fs;

        let fixture_root = fixture();
        let parity_dir = fixture_root.join("n6_fold_parity");
        let manifest_path = parity_dir.join("fold_manifest.json");
        let reference_path = parity_dir.join("python_reference.json");

        let manifest_bytes = fs::read(&manifest_path).expect("frozen N6 fold manifest should load");
        let reference_bytes =
            fs::read(&reference_path).expect("frozen N6 Python reference should load");

        // Lock the fixtures so the expected results cannot change silently.
        let manifest_sha = super::super::sha256_bytes(&manifest_bytes);
        let reference_sha = super::super::sha256_bytes(&reference_bytes);

        assert_eq!(
            manifest_sha, "5d538237a8d0a8ba1446fc586f1f3ab04d91abd44a01669f7475abf6a96978d9",
            "N6 fold manifest changed; update the lock only after deliberate review"
        );
        assert_eq!(
            reference_sha, "0508049721ed5098ef5ceb594a4b2a67782c92c8c9a614fadd08154632696bbe",
            "N6 Python reference changed; update the lock only after deliberate review"
        );

        let manifest: Value =
            serde_json::from_slice(&manifest_bytes).expect("fold manifest should be valid JSON");
        let reference: Value = serde_json::from_slice(&reference_bytes)
            .expect("Python reference should be valid JSON");

        assert_eq!(
            reference
                .get("fold_manifest_sha256")
                .and_then(Value::as_str),
            Some(manifest_sha.as_str()),
            "Python reference must identify this exact fold manifest"
        );

        assert_eq!(manifest["fold_count"].as_u64(), Some(37));
        assert_eq!(reference["fold_count"].as_u64(), Some(37));
        assert_eq!(reference["total_query_fold_instances"].as_u64(), Some(850));

        let manifest_folds = manifest["folds"]
            .as_array()
            .expect("manifest folds should be an array");
        let reference_folds = reference["folds"]
            .as_array()
            .expect("reference folds should be an array");

        assert_eq!(manifest_folds.len(), 37);
        assert_eq!(manifest_folds.len(), reference_folds.len());

        let data = FrozenData::load(fixture_root.join("inputs/benchmark"))
            .expect("frozen benchmark data should load");

        assert_eq!(
            data.cultures.len(),
            reference["reference_culture_count"].as_u64().unwrap() as usize
        );
        assert_eq!(
            data.copies.len(),
            reference["reference_accession_count"].as_u64().unwrap() as usize
        );

        // Build the expensive sequence-distance cache once, then reuse it
        // across every training fold.
        let mut rule = ItsRule::new(&data).expect("N6 ITS distance cache should build");

        let mut total_queries = 0usize;
        let mut missing_all = 0usize;
        let mut missing_main = 0usize;

        for (manifest_fold, reference_fold) in manifest_folds.iter().zip(reference_folds.iter()) {
            let target = manifest_fold["target"]
                .as_str()
                .expect("fold target should be a string");

            assert_eq!(
                reference_fold["target"].as_str(),
                Some(target),
                "fold target differs between manifest and Python reference"
            );
            assert_eq!(
                manifest_fold["fold_signature_sha256"].as_str(),
                reference_fold["fold_signature_sha256"].as_str(),
                "{target}: fold signature mismatch"
            );

            let held: StringSet = json_string_array(
                &manifest_fold["held_cultures"],
                &format!("{target} held cultures"),
            )
            .into_iter()
            .collect();

            let training: StringSet = json_string_array(
                &manifest_fold["training_cultures"],
                &format!("{target} training cultures"),
            )
            .into_iter()
            .collect();

            assert!(
                held.iter()
                    .all(|culture| data.cultures.contains_key(culture)),
                "{target}: manifest contains an unknown held-out culture"
            );

            let expected_training: StringSet = data
                .cultures
                .keys()
                .filter(|culture| !held.contains(*culture))
                .cloned()
                .collect();

            assert_eq!(
                training, expected_training,
                "{target}: training cultures are not exactly the complement of held cultures"
            );

            assert_eq!(
                reference_fold["training_culture_count"].as_u64(),
                Some(training.len() as u64),
                "{target}: training-culture count differs"
            );

            let (threshold_actual, n_actual) = rule
                .threshold(&data, &training, 99)
                .unwrap_or_else(|e| panic!("{target}: threshold calculation failed: {e:#}"));

            let threshold_expected = reference_fold["threshold_pct"].as_f64();

            assert_optional_float_close(
                threshold_expected,
                threshold_actual,
                &format!("{target} threshold"),
            );

            assert_eq!(
                reference_fold["threshold_n"].as_u64(),
                Some(n_actual as u64),
                "{target}: threshold observation count differs"
            );

            let manifest_queries = json_string_array(
                &manifest_fold["query_ids"],
                &format!("{target} query accessions"),
            );

            let reference_queries = reference_fold["queries"]
                .as_array()
                .unwrap_or_else(|| panic!("{target}: reference queries should be an array"));

            assert_eq!(
                manifest_fold["query_ids"].as_array().unwrap().len(),
                reference_queries.len(),
                "{target}: query count differs between fixtures"
            );
            assert_eq!(
                reference_fold["query_count"].as_u64(),
                Some(reference_queries.len() as u64),
                "{target}: recorded query count is inconsistent"
            );

            let reference_accessions: Vec<String> = reference_queries
                .iter()
                .map(|query| {
                    query["accession"]
                        .as_str()
                        .unwrap_or_else(|| panic!("{target}: query accession should be a string"))
                        .to_owned()
                })
                .collect();

            assert_eq!(
                manifest_queries, reference_accessions,
                "{target}: query accessions or their order differ"
            );

            total_queries += reference_queries.len();

            for query in reference_queries {
                let accession = query["accession"]
                    .as_str()
                    .expect("query accession should be a string");

                let copy = data.copies.get(accession).unwrap_or_else(|| {
                    panic!("{target}: query accession {accession} is absent from Rust reference")
                });

                let expected_culture = query["culture"]
                    .as_str()
                    .expect("query culture should be a string");

                assert_eq!(
                    copy.culture, expected_culture,
                    "{target}/{accession}: query culture differs"
                );
                assert!(
                    held.contains(&copy.culture),
                    "{target}/{accession}: query culture is not held out"
                );

                let expected_all = query["distance_all_pct"].as_f64();
                let actual_all = rule
                    .nearest(&data, accession, &training, Pool::All)
                    .unwrap_or_else(|e| {
                        panic!("{target}/{accession}: all-pool calculation failed: {e:#}")
                    });

                if actual_all.is_none() {
                    missing_all += 1;
                }

                assert_optional_float_close(
                    expected_all,
                    actual_all,
                    &format!("{target}/{accession} all-pool distance"),
                );

                let expected_main = query["distance_main_pct"].as_f64();
                let actual_main = rule
                    .nearest(&data, accession, &training, Pool::Main)
                    .unwrap_or_else(|e| {
                        panic!("{target}/{accession}: main-pool calculation failed: {e:#}")
                    });

                if actual_main.is_none() {
                    missing_main += 1;
                }

                assert_optional_float_close(
                    expected_main,
                    actual_main,
                    &format!("{target}/{accession} main-pool distance"),
                );
            }
        }

        assert_eq!(total_queries, 850);
        assert_eq!(
            missing_all,
            reference["missing_distance_all_count"].as_u64().unwrap() as usize,
            "all-pool missing-distance count differs"
        );
        assert_eq!(
            missing_main,
            reference["missing_distance_main_count"].as_u64().unwrap() as usize,
            "main-pool missing-distance count differs"
        );

        println!(
            "N6 fold parity PASS: 37 folds, {total_queries} queries, \
             thresholds and counts matched, all/main distances matched."
        );
    }

    #[test]
    fn full_training_threshold_matches_frozen_python_reference() {
        let data = FrozenData::load(fixture().join("inputs/benchmark"))
            .expect("frozen benchmark data should load");

        // Use every culture as the training set.
        let training: StringSet = data.cultures.keys().cloned().collect();

        let mut rule = ItsRule::new(&data).expect("N6 ITS distance cache should build");

        let (threshold, n) = rule
            .threshold(&data, &training, 99)
            .expect("N6 threshold should calculate");

        // Expected values from the frozen Python implementation.
        assert_eq!(
            n, 755,
            "N6 should use the same number of distance observations as Python"
        );

        let threshold = threshold.expect("threshold should be available");

        assert!(
            (threshold - 5.992549889135255).abs() < 1e-9,
            "Rust threshold {threshold:.15} differs from the Python reference"
        );
    }

    #[test]
    fn ambispora_leptoticha_fold_matches_frozen_python_reference() {
        let data = FrozenData::load(fixture().join("inputs/benchmark"))
            .expect("frozen benchmark data should load");

        let held_out: StringSet = [
            "Ambispora leptoticha | MAFF520090_2501A".to_owned(),
            "Ambispora leptoticha | ON205A".to_owned(),
        ]
        .into_iter()
        .collect();

        let training: StringSet = data
            .cultures
            .keys()
            .filter(|culture| !held_out.contains(*culture))
            .cloned()
            .collect();

        assert_eq!(training.len(), 149);

        let mut rule = ItsRule::new(&data).expect("N6 ITS distance cache should build");

        let (threshold, n) = rule
            .threshold(&data, &training, 99)
            .expect("fold-specific N6 threshold should calculate");

        assert_eq!(n, 745);

        let threshold = threshold.expect("threshold should be available");
        assert!(
            (threshold - 5.993880266075388).abs() < 1e-9,
            "Rust fold threshold {threshold:.15} differs from Python"
        );

        let expected: [(&str, f64); 10] = [
            ("PX207101.1", 10.379241516966),
            ("PX207102.1", 10.379241516966),
            ("PX207103.1", 10.179640718563),
            ("PX207104.1", 9.218436873747),
            ("PX207105.1", 9.819639278557),
            ("PX207106.1", 9.819639278557),
            ("PX215138.1", 9.780439121756),
            ("PX215139.1", 10.379241516966),
            ("PX215140.1", 10.179640718563),
            ("PX215141.1", 9.677419354839),
        ];

        for (accession, expected_distance) in expected {
            let distance = rule
                .nearest(&data, accession, &training, Pool::All)
                .expect("nearest-training distance should calculate")
                .unwrap_or_else(|| panic!("missing distance for {accession}"));

            assert!(
                (distance - expected_distance).abs() < 1e-9,
                "{accession}: Rust distance {distance:.12} differs from Python {expected_distance:.12}"
            );

            let main_distance = rule
                .nearest(&data, accession, &training, Pool::Main)
                .expect("main-pool nearest-training distance should calculate")
                .unwrap_or_else(|| panic!("missing main-pool distance for {accession}"));

            assert!(
                (main_distance - expected_distance).abs() < 1e-9,
                "{accession}: Rust main-pool distance {main_distance:.12} differs from Python {expected_distance:.12}"
            );
        }
    }
    #[test]
    fn apply_matches_ambispora_python_reference_and_preserves_baseline() {
        use super::super::{evaluate_novelty, reconcile, Trainer, ARM_NAMES};

        let data = FrozenData::load(fixture().join("inputs/benchmark"))
            .expect("frozen benchmark data should load");

        let held_out: StringSet = [
            "Ambispora leptoticha | MAFF520090_2501A".to_owned(),
            "Ambispora leptoticha | ON205A".to_owned(),
        ]
        .into_iter()
        .collect();

        let training: StringSet = data
            .cultures
            .keys()
            .filter(|culture| !held_out.contains(*culture))
            .cloned()
            .collect();

        let query_ids: Vec<String> = [
            "PX207101.1",
            "PX207102.1",
            "PX207103.1",
            "PX207104.1",
            "PX207105.1",
            "PX207106.1",
            "PX215138.1",
            "PX215139.1",
            "PX215140.1",
            "PX215141.1",
        ]
        .into_iter()
        .map(str::to_owned)
        .collect();

        let trainers: BTreeMap<String, Trainer> = ARM_NAMES
            .iter()
            .map(|arm| Ok(((*arm).to_owned(), Trainer::new(&data, arm)?)))
            .collect::<Result<_>>()
            .expect("trainers should build");

        let fits: BTreeMap<String, super::super::Fit> = ARM_NAMES
            .iter()
            .map(|arm| Ok(((*arm).to_owned(), trainers[*arm].fit(&data, &training)?)))
            .collect::<Result<_>>()
            .expect("training fits should build");

        let mut rec = BTreeMap::new();
        for accession in &query_ids {
            rec.insert(
                accession.clone(),
                reconcile(&data, accession, &fits, &held_out)
                    .expect("query reconciliation should succeed"),
            );
        }

        let baseline = evaluate_novelty(&data, &query_ids, &fits, &held_out, &rec)
            .expect("baseline novelty evaluation should succeed");

        assert_eq!(baseline.len(), 10);
        assert!(baseline
            .iter()
            .all(|row| { row.get("status").map(String::as_str) == Some("known species") }));

        // Apply N6 to a separate copy: baseline results must remain untouched.
        let mut adjusted = baseline.clone();
        let mut rule = ItsRule::new(&data).expect("N6 rule should build");

        let fired = rule
            .apply(&data, &mut adjusted, &rec, &training, Pool::All, 99)
            .expect("N6 application should succeed");

        assert_eq!(fired, 10);
        assert_eq!(
            baseline
                .iter()
                .filter(|row| { row.get("status").map(String::as_str) == Some("known species") })
                .count(),
            10
        );

        for row in &adjusted {
            let accession = row_value(row, "accession").unwrap();
            assert_eq!(
                row.get("status").map(String::as_str),
                Some("novel species, genus known"),
                "unexpected status for {accession}"
            );
            assert_eq!(
                row.get("novelty_claim").map(String::as_str),
                Some("candidate_only"),
                "unexpected novelty claim for {accession}"
            );
            assert_eq!(
                row.get("its_rule_fired").map(String::as_str),
                Some("1"),
                "N6 should fire for {accession}"
            );
            assert_eq!(row.get("its_pool").map(String::as_str), Some("all"));

            let original_rec = &rec[accession];
            assert_eq!(
                row.get("placement"),
                original_rec.get("call_genera"),
                "placement should use call_genera for {accession}"
            );

            let reason = row.get("reasons").map(String::as_str).unwrap_or("");
            assert!(
                reason.contains("N6 full ITS "),
                "N6 reason missing for {accession}: {reason:?}"
            );

            assert_eq!(row.get("its_threshold_n").map(String::as_str), Some("745"));

            let nearest: f64 = row["its_nearest_training_pct"]
                .parse()
                .expect("nearest distance should be numeric");
            let threshold: f64 = row["its_threshold_pct"]
                .parse()
                .expect("threshold should be numeric");

            assert!(nearest > threshold, "{accession}: {nearest} <= {threshold}");
        }

        // Ineligible statuses must not be demoted, even when their
        // ITS distance exceeds the threshold.
        let mut ineligible = baseline.clone();
        let target = "PX207101.1";

        let row = ineligible
            .iter_mut()
            .find(|row| row_value(row, "accession").ok() == Some(target))
            .expect("target row should exist");

        let before = row.clone();
        row.insert("status".into(), "novel lineage".into());

        let fired = rule
            .apply(&data, &mut ineligible, &rec, &training, Pool::All, 99)
            .expect("N6 application to an ineligible status should succeed");

        assert_eq!(fired, 9, "only the other nine eligible rows should fire");

        let after = ineligible
            .iter()
            .find(|row| row_value(row, "accession").ok() == Some(target))
            .expect("target row should remain present");

        assert_eq!(
            row_value(after, "status").unwrap(),
            "novel lineage",
            "N6 must not demote an ineligible status"
        );
        assert_eq!(
            row_value(after, "its_rule_fired").unwrap(),
            "0",
            "N6 should record that the rule did not fire"
        );

        for field in ["placement", "novelty_claim", "reasons"] {
            assert_eq!(
                after.get(field),
                before.get(field),
                "N6 unexpectedly changed {field} for an ineligible status"
            );
        }

        // Also exercise the application path with main-type targets only.
        let mut main_adjusted = baseline.clone();
        let main_fired = rule
            .apply(&data, &mut main_adjusted, &rec, &training, Pool::Main, 99)
            .expect("N6 main-pool application should succeed");

        assert_eq!(main_fired, 10);
        for row in &main_adjusted {
            assert_eq!(row.get("its_pool").map(String::as_str), Some("main"));
            assert_eq!(row.get("its_rule_fired").map(String::as_str), Some("1"));
        }
    }
}
