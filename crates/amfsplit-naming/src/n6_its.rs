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
