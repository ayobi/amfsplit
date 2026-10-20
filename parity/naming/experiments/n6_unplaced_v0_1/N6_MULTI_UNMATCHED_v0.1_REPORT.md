# N6 multi-unmatched evidence gate: exploratory ablation v0.1

**Date:** 20 October 2026  
**Status:** Development-only experiment; not production behaviour and not a validated biological novelty classifier.  
**Frozen baseline:** M2b Phase 4 reconciliation and novelty outputs; left unchanged.

## Question

The frozen decision tree assigns `novel lineage` when reconciliation rank is `unplaced`, unless a higher-priority artifact, divergent-class, or boundary rule applies. An `unmatched` arm has a non-empty evidence/unit result with no holder among the remaining training cultures; a `no call` arm has no usable assigned unit. These states are not interchangeable, and neither alone proves biological novelty.

N6 tests a narrow guardrail for existing `novel lineage` + `unplaced` rows only:

- retain a `candidate_only` novel-lineage status when no arm is matched and at least two arms are `unmatched`;
- return `insufficient evidence` and `not_asserted` when fewer than two arms are unmatched or an unexpected matched-arm conflict exists;
- preserve every other baseline status and claim.

The rule is identified as `N6_MULTI_UNMATCHED_v0.1`. “Two unmatched marker arms” means multi-arm evidence against the culture bridge, **not statistically independent evidence** and not proof of a taxon absent from the external databases. SSU, ITS, and LSU regions belong to the same linked ribosomal marker system; marker/database-arm agreement must not be overstated as independent-locus replication.

## Input audit

The audit joins the frozen per-copy reconciliation and `baseline_consistent` tables for `physical_culture` and `physical_species`, with 903 unique accessions per mode. It is based on the local `amf-novelty-benchmark` checkout that was at commit `62862a0` when the audit was run. Source tables were read only.

### Existing unplaced calls

| Measure | Physical-culture holdout | Species holdout |
|---|---:|---:|
| Rows | 903 | 903 |
| Unplaced rows | 28 | 234 |
| Existing `novel lineage` among unplaced | 28 | 234 |
| Unplaced rows with all three arms `no call` | 0 | 0 |

Observed VT / SH / LSU patterns among those unplaced rows:

| Pattern | Physical-culture | Species holdout |
|---|---:|---:|
| `unmatched / no call / unmatched` | 20 | 70 |
| `unmatched / unmatched / no call` | 4 | 40 |
| `unmatched / unmatched / unmatched` | 4 | 121 |
| `unmatched / no call / no call` | 0 | 3 |

All unplaced rows in both modes have VT `unmatched`; none has a matched arm. The cohort does not contain an unplaced row with all three arms `no call`, so that edge case is covered by a unit test, not represented as an observed biological case.

## N6 result

Seven unit tests passed. The audit table was evaluated without modifying the frozen source tables or source implementation.

| Outcome | Physical-culture holdout | Species holdout |
|---|---:|---:|
| Copies evaluated | 903 | 903 |
| Baseline `candidate_only` claims | 67 | 617 |
| Experimental `candidate_only` claims | 67 | 614 |
| Changed copies | 0 | 3 |

The three changed copies are all species-holdout cases with one unmatched arm and two no-call arms:

| Accession | Culture | Species | Baseline | N6 v0.1 |
|---|---|---|---|---|
| `PX214751.1` | `4392_1156C` | *Dominikia duoreactiva* | Novel lineage, candidate only | Insufficient evidence, not asserted |
| `PX214765.1` | `3647_963E` | *Microdominikia litorea* | Novel lineage, candidate only | Insufficient evidence, not asserted |
| `PX214767.1` | `3647_963E` | *Microdominikia litorea* | Novel lineage, candidate only | Insufficient evidence, not asserted |

These are three sequence copies from two cultures and two named species, not three independent biological samples. The physical-culture outputs do not change. All baseline decisions outside the exact scope of N6 are preserved.

## Interpretation and limitations

This is a targeted ablation on development outputs that were already inspected when formulating the rule. It is **not independent validation** and does not show that the revised rule improves general taxonomic accuracy. The species holdout excludes target species from the culture bridge but can leave the target represented in the external V18 phylogenetic reference. Therefore it measures culture-bridge behaviour, not novelty detection under complete external-taxon removal.

The outcome `novel lineage` must remain candidate-only. N6 changes only a narrow evidence gate; it does not prove that two unmatched marker arms are sufficient to call a novel species or lineage. The next evaluation should compare sensitivity and abstention across all outcomes, and use a true external-taxon-removal holdout plus independent datasets before stronger biological claims.

## Reproduction

Run from the `amf-novelty-benchmark` checkout, with an output directory outside the checkout:

```bash
python3 audit_novelty_calls.py \
  --benchmark-root "$HOME/projects/amf-novelty-benchmark" \
  --out "$HOME/Downloads/amfsplit_novelty_audit_REPLACE_WITH_UNIQUE_SUFFIX"

python3 -m unittest -v test_unplaced_novelty_rule.py

python3 evaluate_unplaced_novelty_rule.py \
  --audit "$HOME/Downloads/amfsplit_novelty_audit_REPLACE_WITH_UNIQUE_SUFFIX/novelty_decision_audit.tsv" \
  --out "$HOME/Downloads/amfsplit_novelty_rule_n6_REPLACE_WITH_UNIQUE_SUFFIX"
```

The script/test files are for the repository's exploratory experiment folder. Do not replace the frozen files in `parity/naming/frozen_sources/m2b/`.
