# AMFSplit novelty-calling specification v0.1

**Status:** design draft for review; not production behaviour and not a validated biological classifier.  
**Scope:** linked long-read AMF assignments using MaarjAM VT (SSU), UNITE SH (ITS), and LSU phylogenetic evidence.  
**Non-goal:** establishing or formally naming new fungal taxa. AMFSplit can nominate candidates for expert/taxonomic follow-up; sequence evidence alone does not establish a new species.

## 1. Why this specification exists

AMFSplit now has a frozen reproducibility baseline. The Rust M2b naming port passed byte-identical parity for the physical-culture and species-holdout golden outputs, and the independently rebuilt LSU placement workflow reproduced all pinned reference artifacts plus every cell of all 903 query rows. Those results establish reproducibility of the implementation and placement pipeline, not the biological validity of every novelty interpretation.

The next stage must therefore separate the observed evidence, the taxonomic assignment, quality/divergence warnings, and any novelty hypothesis. This prevents a missing database hit, a low-resolution assignment, or a QC warning from being treated as proof of a new lineage.

## 2. Source audit: what the current rules do

The current rule implementation is in `parity/naming/frozen_sources/m2b/novelty_baseline_consistent.py`; the fixed decision block used in the placement comparison is also represented in `decision_tree.py`. Do not silently rewrite those frozen files. New behavior should be implemented as a versioned experimental layer and compared against this baseline.

### 2.1 Reconciliation evidence states

For each arm (VT, SH, LSU), reconciliation distinguishes:

- **matched:** the query's evidence units connect to one or more *training* cultures and therefore have taxonomic holders in that fold;
- **unmatched:** the arm yielded unit(s), but none is supported by another training culture;
- **no call:** the query has no usable unit for that arm under the frozen arm rules.

These are not interchangeable. `unmatched` is evidence of a reference-coverage gap under the defined unit system; `no call` may reflect thresholding, coverage, marker absence, filtering, or other limitations. A `no call` does not establish novelty.

The current reconciliation intersects genera across matched arms and intersects species across matched arms that supply named-species evidence. It reports `species`, `species group`, `genus`, `discordant`, or `unplaced`. Only training culture holders should contribute to a holdout's names.

### 2.2 Existing Phase 4 precedence and reason codes

The frozen Phase 4 rules test a series of operational screens before the final reconciliation/novelty branch:

1. **A1 artifact screen:** compare SSU-flank distance against the nearest main-type copy/copies from other eligible cultures and compare against full-ITS distance. Defaults are SSU distance >= 1.0% and SSU distance >= 0.5 times the ITS distance; all SSU-tied nearest copies must agree. This is a heuristic, not a diagnostic test. Output uses `artifact_claim=screen_flag_only`.
2. **D1 short-5.8S divergent class:** the 5.8S segment is at least 3 bp shorter than the training/reference baseline and at least 3% different under the frozen comparison rule.
3. **D2 paralog-unit screen:** a VT/SH anchor is reached in at least two other cultures by short-5.8S copies but by no main-type copy. This remains a sequence-class screen, not a claim about genomic origin.
4. **D3/D4 and VT-boundary flags:** same-length divergence or a shorter 5.8S below the distance cutoff is flagged for review; it is not automatically called novel. With VT minimum identity 97% and a one-percentage-point boundary, a no-call with best VT identity in [96%, 97%) gets a boundary status. Threshold-adjacent or contradictory cases need coverage/provenance review.
5. **Reconciliation to status:** species/species-group calls are labelled known unless VT or LSU is unmatched; those unmatched cases are demoted to “novel species, genus known.” A genus-level call becomes “novel species, genus known.” Discordant matched arms remain discordant. An unplaced result becomes “novel lineage.” Novel labels carry `novelty_claim=candidate_only`.

The rule tree is deterministic, but the words “novel species” and “novel lineage” are operational candidate labels, not taxonomic conclusions.

## 3. Experimental LSU unit evidence

The independent M2b report describes two prespecified LSU units: (a) distinct credible EPA-ng edge IDs, and (b) credible edges collapsed within maximal pure V18 species-label clades, preserving separate components for non-monophyletic labels and leaving mixed-species edges separate. The credible edge set is the smallest set of highest-LWR edges covering >=0.95 total weight, including cutoff ties.

The archived M2b development report names species-clade grouping as the leading *experimental* option. In its training-only comparison:

| Measure | Frozen baseline | LSU species-clade experiment |
|---|---:|---:|
| Represented-species main copies with compatible known call | 748/755 (99.1%) | 751/755 (99.5%) |
| Exact-species calls among those 755 | 603 (79.9%) | 701 (92.8%) |
| Group calls containing the true species | 145 | 50 |
| Wrong final known names among represented-species copies | 0 | 0 |
| False-known calls among seven single-culture species | 3/48 (6.3%) | 0/48 |
| False-known calls among all named held-out-species main copies | 215/803 (26.8%) | 71/803 (8.8%) |
| Species-weighted false-known rate | 23.5% | 8.5% |
| Culture-weighted false-known rate | 24.7% | 7.3% |

These are development-cohort results, not independent validation. The species-holdout removes species from the culture bridge but **does not remove the target species from the external V18 placement reference**. It is not a full external-taxon removal test.

Important trade-off: 146 copies changed from “novel species, genus known” to “novel lineage” under the LSU species-clade experiment. That is loss of culture-bridge genus resolution, not evidence for genuinely deeper novel lineages. The experiment also has three known-species demotions, while recovering six earlier demotions. Keep this variant versioned and experimental until it passes additional validation.

LWR is relative attachment support conditional on the fitted reference alignment/model. It is not the probability that a taxonomic name is correct. The report notes false-known pairs even where best-edge LWR is >=0.95; do not use a high LWR as a substitute for cross-marker or taxonomic validation.

## 4. Proposed output model: keep separate concepts separate

For an experimental v0.1 interface, add orthogonal fields instead of making one overloaded `status` carry assignment, novelty, and QC simultaneously. Retain legacy columns during the comparison period.

| Field | Suggested values | Meaning |
|---|---|---|
| `assignment_rank` | `species`, `species_group`, `genus`, `discordant`, `unresolved` | What the available matched-arm evidence supports taxonomically |
| `assignment_name` | name/list or empty | Taxon name(s) justified by the assignment rank |
| `novelty_state` | `known_supported`, `novel_species_candidate`, `novel_lineage_candidate`, `unresolved`, `not_assessed` | Hypothesis about reference novelty, with candidates explicitly provisional |
| `quality_flags` | zero or more reason codes | Artifact screen, short 5.8S signal, paralog-unit signal, VT boundary, insufficient coverage, marker conflict, etc. |
| `arm_evidence` | per-arm state plus units/holder counts | Auditable evidence for each marker |
| `decision_rule` | stable rule ID(s) | Explains which rule produced the output |
| `algorithm_version` | immutable version tag | Enables comparing variants without overwriting the frozen baseline |

The existing `status`, `novelty_claim`, and `artifact_claim` fields remain as compatibility fields until downstream callers have migrated. `artifact_screen` must never itself imply `novel_lineage_candidate`.

## 5. Decision principles for v0.1

These are the principles the implementation and tests should enforce. Exact thresholds and any newly proposed combinations are hypotheses to be evaluated, not assumed facts.

1. **Evidence before label.** Report per-arm `matched`/`unmatched`/`no call`, alignment eligibility/retention, and key placement measurements before deriving a novelty category.
2. **Missing evidence is not novelty.** If a result is unplaced because all arms are `no call`, or evidence is below QC/retention requirements, prefer `unresolved`/`not_assessed`; do not equate it with a candidate novel lineage. Evaluate how this refinement changes historical outputs.
3. **An unmatched informative unit is a candidate signal, not proof.** It indicates no supporting holder in the training culture bridge under that arm's unit definition. It may reflect incomplete sampling, database structure, a known but unrepresented species, or a genuinely divergent lineage.
4. **Do not allow one low-confidence marker to overrule all others invisibly.** Any VT/LSU unmatched veto of a species/group assignment must be surfaced as a named, testable rule and compared against alternatives, including SH-only mismatch and marker-specific no-calls.
5. **Cross-marker disagreement is a separate outcome.** A matched-arm genus intersection that is empty should report `discordant`; do not convert it straight into novelty. If only one arm matches, preserve that fact and use `unresolved` unless the predeclared evidence rule supports a candidate status.
6. **Known group versus novel species.** A `species group` call containing the correct species is successful coarse resolution on the negative control. When the true species is held out, any species or species-group call is a false-known outcome, even if the named group is plausible.
7. **Sequence type/QC flags are orthogonal.** Short-5.8S and possible artifact flags do not determine taxonomy without their own validated criteria. A same-length divergent marker should remain a warning unless supported by further evidence.
8. **LSU clade geometry informs compatibility, not a name by itself.** Training culture holders provide taxonomic names. External reference clade labels are useful for unit construction and audits but should not be substituted directly for a validated culture prediction.
9. **No tuning on the reported holdout outcomes.** Treat current holdouts as development benchmarks; tune only inside training folds or on separate development splits. Retain an untouched validation set for final evaluation.
10. **No taxonomic claim.** The tool reports “candidate” or “unresolved” outputs requiring taxonomic review; it does not name new species.

## 6. Proposed provisional decision table

| Evidence case | Assignment output | Novelty output for v0.1 | Reason/guardrail |
|---|---|---|---|
| Species call; VT and LSU not unmatched; sufficient quality | `species` | `known_supported` (provisional) | Known compatibility on the culture bridge; still subject to holdout performance |
| Species-group call; VT and LSU not unmatched; sufficient quality | `species_group` | `known_supported` at group level | Preserve uncertainty; do not render as an exact species |
| Species/group call but informative VT or LSU evidence is unmatched | Keep reconciled assignment in evidence fields; downgrade final taxonomic confidence | `novel_species_candidate` only as a tested rule variant | Existing N1/N2 behaviour; measure known-species sensitivity loss and false-known reduction |
| Genus call; one or more matched arms support a genus but no species survives | `genus` | `novel_species_candidate` only as a tested rule variant | Existing N3; broad genus support does not alone prove a novel species |
| Matched arms have no common genus | `discordant` | `unresolved` | X1; investigate marker conflict before novelty interpretation |
| No matched arms, but one or more adequately observed arms are explicitly unmatched | `unresolved` taxonomically | `novel_lineage_candidate` only if a prespecified candidate rule is met | Separate explicit reference mismatch from absent/no-call evidence; require audit trail |
| All arms no-call, low-quality/non-retained sequence, or insufficient marker coverage | `unresolved` | `unresolved` / `not_assessed` | Never infer novelty from failure to obtain usable evidence |
| D1/D2 divergent-type screen or A1 artifact screen | Assignment evaluated separately where possible | Novelty not asserted by the QC flag itself | Preserve QC flag and rule code; require additional evidence |
| VT threshold-boundary signal | `unresolved` pending coverage/provenance check | `unresolved` | Do not force a known/novel category just below a threshold |

This table is intentionally conservative and exposes cases where the current code is more aggressive than the proposed interpretation. The table is a design target for a *new experimental version*, not an instruction to edit the frozen implementation in place.

## 7. Evaluation design

### 7.1 Maintain separate tests

- **Gate 0 regression:** Rust output remains byte-identical to its frozen Python baseline. No novelty experiment may change this test.
- **Physical-culture holdout:** exclude the held culture and any co-cultured lineages under the benchmark protocol. Tests recognition of represented species without culture self-leakage.
- **Species holdout:** exclude all cultures of the target species and co-cultured lineages from culture holders. In the current LSU experiment, external V18 placements still contain their external sequence/taxonomy, so this tests the culture bridge, not external-reference absence.
- **True external-taxon holdout (future):** remove the target species from external placement/reference inputs as feasible; freeze the removal manifest before evaluation. This is necessary for stronger claims about detection of taxa absent from the reference framework.
- **Environmental validation (future):** use independent samples/replicates and report candidate clusters, not ground-truth novelty, unless independent taxonomy supports it.

### 7.2 Required metrics

Report row/copy-level counts and macro-averages by species and culture. Do not treat copies from one culture as independent biological replicates.

- False-known: known species or known species-group call for a held-out species.
- Known compatible retention: represented-species query still receives a known assignment that contains the truth.
- Exact species accuracy and containing-group rate among represented species.
- Wrong known name rate.
- Genus-resolution retention (to detect the observed loss of genus calls).
- Novelty-candidate precision/recall where the holdout label actually justifies that interpretation; otherwise report false-known as the cleaner measured endpoint.
- Unresolved/no-call rate and coverage/retention-stratified results.
- Short-5.8S recognition, artifact screen count, and false novelty claims on known divergent types.
- Paired per-accession changes from baseline, with reason code and decision rule.

### 7.3 Acceptance criteria before promotion

Set numerical acceptance thresholds before running any new tuning. At minimum, a candidate variant must:

- preserve Gate 0 byte-identical output;
- demonstrate no unreported known-species regression among the represented-species controls;
- reduce false-known calls without hiding errors by turning everything into unresolved/no-call;
- report the sensitivity/resolution cost and per-species/culture macro metrics;
- correctly handle no-call vs unmatched vs discordant states in targeted tests;
- preserve training-only holder construction and pass leakage assertions;
- be reproducible from pinned source, inputs, parameters, and environment metadata.

If a variant improves false-known performance but sacrifices genus resolution or known-species retention, report that trade-off rather than silently choosing it. Final promotion needs an independent validation design; development holdouts alone are insufficient.

## 8. Implementation plan

1. **Audit and snapshot:** preserve the frozen `novelty_baseline_consistent.py`, `decision_tree.py`, M2b report, baseline rows, and placement-clade comparison. Record their checksums and the exact benchmark repository commit.
2. **Reproduce the development comparison locally:** run the baseline parity comparison and species-clade variant on the pinned 903-copy cohort; check current local files/results before treating the external report as a newly rerun result.
3. **Build a case audit table:** all 903 copies with per-arm state, reconcile rank, legacy status/reason, species/culture holdout class, LSU compatible units, quality metrics, and new proposed output axes.
4. **Freeze decision scenarios:** explicitly enumerate every combination of per-arm state (`matched`, `unmatched`, `no call`), rank, QC flag, and VT/LSU veto. Unit-test precedence and contradictory cases.
5. **Run ablations:** compare (a) legacy Phase 4, (b) orthogonal-output-only with unchanged decisions, (c) no-call/unmatched conservative distinction, (d) current VT/LSU unmatched veto, (e) alternatives to the veto, and (f) species-clade versus edge LSU units. No threshold fitting on the evaluation folds.
6. **Review difficult cases:** include the six recovered known-species copies, three new species-clade demotions, remaining six false-known sister-species pairs, 146 genus-to-lineage resolution losses, all `unplaced`/no-call cases, and the flagged artifact/divergent examples.
7. **Add independent external-taxon validation:** specify exactly what is removed from the LSU reference and marker databases, then lock and run the test without using outcomes to tune rules.
8. **Implement separately in Rust:** only after the decision table and tests are accepted; keep Gate 0 frozen and add a separate versioned novelty test suite.

## 9. Reproducibility and provenance

The audited M2b report and comparison outputs were read from `ayobi/amf-novelty-benchmark` at commit `02c0844604f958d6326645829e20e32ca0996b03` (commit date 2026-10-16). Record this commit explicitly when reproducing the report; do not rely on a moving `main` branch. Keep the existing Gate 0 baseline and LSU reconstruction as separate checks. Also record binary hashes and versions for MAFFT, RAxML-NG, and EPA-ng: an exact output match establishes the artifacts were reconstructed in the successful run, but does not by itself prove every installed binary hash matches the original provenance lock.

**Runner-layout caveat:** the copy of `compare_holdouts.py` in `ayobi/amfsplit/parity/naming/frozen_sources/m2b/` imports from adjacent `code/baseline` and expects `inputs/previous_results/`, which are not part of the standalone AMFSplit tree. Do not assume that frozen source copy can run in isolation. Use the complete pinned `amf-novelty-benchmark` checkout for the development comparison, or create a separately packaged runner with all inputs and imports made explicit; do not silently fill missing files from a different version.

## 10. References within the project

- `parity/naming/frozen_sources/m2b/PROTOCOL.md`
- `parity/naming/frozen_sources/m2b/reconcile.py`
- `parity/naming/frozen_sources/m2b/novelty_baseline_consistent.py`
- `parity/naming/frozen_sources/m2b/decision_tree.py`
- `parity/naming/frozen_sources/m2b/compare_holdouts.py`
- `parity/naming/README.md`
- `ayobi/amf-novelty-benchmark`: `benchmarks/m4_external_holdout/frozen_m2b/AMF_M2b_Report.md`, protocol, per-copy comparison outputs, and software provenance

**Decision requested before coding:** accept or amend the proposed separation between (1) taxonomic assignment, (2) novelty state, and (3) independent QC/divergence flags; and decide whether the no-call/unmatched distinction should become a behaviour-changing experimental variant in the next benchmark run.
