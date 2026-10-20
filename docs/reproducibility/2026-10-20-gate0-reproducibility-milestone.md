# AMFSplit Reproducibility Milestone Report and Roadmap

**Milestone:** Gate 0 naming parity + independent LSU placement reconstruction  
**Status date:** 20 October 2026  
**Overall status:** **PASS — freeze recommended**  
**Project checkout:** `~/projects/amfsplit`  
**Last known implementation commit:** `a827275` — `Add Rust frozen M2b naming parity gate`

---

## 1. Executive summary

AMFSplit has reached an important reproducibility checkpoint. Two different parts of the workflow have now been independently checked against frozen outputs:

1. **Rust naming parity (Gate 0):** the Rust implementation reproduced the archived Python M2b per-copy naming outputs byte-for-byte for both the physical-culture holdout and physical-species holdout benchmarks. All rows and all compared cells matched; no differing cells were reported.
2. **Independent LSU placement reconstruction:** the external LSU phylogenetic-placement workflow was rerun using the frozen Delavaux V18 reference inputs and an archive-verified RAxML-NG executable. The fitted tree, fitted model, query alignment, edge annotations, and all 903 per-query placement rows exactly matched the frozen artifacts.

The Rust workspace release test suite also passed: **34 tests passed**. These results support freezing the current baseline before introducing additional novelty-detection behaviour.

The next work should preserve this baseline as an immutable regression target, bring the verification evidence and helper scripts into version control, and only then progress into the next scientific/algorithmic stage. Exact reproduction establishes that the current implementation and archived workflow can be reproduced under the recorded conditions; it does **not** by itself establish that every biological classification or novelty claim is correct.

## 2. Goal and reasoning

### 2.1 Why Gate 0 came first

The immediate engineering risk was that changes to the implementation could be confused with changes in the scientific output. Gate 0 therefore asked a narrow question: can the Rust naming implementation reproduce the already-frozen Python baseline exactly?

Byte-for-byte comparison is useful here because it avoids choosing arbitrary tolerances for identifiers, classifications, missing-value handling, ordering, or serialized values. Any difference becomes a specific regression to explain instead of an ambiguous change in aggregate accuracy.

### 2.2 Why LSU placement was checked independently

LSU placement depends on external reference data and several bioinformatics tools, not only AMFSplit's Rust code. A naming benchmark can appear to differ even when the naming logic has not changed if the upstream alignment, fitted tree/model, or phylogenetic placement changes.

We therefore rebuilt the frozen placement workflow independently and compared its intermediate and final artefacts with the frozen outputs. This separates toolchain/reproducibility problems from naming-logic problems and provides a traceable upstream input for later experiments.

### 2.3 Why we are freezing now

With both the naming output and the upstream LSU placements reproducing exactly, there is now a stable reference point against which future code changes can be assessed. Novelty rules should be developed as an explicitly testable layer on top of this baseline, not mixed into the parity implementation in a way that makes it harder to tell what changed.

## 3. Workstream A — Rust frozen naming parity (Gate 0)

### 3.1 Frozen reference set

The frozen culture reference fingerprint reported by the gate is:

`181938d695f3df40`

The retained reference set contains **903 sequence copies**, representing **151 cultures** and **37 named species**. The project notes also distinguish **852 main-type copies** within this set.

### 3.2 Exact parity results

| Benchmark | Holdout rounds | Rows compared | Result | Generated-output SHA-256 |
|---|---:|---:|---|---|
| `physical_culture` | 148 | 903/903 | Byte-identical; 0 differing cells | `6cc5412acf147ba7644d767851872517f78551b372e614c3f07297227d539f34` |
| `physical_species` | 50 | 903/903 | Byte-identical; 0 differing cells | `906420d70e9d8d24ab5b115ecd63e72451c396b4f1e7c60cfccf7b8d6a8a0800` |

The gate concluded:

```text
GATE 0: PASS (exact archived per-copy output parity)
```

This is stronger than agreement on a small set of aggregate metrics: the archived per-copy outputs and every compared cell matched for both holdout designs.

### 3.3 Rust test-suite result

Command used:

```bash
cargo test --workspace --release
```

**Result: all 34 tests passed.** The suite included the new naming parity coverage and existing stage-1, triage, reference, Edlib, and region-related tests. Some non-fatal compiler warnings were observed (unused imports, variables, or fields); they did not cause test failures and can be addressed separately from the frozen scientific baseline.

### 3.4 Source milestone

The implementation commit recorded in the project history is:

```text
a827275 Add Rust frozen M2b naming parity gate
```

At the last recorded check, the source changes were committed, while `AMFSplit_Gate0_Reproduction_and_Wiring_Report.md` and `amfsplit_gate0_wiring.patch` were present as untracked files. Re-check `git status --short` before finalizing the documentation commit. Inspect the patch and report rather than automatically committing or deleting either one.

## 4. Workstream B — Independent LSU placement reconstruction

### 4.1 Frozen inputs and reference provenance

The rebuild used the Delavaux V18 reference checkout at:

`~/Downloads/AMF_test/ref/delavaux/repo/2024_AMFPipeline/2024_AMFPipeline_ASV`

The verified V18 repository revision was:

`f4a0014336f49be6aea48e8374dae052360a88bd`

The frozen placement script was SHA-256 checked by the verification step. The query input was the project's frozen `LSU_LROR-FLR2.fasta` file. The verifier also confirmed the reference input files and the query-alignment files against the frozen checksums.

### 4.2 Rebuild procedure

The frozen script performed the following stages:

1. **MAFFT:** construct the reference-only alignment.
2. **RAxML-NG:** evaluate the fixed backbone topology, estimate the model/tree parameters under `GTR+G4`, and write the fitted reference tree and model.
3. **MAFFT:** add the LSU query fragments to the reference alignment while retaining the reference alignment length.
4. **EPA-ng:** place query sequences on the fitted reference tree and generate the placement/edge summary, using credible edge mass `0.95`.
5. **Independent verification:** compare the recreated artefacts and per-query table with the frozen benchmark.

The completed rebuild reported **903 queries**, **495 reference tips**, **987 tree edges**, **0.95 credible mass**, and **903/903 alignment-eligible queries**.

### 4.3 Toolchain and executable evidence

| Component | Version / identity | Evidence recorded |
|---|---|---|
| MAFFT | 7.526 | Version checked; reference/query alignments reproduced exactly |
| RAxML-NG | 1.2.2 (`1.2.2-master`, released 30 April 2024) | Pinned executable SHA-256 matched exactly |
| EPA-ng | 0.3.8 | Version checked; placement artefacts reproduced exactly |
| Python | 3.12.14 | Runtime used for the rebuild in `amf-m4` |
| Biopython | 1.88 | Import verified |
| Edlib | 1.3.9.post1 | Import verified |

Pinned RAxML-NG archive SHA-256:

`962be94d9b7469757fa27af05c32fad6c7605408e8975b13aa58e6ea91a48d81`

Pinned RAxML-NG executable SHA-256:

`b5038949eb005d14c48efe47f1ddb4fa750e89dbafeb30c46e05c651bacba6d5`

The executable was located at:

`~/Downloads/amfsplit_pinned_tools/raxml-ng-1.2.2/raxml-ng`

**Provenance caveat:** earlier checks found that the installed MAFFT launcher and EPA-ng executable hashes did not match the hashes recorded in the frozen software-provenance manifest, although their reported versions were as expected. The verification run demonstrates exact output reproduction for this particular run, but the precise origin of those two executables should be reconciled before making a stronger claim that every executable byte matches the provenance manifest. The RAxML-NG executable, by contrast, was independently checksum-verified against the expected pinned binary.

### 4.4 Exact artefact checks

The independent verifier reported `PASS` for each of the following:

| Artefact | SHA-256 |
|---|---|
| `reference/reference.unaligned.fasta` | `874308b806b60ddb4ddf917091893a2138c3538bda304d95703ae3b139c48d66` |
| `reference/reference.aligned.fasta` | `1cd6426557f842d9e0085b95f02f3156406d43fb71f840a2c6891bece769ce33` |
| `reference/reference_metadata.tsv` | `5c38fd7fddd3896492644503359917d5fbd1ac3d4c818c0b4c7fc8a49ff5591d` |
| `reference/backbone.topology.newick` | `a6ae1cb91d1ca04021028d066a4b30f2a77c1235d13cfcfdafa9598eb5f0fa66` |
| `reference/fit.raxml.bestTree` | `d3da113b83e3d185ebfe10fd1958aa121e401c0dd035c74350ef607a74abf11b` |
| `reference/fit.raxml.bestModel` | `2e7a4dd93fd75c60586bad52a2fd4e4c5a73043289e89ad87f686df21e2a6194` |
| `queries.fasta` | `18679043a974d1bff4fdf798c0f74655d7444e1b0916a5c0c8356c796a9cca2f` |
| `queries.aligned.fasta` | `c8f5d1e41e049aad46e9197396cfe24ec6accaa61deafa38a56c930fe739755f` |
| `summary/edge_annotations.json` | `d54eca9e489977de8c6dead886c48cf83028a4e61cdb67a766e7a7b67bd25dd6` |
| `summary/lsu_placements.copies.tsv` | `1de103b6ae45f0372dda9ab0799230c4d1243cb98add943e5bab67785e5e910d` |

The frozen placement script SHA-256 recorded in the benchmark provenance is:

`359f297a4c7038327fc94edea19ec4b8dab9f91a723b64b7c5a09ca4810610cc`

The final verifier output stated:

```text
PASS summary.queries: 903
PASS summary.reference_tips: 495
PASS summary.tree_edges: 987
PASS summary.credible_mass: 0.95
PASS summary.alignment_eligible: 903/903
PASS Per-cell comparison: exact (all 903 rows and every column match)

PLACEMENT REBUILD: PASS (all pinned artifacts and all 903 per-query summary rows match exactly)
```

### 4.5 What the first failed comparison taught us

An earlier rebuild using the environment-installed RAxML-NG binary reproduced the input alignments but did **not** reproduce the fitted tree/model or final placement output. It showed thousands of per-cell differences and a best-edge difference for query `PX207259.1` (edge 182 in that run versus edge 181 in the frozen output).

After replacing that executable with the archive-verified RAxML-NG binary, the fitted tree, fitted model, edge annotations, and full placement summary all matched exactly. Since the reference/query inputs and the other recorded tool versions were retained, this is strong evidence that the previous RAxML executable/build was the source of the earlier discrepancy. It is not a general proof that tool versions alone uniquely determine a result; the executable provenance caveat above remains worth addressing.

### 4.6 Preserved run locations

Successful rebuild directory:

`~/Downloads/amfsplit_lsu_pinned_raxml_20261020_retry1`

Build log:

`~/Downloads/amfsplit_lsu_pinned_raxml_20261020_retry1.log`

Verification helper currently used:

`~/Downloads/amfsplit_lsu_check/verify_lsu_placement_rebuild.py`

These paths refer to the user's local machine. They are not automatically included in the downloadable report or in the Git repository.

### 4.7 Operational notes during the rebuild

The first background launch was inadvertently run under the base Conda environment and stopped immediately with `ModuleNotFoundError: No module named 'Bio'`. The environment was corrected by activating `amf-m4`; Biopython 1.88 and Edlib then imported successfully, and the subsequent detached rebuild completed.

During debugging, running `ps -fp 557170` was observed to close the terminal window unexpectedly. The cause was not established and is not required to explain the successful placement run. The successful run was launched with `nohup` and its output was preserved in the log listed above. Treat the terminal behaviour as a separate host/environment issue to investigate if it recurs, not as evidence of a placement error.

## 5. Interpretation: what has and has not been established

### Established by the recorded checks

- The Rust naming implementation exactly reproduces the archived naming outputs for the two recorded holdout designs.
- The workspace release tests passed at the reported checkpoint.
- The frozen LSU placement workflow can be reconstructed on this host with exact matches across reference inputs, fitted tree/model, query alignment, edge annotations, and all 903 per-query rows.
- A non-identical RAxML executable was associated with a differing fitted tree/model and placement summary; using the checksum-verified executable restored exact matches.

### Not established by these checks alone

- That the biological taxonomy or all species-level assignments are correct in every dataset.
- That novelty detection generalizes to unseen species, cultures, geographic regions, or sequencing protocols.
- That the novelty threshold or marker-consensus policy is optimally calibrated.
- That MAFFT and EPA-ng executable hashes have been reconciled with every value in the frozen software-provenance manifest.
- That the successful local placement run is fully captured as a reproducible, repository-managed regression procedure. The verification helper and raw run output are still outside the project checkout at the recorded paths.

These are follow-on work items, not failures of the two exact-reproduction checks.

## 6. Freeze checklist — do this before further algorithm changes

### 6.1 Preserve the successful evidence

- Keep the successful run directory and log unchanged at the paths recorded above.
- Keep the frozen query FASTA, expected summaries, source revision, and software provenance intact.
- Treat the archived Rust naming outputs and LSU placement outputs as regression fixtures. Any intentional update should require an explicit reason and regenerated evidence.

### 6.2 Commit the report and inspect the working tree

After downloading this report, place it in the repository, for example at:

`docs/reproducibility/AMFSplit_Reproducibility_Milestone_Report.md`

Then, from `~/projects/amfsplit`, inspect before staging:

```bash
cd ~/projects/amfsplit
git status --short
git diff --check
```

Review the existing `AMFSplit_Gate0_Reproduction_and_Wiring_Report.md` and `amfsplit_gate0_wiring.patch` if they are still untracked. Do not add the patch blindly, and do not delete it until its purpose and whether it is redundant are clear. Stage the milestone report explicitly and make a separate documentation commit, for example:

```bash
mkdir -p docs/reproducibility
# Copy the downloaded Markdown report into docs/reproducibility/ first.
git add docs/reproducibility/AMFSplit_Reproducibility_Milestone_Report.md
git diff --cached --check
git diff --cached --stat
git commit -m "docs: freeze Gate 0 and LSU placement reproducibility milestone"
```

The commit command is a suggested local action; the report has not itself been committed to the user's local repository from this chat. After the documentation commit, check `git status --short` again. A tag can be considered only once the intended files are committed and the working tree is understood.

### 6.3 Make the evidence repeatable by other developers

Before treating this as a permanent project regression test, move or copy the useful preflight and verification helpers from `~/Downloads/amfsplit_lsu_check/` into a documented repository location (for example, `parity/naming/scripts/`). Review and test the scripts before committing them. Store small text manifests, commands, checksums, tool versions and expected summaries in version control; do not commit large generated intermediate files unless the repository's data policy intends that.

## 7. Proposed roadmap from here

The order below deliberately puts reproducibility infrastructure before new novelty behaviour.

### Phase 1 — Finish the freeze and make verification portable

**Goal:** a clean, reviewable record of what passed, with the evidence runnable again by someone else.

1. Commit this milestone report and resolve the previously untracked report/patch files intentionally.
2. Bring the LSU preflight/verification scripts into the repository, with clear inputs and non-zero exit codes on mismatch.
3. Record exact tool provenance for MAFFT and EPA-ng and reconcile their executable/package identities against `software_provenance.json`.
4. Add a short top-level reproducibility guide that distinguishes the fast Rust Gate 0 test from the longer external LSU rebuild.
5. Verify the clean checkout and commit state; optionally tag the milestone after review.

**Exit criteria:** a clean Git tree; a durable milestone commit; Gate 0 and the external verifier documented and repeatable; remaining provenance differences either resolved or explicitly recorded.

### Phase 2 — Define the naming and novelty contract before implementation changes

**Goal:** make the intended behaviour precise enough to test independently from code.

Write down the inputs, outputs, decision levels, and abstention behaviour for VT, SH, and LSU evidence. Define the distinction between an exact species call, group-level call, novel species candidate, novel lineage candidate, divergent/unresolved result, and no-call. Document how marker evidence is reconciled and which conflicts should result in abstention rather than an overconfident label.

Do not change frozen Gate 0 behaviour while writing this contract. It should describe the intended new layer and explicitly identify any intentional departures from the archived baseline.

**Exit criteria:** a reviewed decision table and a compact set of hand-checkable examples covering known matches, ambiguous matches, conflicting markers, missing marker evidence, and candidate novelty.

### Phase 3 — Implement novelty as a separate, testable layer

**Goal:** add novelty decisions without weakening baseline reproducibility.

1. Preserve the frozen naming path and its exact fixtures.
2. Put novelty scoring/rules and their thresholds behind a separately testable interface or mode.
3. Add unit tests for threshold boundaries, missing data, ties, marker conflicts, and abstention.
4. Keep novelty status explanations explicit so that a result can be traced to its evidence and rule path.
5. Version any thresholds or policy changes so result tables can identify which rule set produced them.

**Exit criteria:** Gate 0 stays exact and green; novelty tests pass independently; each result has an interpretable evidence trail.

### Phase 4 — Benchmark biological performance without leakage

**Goal:** estimate how well the novelty layer distinguishes known taxa from genuinely held-out taxa.

Use physical-culture holdouts and species-level holdouts as separate evaluation questions. Ensure that references, thresholds, or feature selection used for each fold do not leak information from the held-out culture/species. Compare marker-specific evidence and combined decisions. Report both coverage and correctness, not accuracy alone.

Recommended metrics include exact-species precision/recall, group-level accuracy, false species assignment rate, novel-candidate precision/recall where labels permit, abstention/no-call rate, and performance broken down by marker and holdout design. Use uncertainty intervals or repeated folds where appropriate. Do not tune thresholds on the same held-out predictions later used as the headline evaluation.

**Exit criteria:** documented fold construction and leakage checks; per-fold and aggregated tables; error analysis of false species calls, false novelty calls, and abstentions; a defensible threshold-selection policy.

### Phase 5 — External validation and robustness

**Goal:** test whether results hold beyond the exact frozen benchmark.

Evaluate representative external datasets and differing read quality/length conditions. Test sensitivity to incomplete marker evidence, reference imbalance, divergent sequences, and alternative but valid placement uncertainty. Clarify which conclusions are specific to the current reference collection and which transfer to new data. Keep external validations separate from the frozen parity tests.

**Exit criteria:** an external validation report, documented failure modes, and explicit limits on the intended scope of use.

### Phase 6 — Performance, interface, and publication readiness

**Goal:** turn the tested method into a reliable research tool and defensible scientific result.

After correctness is stable, measure runtime and memory on defined hardware, document CLI inputs/outputs and reproducible environment setup, add examples, and prepare release-quality documentation. For a manuscript, present the frozen-baseline parity as software validation, then present novelty performance and biological interpretation as separate evidence. Avoid using exact reproducibility alone as evidence of biological novelty detection accuracy.

**Exit criteria:** reproducible example run; clear method/version provenance; stable output schema; benchmark tables and limitations suitable for review.

## 8. Immediate next actions

The recommended next actions, in order, are:

1. **Commit this report** into `docs/reproducibility/` after checking the working tree and resolving the status of the existing untracked files.
2. **Promote the LSU verifier/preflight scripts into the repository** and make their prerequisites and expected outputs explicit.
3. **Resolve or document the MAFFT/EPA-ng executable provenance discrepancy** without changing the frozen successful run.
4. **Write the novelty decision contract** and sample decision table before editing novelty code.
5. **Implement novelty separately and evaluate with leakage-safe holdouts**, keeping Gate 0 as a permanent regression gate.

## 9. Milestone conclusion

As of 20 October 2026, AMFSplit has an exactly reproduced frozen Rust naming baseline and an independently reconstructed LSU placement baseline. This removes a major class of ambiguity from the next stage: future differences can be checked against stable archived outputs and examined as intentional algorithm changes rather than unexplained changes in upstream placement.

**Recommendation: freeze the evidence and repository state now, then proceed to the portable verification harness and a formally specified novelty-decision layer.**
