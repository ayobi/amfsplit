# M2b development protocol, specified before placements were inspected

- Freeze V18 repository commit f4a0014336f49be6aea48e8374dae052360a88bd; verify AMF-only digest against the uploaded evidence.
- Start with the 903 LROR–FLR2 windows. D2 is outside this first run.
- Align external references alone with MAFFT L-INS-i (localpair, 1000 refinement iterations, 8 threads with single-threaded iterative refinement for reproducibility). Preserve the published backbone topology, re-estimate branch lengths and GTR+G4 model parameters using only those references.
- Add each query as a fragment into the fixed reference alignment with MAFFT --addfragments --keeplength --mapout. Verify the reference alignment is unchanged. Record deleted query insertions and ambiguous bases; require >=90% query retention for the experimental taxonomy bridge. Placement itself is attempted for every nonempty aligned query.
- Use EPA-ng 0.3.8 with exhaustive edge evaluation (--no-heur), no premasking, and retain all placement weights. Likelihood weight ratios (LWRs) measure relative support conditional on this reference/alignment/model; they are not probabilities that a species name is correct.
- Primary compatible set: minimum set of highest-LWR edges containing >=0.95 mass, including every edge tied at the cutoff. No threshold is selected using held-out results.
- Report raw best edge, LWR, entropy, pendant length, credible edges, genus/species purity of supported reference clades, alignment retention and provisional taxonomy.
- Benchmark two explicitly experimental LSU unit definitions: (1) credible edge IDs; (2) credible edges collapsed within maximal pure V18 species-label clades. A polyphyletic reference species remains split. Mixed-species edges stay separate; they are not assigned the name of a nearest tip.
- Build each training culture's unit set using the existing intersection-if-nonempty, otherwise union rule over eligible main copies. Learn VT and SH merges/margins only from training cultures. D1/D2 and the Phase 4 rules remain fixed. No held-out culture names/classes contribute to training pools.
- Use physical-culture exclusions, including other lineages sharing the same culture label; named species holdout excludes all cultures of the species and their co-cultured lineages.
- Compare against the frozen training-only baseline. Evaluate false-known errors, known-species retention, exact/group/genus resolution, short-type recognition, and the six earlier demotions. Do not select or promote a method solely because it gives fewer false-known calls; loss of known-species sensitivity must be shown.
- Retain external databases unchanged. This is a culture-bridge development experiment, not a sealed test or full external-taxon removal experiment.
