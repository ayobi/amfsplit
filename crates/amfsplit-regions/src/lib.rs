//! Region extraction for long-read SSU–ITS–LSU sequences of arbuscular mycorrhizal fungi.
//!
//! This crate is the Rust port of the Stage 1 prototype (`amfsplit.py` v0.3.1, in the
//! amf-novelty-benchmark repository): ITS sub-regions and flanks come from ITSxRust's per-read
//! anchors, and the legacy barcode windows (MaarjAM NS31–AM1, LSU LR0R–FLR2 and others) are cut
//! at their original primer sites, so each window stays comparable with the database built on
//! it. Outputs are byte-identical to the prototype's; see `parity/stage1` in the workspace.

pub mod anchors;
pub mod concat;
pub mod panel;
pub mod pycompat;
pub mod search;
pub mod seqio;
pub mod stage1;
pub mod triage;
pub mod window;

pub use stage1::{Params, RunReport, STAGE1_ALGORITHM_VERSION, run};
