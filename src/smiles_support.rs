//! Public `smiles-parser` adapter surface for RDKit-parity fingerprinting.

pub use crate::smiles_support_impl::{
    SmilesPreparationError, SmilesRdkitGraph, SmilesRdkitScratch,
};

/// Backward-compatible alias for the RDKit-normalized `smiles-parser` graph view.
pub type SmilesEcfpGraph<'a> = SmilesRdkitGraph<'a>;

/// Backward-compatible alias for the reusable RDKit-normalization scratch state.
pub type SmilesEcfpScratch = SmilesRdkitScratch;
