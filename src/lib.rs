#![no_std]
#![doc = include_str!("../README.md")]

extern crate alloc;
#[cfg(test)]
extern crate std;

pub mod bit_fingerprint;
pub mod count_fingerprint;
pub mod fingerprint;
pub mod fingerprints;
pub mod traits;

#[cfg(feature = "smarts-support")]
pub use crate::fingerprints::MaccsFingerprint;
#[cfg(feature = "smarts-support")]
pub use crate::maccs_support::{
    MACCS_KEY_COUNT, MaccsBuildError, MaccsKeyDefinition, MaccsSpecialCase,
    RDKIT_MACCS_THRESHOLD_KEY_IDS, compile_rdkit_maccs_queries, has_multiple_aromatic_rings,
    has_multiple_fragments, rdkit_maccs_keys,
};
#[cfg(feature = "smiles-support")]
pub use crate::smiles_support::{
    SmilesEcfpGraph, SmilesEcfpScratch, SmilesPreparationError, SmilesRdkitGraph,
    SmilesRdkitScratch,
};
pub use crate::{
    bit_fingerprint::BitFingerprint,
    count_fingerprint::{CountFingerprint, LayeredCountFingerprint},
    fingerprint::Fingerprint,
    fingerprints::{
        AtomPairFingerprint, CountEcfpFingerprint, EcfpFingerprint, LayeredCountEcfpFingerprint,
        TopologicalTorsionFingerprint,
    },
    traits::{
        AtomPairGraph, EcfpGraph, MolecularAtom, MolecularBond, MolecularGraph,
        TopologicalTorsionGraph,
    },
};

/// Common imports for working with this crate.
pub mod prelude {
    #[cfg(feature = "smarts-support")]
    pub use crate::MaccsFingerprint;
    pub use crate::{
        AtomPairFingerprint, AtomPairGraph, BitFingerprint, CountEcfpFingerprint, CountFingerprint,
        EcfpFingerprint, EcfpGraph, Fingerprint, LayeredCountEcfpFingerprint,
        LayeredCountFingerprint, MolecularAtom, MolecularBond, MolecularGraph,
        TopologicalTorsionFingerprint, TopologicalTorsionGraph,
    };
    #[cfg(feature = "smiles-support")]
    pub use crate::{
        SmilesEcfpGraph, SmilesEcfpScratch, SmilesPreparationError, SmilesRdkitGraph,
        SmilesRdkitScratch,
    };
}

#[cfg(feature = "smarts-support")]
pub mod maccs_support;
#[cfg(feature = "smiles-support")]
pub mod smiles_support;
#[cfg(any(test, feature = "smiles-support"))]
mod smiles_support_impl;
#[cfg(test)]
mod test_fixtures;
