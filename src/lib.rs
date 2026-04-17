#![no_std]
#![doc = include_str!("../README.md")]

extern crate alloc;
#[cfg(test)]
extern crate std;

pub mod bit_fingerprint;
pub mod fingerprint;
pub mod fingerprints;
pub mod traits;

#[cfg(feature = "smiles-support")]
pub use crate::smiles_support::{SmilesEcfpGraph, SmilesEcfpScratch};
pub use crate::{
    bit_fingerprint::BitFingerprint,
    fingerprint::Fingerprint,
    fingerprints::EcfpFingerprint,
    traits::{EcfpGraph, MolecularAtom, MolecularBond, MolecularGraph},
};

/// Common imports for working with this crate.
pub mod prelude {
    pub use crate::{
        BitFingerprint, EcfpFingerprint, EcfpGraph, Fingerprint, MolecularAtom, MolecularBond,
        MolecularGraph,
    };
    #[cfg(feature = "smiles-support")]
    pub use crate::{SmilesEcfpGraph, SmilesEcfpScratch};
}

#[cfg(feature = "smiles-support")]
pub mod smiles_support;
#[cfg(any(test, feature = "smiles-support"))]
mod smiles_support_impl;
#[cfg(test)]
mod test_fixtures;
