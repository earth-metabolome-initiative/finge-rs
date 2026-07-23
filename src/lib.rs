#![no_std]
#![doc = include_str!("../README.md")]

extern crate alloc;
#[cfg(test)]
extern crate std;

pub use minhash_rs::prelude::{MinHash, MinHasher};

pub mod atom_invariant_fields;
pub mod bit_fingerprint;
pub mod count_fingerprint;
pub mod fingerprint;
pub mod fingerprints;
pub mod lsh;
pub mod mutations;
pub mod traits;

#[cfg(feature = "smarts-support")]
pub use crate::fingerprints::MaccsFingerprint;
#[cfg(feature = "smarts-support")]
pub use crate::maccs_support::{
    MACCS_KEY_COUNT, MaccsBuildError, MaccsKeyDefinition, MaccsSpecialCase,
    RDKIT_MACCS_THRESHOLD_KEY_IDS, compile_rdkit_maccs_queries, has_multiple_aromatic_rings,
    has_multiple_fragments, rdkit_maccs_keys,
};
pub use crate::smiles_support::{
    SmilesEcfpGraph, SmilesEcfpScratch, SmilesPreparationError, SmilesRdkitGraph,
    SmilesRdkitScratch,
};
pub use crate::{
    atom_invariant_fields::{AtomInvariantFields, morgan_hash_combine, morgan_hash_u32_sequence},
    bit_fingerprint::BitFingerprint,
    count_fingerprint::{CountFingerprint, LayeredCountFingerprint},
    fingerprint::Fingerprint,
    fingerprints::{
        AtomPairFingerprint, CountEcfpFingerprint, CountMap4Fingerprint, DEFAULT_MAP4_FP_SIZE,
        DEFAULT_MAP4_RADIUS, EcfpFingerprint, LayeredCountEcfpFingerprint,
        MAP4_DISCONNECTED_DISTANCE, Map4Fingerprint, RdkFingerprint, TopologicalTorsionFingerprint,
    },
    lsh::{LshIndex, Sketcher},
    mutations::{
        AtomFieldOverride, HypervalentMutator, HypervalentPredicate, ImpossibleAtomicNumberMutator,
        ImpossibleAtomicNumberPredicate, ImpossibleBondTypeMutator, ImpossibleBondTypePredicate,
        ImpossibleChargeMutator, ImpossibleChargePredicate, ImpossibleHCountMutator,
        ImpossibleHCountPredicate, ImpossibleIsotopeMutator, ImpossibleIsotopePredicate,
        ImpossibleRingFlagMutator, ImpossibleRingFlagPredicate, InvalidatedGraph, Mutator,
        MutatorError, MutatorMix, PredicateClass, TopologicalPathologyMutator,
        TopologicalPathologyPredicate, ViolationClass, ViolationLabel, ViolationPredicate,
        has_impossible_bond_type, has_impossible_charge, has_impossible_hydrogen_count,
        has_impossible_isotope, has_impossible_ring_flag, has_topological_pathology,
        is_hypervalent, is_impossible_atomic_number, max_natural_valence,
    },
    traits::{
        AtomPairGraph, EcfpGraph, Map4Graph, MolecularAtom, MolecularBond, MolecularGraph,
        RdkFingerprintGraph, TopologicalTorsionGraph,
    },
};

/// Common imports for working with this crate.
pub mod prelude {
    #[cfg(feature = "smarts-support")]
    pub use crate::MaccsFingerprint;
    pub use crate::{
        AtomFieldOverride, AtomInvariantFields, AtomPairFingerprint, AtomPairGraph, BitFingerprint,
        CountEcfpFingerprint, CountFingerprint, CountMap4Fingerprint, EcfpFingerprint, EcfpGraph,
        Fingerprint, HypervalentMutator, HypervalentPredicate, ImpossibleAtomicNumberMutator,
        ImpossibleAtomicNumberPredicate, ImpossibleBondTypeMutator, ImpossibleBondTypePredicate,
        ImpossibleChargeMutator, ImpossibleChargePredicate, ImpossibleHCountMutator,
        ImpossibleHCountPredicate, ImpossibleIsotopeMutator, ImpossibleIsotopePredicate,
        ImpossibleRingFlagMutator, ImpossibleRingFlagPredicate, InvalidatedGraph,
        LayeredCountEcfpFingerprint, LayeredCountFingerprint, Map4Fingerprint, Map4Graph,
        MolecularAtom, MolecularBond, MolecularGraph, Mutator, MutatorError, MutatorMix,
        PredicateClass, RdkFingerprint, RdkFingerprintGraph, TopologicalPathologyMutator,
        TopologicalPathologyPredicate, TopologicalTorsionFingerprint, TopologicalTorsionGraph,
        ViolationClass, ViolationLabel, ViolationPredicate, has_impossible_bond_type,
        has_impossible_charge, has_impossible_hydrogen_count, has_impossible_isotope,
        has_impossible_ring_flag, has_topological_pathology, is_hypervalent,
        is_impossible_atomic_number, max_natural_valence,
    };
    pub use crate::{
        SmilesEcfpGraph, SmilesEcfpScratch, SmilesPreparationError, SmilesRdkitGraph,
        SmilesRdkitScratch,
    };
}

#[cfg(feature = "smarts-support")]
pub mod maccs_support;
pub mod smiles_support;
mod smiles_support_impl;
#[cfg(test)]
mod test_fixtures;
