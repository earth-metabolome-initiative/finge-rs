//! ECFP-detectable mutation operators for auxiliary negative-sampling tasks.
//!
//! See `docs/plan/` and `docs/research/ecfp-negative-sampling-survey.md` in the
//! repository for the design rationale.

pub mod invalidated_graph;
pub mod mutator;
pub mod mutators;
pub mod predicate;
pub mod violation_class;

pub use invalidated_graph::{AtomFieldOverride, InvalidatedGraph};
pub use mutator::{Mutator, MutatorError};
pub use mutators::{
    HypervalentMutator, ImpossibleAtomicNumberMutator, ImpossibleBondTypeMutator,
    ImpossibleChargeMutator, ImpossibleHCountMutator, ImpossibleIsotopeMutator,
    ImpossibleRingFlagMutator, TopologicalPathologyMutator,
};
pub use predicate::{
    HypervalentPredicate, ImpossibleAtomicNumberPredicate, ImpossibleBondTypePredicate,
    ImpossibleChargePredicate, ImpossibleHCountPredicate, ImpossibleIsotopePredicate,
    ImpossibleRingFlagPredicate, PredicateClass, TopologicalPathologyPredicate, ViolationPredicate,
    has_impossible_bond_type, has_impossible_charge, has_impossible_hydrogen_count,
    has_impossible_isotope, has_impossible_ring_flag, has_topological_pathology, is_hypervalent,
    is_impossible_atomic_number, max_natural_valence,
};
pub use violation_class::{ViolationClass, ViolationLabel};
