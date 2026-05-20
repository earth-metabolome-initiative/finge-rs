//! The [`Mutator`] trait and the shared [`MutatorError`] returned when a
//! mutator declines to operate on a given graph.

use rand_core::RngCore;

use crate::{
    mutations::violation_class::ViolationClass,
    traits::{EcfpGraph, MolecularAtom},
};

/// A mutation operator that perturbs an [`EcfpGraph`] in a way that is, by
/// construction, visible to the ECFP output.
///
/// Mutators consume the input graph and return either a perturbed wrapper or a
/// [`MutatorError`] indicating that no eligible mutation exists for this input.
///
/// The trait method takes `rng: &mut dyn RngCore` rather than a generic `R`
/// parameter so that mutators remain object-safe and can be stored together in
/// a heterogeneous collection (see `MutatorMix`).
pub trait Mutator<G>
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    /// Concrete output graph type. Every built-in mutator returns
    /// `InvalidatedGraph<G>`; user-defined mutators may pick another type.
    type Output: EcfpGraph<NodeSymbol = G::NodeSymbol>;

    /// Applies the mutation, consuming `graph`.
    ///
    /// # Errors
    ///
    /// Returns a [`MutatorError`] when the input graph contains no candidate
    /// the mutator could target without producing a no-op.
    fn mutate(&self, graph: G, rng: &mut dyn RngCore) -> Result<Self::Output, MutatorError>;

    /// Returns the violation class this mutator is designed to introduce.
    fn primary_class(&self) -> ViolationClass;
}

/// Why a [`Mutator`] declined to operate on a given input.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MutatorError {
    /// No atom in the graph is a valid target for this mutator.
    NoEligibleAtom,
    /// No bond in the graph is a valid target for this mutator.
    NoEligibleBond,
    /// The mutator's target violation is already present in the input.
    AlreadyViolated,
    /// The graph is below the minimum size required by this mutator
    /// (e.g. mutators that need at least one bond reject single-atom graphs).
    GraphTooSmall,
}

impl core::fmt::Display for MutatorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoEligibleAtom => f.write_str("no eligible atom for mutation"),
            Self::NoEligibleBond => f.write_str("no eligible bond for mutation"),
            Self::AlreadyViolated => f.write_str("input already exhibits the target violation"),
            Self::GraphTooSmall => f.write_str("graph too small for this mutator"),
        }
    }
}

impl core::error::Error for MutatorError {}
