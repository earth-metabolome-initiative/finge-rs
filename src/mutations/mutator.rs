//! The [`Mutator`] trait and the shared [`MutatorError`] returned when a
//! mutator declines to operate on a given graph.

use rand_core::RngCore;

use crate::{
    mutations::{invalidated_graph::InvalidatedGraph, violation_class::ViolationClass},
    traits::{EcfpGraph, MolecularAtom},
};

/// A mutation operator that perturbs an [`EcfpGraph`] in a way that is, by
/// construction, visible to the ECFP output.
///
/// Mutators write field-level overrides into an [`InvalidatedGraph`] wrapper
/// shared across one or more invocations, enabling [`crate::MutatorMix`] to
/// compose several mutations on the same molecule (possibly on the same atom).
/// They return either `Ok(())` or a [`MutatorError`] when no eligible target
/// exists.
///
/// The trait method takes `rng: &mut dyn RngCore` rather than a generic `R`
/// parameter so that mutators remain object-safe and can be stored together in
/// a heterogeneous collection (see [`crate::MutatorMix`]).
pub trait Mutator<G>
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    /// Applies the mutation in place, writing overrides into `wrapper`.
    ///
    /// Mutators observe the wrapper's *current* (possibly already-mutated)
    /// state when selecting targets, so composing several mutations means
    /// each mutator picks against the cumulative perturbation.
    ///
    /// # Errors
    ///
    /// Returns a [`MutatorError`] when the wrapper contains no candidate the
    /// mutator could target without producing a no-op.
    fn mutate_in_place(
        &self,
        wrapper: &mut InvalidatedGraph<G>,
        rng: &mut dyn RngCore,
    ) -> Result<(), MutatorError>;

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
    /// Returned by [`crate::MutatorMix::sample`] when every successfully
    /// applied mutation in a composition was clobbered by a later mutation
    /// on the same atom-channel, leaving the wrapper pre-hash-identical to
    /// the inner graph. The caller should retry with a fresh seed or input.
    CompositionCancelled,
    /// Returned by [`crate::MutatorMix::sample`] when the collision filter
    /// (configured via
    /// [`crate::MutatorMix::with_collision_filter`]) detected that the
    /// wrapper's folded ECFP is byte-equal to the baseline's at the
    /// configured `(radius, fp_size)`. The mutation changed pre-hash atom
    /// invariants but the change collapsed back to identical bins under
    /// folding — the sample carries no signal at the chosen `fp_size`.
    /// The caller should retry with a fresh seed or input.
    FingerprintCollision,
}

impl core::fmt::Display for MutatorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoEligibleAtom => f.write_str("no eligible atom for mutation"),
            Self::NoEligibleBond => f.write_str("no eligible bond for mutation"),
            Self::AlreadyViolated => f.write_str("input already exhibits the target violation"),
            Self::GraphTooSmall => f.write_str("graph too small for this mutator"),
            Self::CompositionCancelled => {
                f.write_str("composed mutations cancelled out, no effective perturbation")
            }
            Self::FingerprintCollision => {
                f.write_str("wrapper ECFP folded to the same bins as the baseline")
            }
        }
    }
}

impl core::error::Error for MutatorError {}

#[cfg(test)]
mod tests {
    use alloc::string::ToString;

    use super::MutatorError;

    #[test]
    fn display_covers_every_error_variant() {
        assert_eq!(
            MutatorError::NoEligibleAtom.to_string(),
            "no eligible atom for mutation",
        );
        assert_eq!(
            MutatorError::NoEligibleBond.to_string(),
            "no eligible bond for mutation",
        );
        assert_eq!(
            MutatorError::AlreadyViolated.to_string(),
            "input already exhibits the target violation",
        );
        assert_eq!(
            MutatorError::GraphTooSmall.to_string(),
            "graph too small for this mutator",
        );
        assert_eq!(
            MutatorError::CompositionCancelled.to_string(),
            "composed mutations cancelled out, no effective perturbation",
        );
        assert_eq!(
            MutatorError::FingerprintCollision.to_string(),
            "wrapper ECFP folded to the same bins as the baseline",
        );
    }
}
