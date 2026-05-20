//! [`ViolationPredicate`] — pure checks for the eight ECFP-detectable
//! chemistry violations.
//!
//! This module exposes only the two trait shapes ([`PredicateClass`] +
//! [`ViolationPredicate`]) here. The eight concrete per-class
//! [`ViolationPredicate`] impls — and the matching standalone validity
//! free functions — are added in a follow-up commit.

use crate::{
    mutations::violation_class::ViolationClass,
    traits::{EcfpGraph, MolecularAtom},
};

/// Static identity of the [`ViolationClass`] a predicate detects.
///
/// Split out from [`ViolationPredicate`] (which is generic over `G`) so that
/// callers — and the `MutatorMix` builder — can recover the class without
/// having to pick a concrete graph type for type inference.
pub trait PredicateClass {
    /// Returns the violation class this predicate detects.
    fn class(&self) -> ViolationClass;
}

/// A pure-function detector for one [`ViolationClass`].
pub trait ViolationPredicate<G>: PredicateClass
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    /// Returns whether `graph` exhibits this predicate's violation.
    fn check(&self, graph: &G) -> bool;
}
