//! The eight built-in [`crate::Mutator`] implementations, one per
//! [`crate::ViolationClass`].
//!
//! Each mutator targets exactly one ECFP invariant channel (atomic number,
//! degree, hydrogen count, formal charge, isotope, ring flag, bond type, or a
//! composite topological violation) by writing a field-level override to the
//! resulting [`crate::InvalidatedGraph`]. Standard ECFP then runs over the
//! perturbed fields, so the resulting fingerprint differs from the original
//! and lives in a chemistry-shaped region of invariant space (rather than
//! carrying a synthetic marker).

pub mod atomic_number;
pub mod bond_type;
pub mod formal_charge;
pub mod h_count;
pub mod hypervalent;
pub mod isotope;
pub mod ring_flag;
pub mod topology;

pub use atomic_number::ImpossibleAtomicNumberMutator;
pub use bond_type::ImpossibleBondTypeMutator;
pub use formal_charge::ImpossibleChargeMutator;
pub use h_count::ImpossibleHCountMutator;
pub use hypervalent::HypervalentMutator;
pub use isotope::ImpossibleIsotopeMutator;
pub use ring_flag::ImpossibleRingFlagMutator;
pub use topology::TopologicalPathologyMutator;

use alloc::vec::Vec;

use rand_core::RngCore;

use crate::traits::{EcfpGraph, MolecularAtom};

/// Picks a uniform random element from `candidates`.
///
/// Returns `None` when `candidates` is empty.
#[inline]
pub(crate) fn pick_from_slice<'a, T>(rng: &mut dyn RngCore, candidates: &'a [T]) -> Option<&'a T> {
    if candidates.is_empty() {
        None
    } else {
        Some(&candidates[(rng.next_u32() as usize) % candidates.len()])
    }
}

/// Collects every `atom_id` in `0..atom_count` for which `predicate` returns
/// `true`.
#[inline]
pub(crate) fn collect_atoms_matching<F>(atom_count: usize, predicate: F) -> Vec<usize>
where
    F: Fn(usize) -> bool,
{
    (0..atom_count).filter(|&id| predicate(id)).collect()
}

/// Picks a uniform random atom id from the subset of atoms that ECFP will
/// actually emit features for — i.e. those for which
/// [`EcfpGraph::ecfp_atom_is_ignored`] returns `false`.
///
/// Atom-channel mutators must use this instead of a raw `0..atom_count`
/// pick: an ignored atom (typically a collapsed plain explicit `[H]` in
/// `SmilesRdkitGraph`) is skipped by the ECFP radius-0 loop, so a field
/// override on it is silently invisible — ECFP and every fingerprint
/// derived from it stay byte-identical to the baseline. Fuzz crashes
/// `mutator_atomic_number/crash-559b89bb…`, `…/crash-75a084e2…`,
/// `…/crash-bd77b81f…`, `…/crash-fcec3d66…`, `mutator_isotope/crash-7c50603e…`,
/// and `…/crash-a0439b0c…` all reduced to this single pattern.
///
/// Returns `None` when there are no unignored atoms (empty graph, or every
/// atom is ignored). The caller should map that to
/// [`crate::MutatorError::NoEligibleAtom`].
#[inline]
pub(crate) fn pick_unignored_atom<G>(rng: &mut dyn RngCore, graph: &G) -> Option<usize>
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    let candidates: Vec<usize> = (0..graph.atom_count())
        .filter(|&id| !graph.ecfp_atom_is_ignored(id))
        .collect();
    pick_from_slice(rng, &candidates).copied()
}
