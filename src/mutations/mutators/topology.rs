//! [`TopologicalPathologyMutator`] — flips an aromatic-ring atom's `in_ring`
//! to `false`, breaking the cross-atom consistency the predicate detects.

use rand_core::RngCore;

use crate::{
    mutations::{
        invalidated_graph::InvalidatedGraph,
        mutator::{Mutator, MutatorError},
        mutators::{collect_atoms_matching, pick_from_slice},
        violation_class::ViolationClass,
    },
    traits::{EcfpGraph, MolecularAtom, MolecularBond},
};

/// Returns whether flipping `atom_id`'s `in_ring` to `false` will trigger
/// [`crate::mutations::predicate::has_topological_pathology`].
///
/// The predicate detects an aromatic bond (invariant `12`) whose endpoints
/// disagree on `in_ring`. To guarantee a fire we need `atom_id` to be
/// in-ring *and* to have an aromatic bond to a neighbour that is also
/// in-ring: flipping `atom_id` then creates the required disagreement on
/// that bond. A previous, looser check ("any aromatic bond") let the fuzzer
/// find SMILES `"nnB1BBnBcBF1:B"` with seed 33, where the chosen atom's
/// aromatic neighbour already had `in_ring = false` — the flip left both
/// endpoints non-ring and the predicate did not fire.
#[inline]
fn flipping_in_ring_creates_aromatic_inconsistency<G>(graph: &G, atom_id: usize) -> bool
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    if !graph.ecfp_atom_invariant_fields(atom_id).in_ring {
        return false;
    }
    graph.bonds(atom_id).any(|bond| {
        if graph.ecfp_bond_invariant(&bond, true) != 12 {
            return false;
        }
        let other = if bond.source() == atom_id {
            bond.target()
        } else {
            bond.source()
        };
        graph.ecfp_atom_invariant_fields(other).in_ring
    })
}

/// Mutator targeting [`ViolationClass::TopologicalPathology`].
///
/// Picks an aromatic-ring atom and overrides its `in_ring` to `false`.
/// The aromatic bond(s) touching this atom now have endpoints that disagree
/// on ring membership — exactly what
/// [`crate::mutations::predicate::has_topological_pathology`] detects.
#[derive(Clone, Copy, Debug, Default)]
pub struct TopologicalPathologyMutator;

impl<G> Mutator<G> for TopologicalPathologyMutator
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    type Output = InvalidatedGraph<G>;

    fn mutate(&self, graph: G, rng: &mut dyn RngCore) -> Result<Self::Output, MutatorError> {
        let atom_count = graph.atom_count();
        if atom_count == 0 {
            return Err(MutatorError::GraphTooSmall);
        }

        let candidates = collect_atoms_matching(atom_count, |atom_id| {
            !graph.ecfp_atom_is_ignored(atom_id)
                && flipping_in_ring_creates_aromatic_inconsistency(&graph, atom_id)
        });
        let target = *pick_from_slice(rng, &candidates).ok_or(MutatorError::NoEligibleAtom)?;

        let mut wrapper = InvalidatedGraph::new(graph);
        wrapper.set_in_ring_override(target, false);
        Ok(wrapper)
    }

    #[inline]
    fn primary_class(&self) -> ViolationClass {
        ViolationClass::TopologicalPathology
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};

    use super::TopologicalPathologyMutator;
    use crate::{
        EcfpFingerprint, Fingerprint, Mutator, MutatorError, ViolationClass,
        mutations::predicate::{
            TopologicalPathologyPredicate, ViolationPredicate, has_topological_pathology,
        },
        smiles_support_impl::SmilesRdkitScratch,
    };

    fn prepared(smiles: &str) -> (SmilesRdkitScratch, smiles_parser::smiles::Smiles) {
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("fixture SMILES should parse");
        (SmilesRdkitScratch::default(), parsed)
    }

    #[test]
    fn primary_class_is_topological_pathology() {
        assert_eq!(
            <TopologicalPathologyMutator as Mutator<
                crate::smiles_support_impl::SmilesRdkitGraph<'_>,
            >>::primary_class(&TopologicalPathologyMutator),
            ViolationClass::TopologicalPathology,
        );
    }

    #[test]
    fn mutator_succeeds_on_aromatic_molecule() {
        let (mut scratch, parsed) = prepared("c1ccccc1");
        let inner = scratch.prepare(&parsed);
        let mut rng = ChaCha8Rng::seed_from_u64(0x60);
        let mutated = TopologicalPathologyMutator
            .mutate(inner, &mut rng)
            .expect("benzene has aromatic ring atoms");
        assert!(TopologicalPathologyPredicate.check(&mutated));
        assert!(has_topological_pathology(&mutated));
    }

    #[test]
    fn mutator_fails_on_aliphatic_molecule() {
        // CCO has no aromatic bonds; no candidate ring atom exists.
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut rng = ChaCha8Rng::seed_from_u64(0x61);
        let err = TopologicalPathologyMutator
            .mutate(inner, &mut rng)
            .expect_err("CCO is aliphatic");
        assert_eq!(err, MutatorError::NoEligibleAtom);
    }

    #[test]
    fn mutator_fails_on_saturated_ring() {
        // Cyclohexane has ring atoms but no aromatic bonds.
        let (mut scratch, parsed) = prepared("C1CCCCC1");
        let inner = scratch.prepare(&parsed);
        let mut rng = ChaCha8Rng::seed_from_u64(0x62);
        let err = TopologicalPathologyMutator
            .mutate(inner, &mut rng)
            .expect_err("cyclohexane has no aromatic bonds");
        assert_eq!(err, MutatorError::NoEligibleAtom);
    }

    #[test]
    fn ecfp_changes_at_radius_two() {
        let (mut scratch, parsed) = prepared("c1ccccc1");
        let inner = scratch.prepare(&parsed);
        let baseline = EcfpFingerprint::new(2, 2048).compute(&inner);
        let mut rng = ChaCha8Rng::seed_from_u64(0x63);
        let mutated = TopologicalPathologyMutator
            .mutate(inner, &mut rng)
            .expect("mutation should succeed");
        let mutated_fp = EcfpFingerprint::new(2, 2048).compute(&mutated);
        assert_ne!(
            baseline.active_bits().collect::<Vec<_>>(),
            mutated_fp.active_bits().collect::<Vec<_>>(),
        );
    }

    // -----------------------------------------------------------------
    // Fuzz regression — see
    // `fuzz/artifacts/mutator_topology/crash-113b2b20...`.
    //
    // The mutator either declines the input or produces a wrapper for which
    // `TopologicalPathologyPredicate.check(&wrapper)` must return `true`.
    // The earlier fix (require aromatic-bonded neighbours to be in-ring)
    // wasn't enough for this case — surfaced by the fuzzer after the first
    // patch landed.
    // -----------------------------------------------------------------

    #[test]
    fn fuzz_regression_113b2b20_predicate_does_not_fire() {
        let parsed: smiles_parser::smiles::Smiles = "SSS:B7B-SSB-B7"
            .parse()
            .expect("fuzz-regression SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let inner = scratch.prepare(&parsed);

        let mut rng = ChaCha8Rng::seed_from_u64(151);
        let result = TopologicalPathologyMutator.mutate(inner, &mut rng);
        if let Ok(mutated) = result {
            // Whenever the mutator accepts an input, the matching predicate
            // must fire on the output. (If it declines, the regression is
            // trivially satisfied — see the existing "fails_on_aliphatic"
            // tests.)
            assert!(
                TopologicalPathologyPredicate.check(&mutated),
                "TopologicalPathologyMutator accepted SSS:B7B-SSB-B7 (seed 151) \
                 but the matching predicate did not fire",
            );
        }
    }
}
