//! [`ImpossibleRingFlagMutator`] — flips `in_ring = true` on an atom that
//! cannot possibly be in a ring (heavy degree `< 2`).

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

#[inline]
fn heavy_neighbour_count<G>(graph: &G, atom_id: usize) -> usize
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    graph
        .bonds(atom_id)
        .filter(|bond| {
            let other = if bond.source() == atom_id {
                bond.target()
            } else {
                bond.source()
            };
            graph.ecfp_atom_invariant_fields(other).atomic_number != 1
        })
        .count()
}

/// Mutator targeting [`ViolationClass::ImpossibleRingFlag`].
///
/// Picks an atom with fewer than two heavy neighbours (a topological
/// terminus that cannot be in a cycle) and forces its `in_ring` field to
/// `true`. The matching predicate then fires.
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleRingFlagMutator;

impl<G> Mutator<G> for ImpossibleRingFlagMutator
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
            if graph.ecfp_atom_is_ignored(atom_id) {
                return false;
            }
            let fields = graph.ecfp_atom_invariant_fields(atom_id);
            !fields.in_ring && heavy_neighbour_count(&graph, atom_id) < 2
        });
        let target = *pick_from_slice(rng, &candidates).ok_or(MutatorError::NoEligibleAtom)?;

        let mut wrapper = InvalidatedGraph::new(graph);
        wrapper.set_in_ring_override(target, true);
        Ok(wrapper)
    }

    #[inline]
    fn primary_class(&self) -> ViolationClass {
        ViolationClass::ImpossibleRingFlag
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};

    use super::ImpossibleRingFlagMutator;
    use crate::{
        EcfpFingerprint, Fingerprint, Mutator, MutatorError, ViolationClass,
        mutations::predicate::{
            ImpossibleRingFlagPredicate, ViolationPredicate, has_impossible_ring_flag,
        },
        smiles_support_impl::SmilesRdkitScratch,
    };

    fn prepared(smiles: &str) -> (SmilesRdkitScratch, smiles_parser::smiles::Smiles) {
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("fixture SMILES should parse");
        (SmilesRdkitScratch::default(), parsed)
    }

    #[test]
    fn primary_class_is_impossible_ring_flag() {
        assert_eq!(
            <ImpossibleRingFlagMutator as Mutator<
                crate::smiles_support_impl::SmilesRdkitGraph<'_>,
            >>::primary_class(&ImpossibleRingFlagMutator),
            ViolationClass::ImpossibleRingFlag,
        );
    }

    #[test]
    fn mutator_succeeds_on_molecule_with_terminal_atom() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut rng = ChaCha8Rng::seed_from_u64(0x40);
        let mutated = ImpossibleRingFlagMutator
            .mutate(inner, &mut rng)
            .expect("CCO has terminal atoms; mutation should succeed");
        assert!(ImpossibleRingFlagPredicate.check(&mutated));
        assert!(has_impossible_ring_flag(&mutated));
    }

    #[test]
    fn mutator_fails_when_no_terminal_atom_exists() {
        // Every atom in benzene has heavy_neighbour_count == 2.
        let (mut scratch, parsed) = prepared("c1ccccc1");
        let inner = scratch.prepare(&parsed);
        let mut rng = ChaCha8Rng::seed_from_u64(0x41);
        let err = ImpossibleRingFlagMutator
            .mutate(inner, &mut rng)
            .expect_err("benzene has no terminus");
        assert_eq!(err, MutatorError::NoEligibleAtom);
    }

    // -----------------------------------------------------------------
    // Fuzz regression — see
    // `fuzz/artifacts/mutator_ring_flag/crash-cb578195...`.
    //
    // Same root cause as the atomic_number regressions: the mutator's
    // candidate filter didn't exclude `ecfp_atom_is_ignored` atoms, so it
    // could land on a collapsed `[H]` whose ring-flag override was
    // invisible to ECFP. Fixed by extending the
    // `collect_atoms_matching` filter with `!graph.ecfp_atom_is_ignored`.
    // -----------------------------------------------------------------

    #[test]
    fn fuzz_regression_cb578195_ignored_h_atom() {
        use crate::CountEcfpFingerprint;
        let parsed: smiles_parser::smiles::Smiles = "bFFFFCsBF.cCsFsCI2NF2P[H]"
            .parse()
            .expect("fuzz-regression SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let inner = scratch.prepare(&parsed);
        let baseline = CountEcfpFingerprint::new(2, 65_536).compute(&inner);
        let mut rng = ChaCha8Rng::seed_from_u64(8_092_465_109_790_371_496);
        // Mutator may decline if no eligible atom — but if it accepts, ECFP
        // must differ.
        if let Ok(mutated) = ImpossibleRingFlagMutator.mutate(inner, &mut rng) {
            let mutated_fp = CountEcfpFingerprint::new(2, 65_536).compute(&mutated);
            assert_ne!(baseline, mutated_fp);
        }
    }

    #[test]
    fn ecfp_changes_at_radius_two() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline = EcfpFingerprint::new(2, 2048).compute(&inner);
        let mut rng = ChaCha8Rng::seed_from_u64(0x42);
        let mutated = ImpossibleRingFlagMutator
            .mutate(inner, &mut rng)
            .expect("mutation should succeed");
        let mutated_fp = EcfpFingerprint::new(2, 2048).compute(&mutated);
        assert_ne!(
            baseline.active_bits().collect::<Vec<_>>(),
            mutated_fp.active_bits().collect::<Vec<_>>(),
        );
    }
}
