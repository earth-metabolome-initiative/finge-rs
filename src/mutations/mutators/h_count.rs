//! [`ImpossibleHCountMutator`] — overrides a random atom's `total_hydrogens`
//! past the predicate threshold of `max_natural_valence(Z) + 1`.

use rand_core::RngCore;

use crate::{
    mutations::{
        invalidated_graph::InvalidatedGraph,
        mutator::{Mutator, MutatorError},
        mutators::pick_unignored_atom,
        predicate::max_natural_valence,
        violation_class::ViolationClass,
    },
    traits::{EcfpGraph, MolecularAtom},
};

/// Mutator targeting [`ViolationClass::ImpossibleHCount`].
///
/// Overrides the chosen atom's `total_hydrogens` to
/// `max_natural_valence(Z) + Δ` with `Δ ∈ {2, 3, 4, 5}` so the predicate
/// (`total_hydrogens > max + 1`) reliably fires.
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleHCountMutator;

impl<G> Mutator<G> for ImpossibleHCountMutator
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    type Output = InvalidatedGraph<G>;

    fn mutate(&self, graph: G, rng: &mut dyn RngCore) -> Result<Self::Output, MutatorError> {
        let target = pick_unignored_atom(rng, &graph).ok_or(MutatorError::NoEligibleAtom)?;
        let fields = graph.ecfp_atom_invariant_fields(target);
        let bump = 2 + (rng.next_u32() % 4);
        let mut override_h = max_natural_valence(fields.atomic_number).saturating_add(bump);
        // Guarantee `override != current` — same defensive measure as in
        // `HypervalentMutator` for pathological inputs where the natural
        // value coincides with `max + bump`.
        if override_h == fields.total_hydrogens {
            override_h = override_h.saturating_add(1);
        }

        let mut wrapper = InvalidatedGraph::new(graph);
        wrapper.set_total_hydrogens_override(target, override_h);
        Ok(wrapper)
    }

    #[inline]
    fn primary_class(&self) -> ViolationClass {
        ViolationClass::ImpossibleHCount
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};

    use super::ImpossibleHCountMutator;
    use crate::{
        EcfpFingerprint, Fingerprint, Mutator, ViolationClass,
        mutations::predicate::{
            ImpossibleHCountPredicate, ViolationPredicate, has_impossible_hydrogen_count,
        },
        smiles_support_impl::SmilesRdkitScratch,
    };

    fn prepared(smiles: &str) -> (SmilesRdkitScratch, smiles_parser::smiles::Smiles) {
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("fixture SMILES should parse");
        (SmilesRdkitScratch::default(), parsed)
    }

    #[test]
    fn primary_class_is_impossible_h_count() {
        assert_eq!(
            <ImpossibleHCountMutator as Mutator<
                crate::smiles_support_impl::SmilesRdkitGraph<'_>,
            >>::primary_class(&ImpossibleHCountMutator),
            ViolationClass::ImpossibleHCount,
        );
    }

    #[test]
    fn mutator_returns_ok_and_predicate_fires() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut rng = ChaCha8Rng::seed_from_u64(0x10);
        let mutated = ImpossibleHCountMutator
            .mutate(inner, &mut rng)
            .expect("mutation should succeed");
        assert!(ImpossibleHCountPredicate.check(&mutated));
        assert!(has_impossible_hydrogen_count(&mutated));
    }

    // -----------------------------------------------------------------
    // Fuzz regression — see
    // `fuzz/artifacts/mutator_h_count/crash-3df93ad3...`.
    //
    // Same root cause as the atomic_number regressions: the mutator could
    // pick a target that ECFP later skipped via `ecfp_atom_is_ignored`
    // (the plain explicit `[H]`). Fixed by `pick_unignored_atom`.
    // -----------------------------------------------------------------

    #[test]
    fn mutator_returns_no_eligible_atom_when_all_atoms_ignored() {
        use crate::{InvalidatedGraph, MutatorError, traits::MolecularGraph as _};
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut all_ignored = InvalidatedGraph::new(inner);
        for atom_id in 0..all_ignored.atom_count() {
            all_ignored.set_atom_is_ignored_override(atom_id, true);
        }
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let err = ImpossibleHCountMutator
            .mutate(all_ignored, &mut rng)
            .expect_err("every atom is ignored; mutator should decline");
        assert_eq!(err, MutatorError::NoEligibleAtom);
    }

    #[test]
    fn fuzz_regression_3df93ad3_ignored_h_atom() {
        use crate::CountEcfpFingerprint;
        let parsed: smiles_parser::smiles::Smiles =
            "[H]s".parse().expect("fuzz-regression SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let inner = scratch.prepare(&parsed);
        let baseline = CountEcfpFingerprint::new(2, 65_536).compute(&inner);
        let mut rng = ChaCha8Rng::seed_from_u64(5_003_012_597_845_805_159);
        let mutated = ImpossibleHCountMutator
            .mutate(inner, &mut rng)
            .expect("mutator returned Err on a fuzz-regression input");
        let mutated_fp = CountEcfpFingerprint::new(2, 65_536).compute(&mutated);
        assert_ne!(baseline, mutated_fp);
    }

    #[test]
    fn ecfp_changes_at_radius_two() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline = EcfpFingerprint::new(2, 2048).compute(&inner);
        let mut rng = ChaCha8Rng::seed_from_u64(0x11);
        let mutated = ImpossibleHCountMutator
            .mutate(inner, &mut rng)
            .expect("mutation should succeed");
        let mutated_fp = EcfpFingerprint::new(2, 2048).compute(&mutated);
        assert_ne!(
            baseline.active_bits().collect::<Vec<_>>(),
            mutated_fp.active_bits().collect::<Vec<_>>(),
        );
    }
}
