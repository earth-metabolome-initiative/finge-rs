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

/// Ceiling of the OOD hydrogen-count range. Predicate fires for any value
/// `> max_natural_valence(Z) + 1`; the lower bound `max + 2` preserves a
/// one-step gap past that threshold.
const HYDROGENS_CEILING: u32 = 200;

/// Mutator targeting [`ViolationClass::ImpossibleHCount`].
///
/// Overrides the chosen atom's `total_hydrogens` to a uniform random integer
/// in `[max_natural_valence(Z) + 2, HYDROGENS_CEILING]`. The wide upper
/// bound spreads the OOD hash region so the model can't memorise a small
/// set of magic counts.
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleHCountMutator;

impl<G> Mutator<G> for ImpossibleHCountMutator
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    fn mutate_in_place(
        &self,
        wrapper: &mut InvalidatedGraph<G>,
        rng: &mut dyn RngCore,
    ) -> Result<(), MutatorError> {
        let target = pick_unignored_atom(rng, wrapper).ok_or(MutatorError::NoEligibleAtom)?;
        let fields = wrapper.ecfp_atom_invariant_fields(target);
        let lower = max_natural_valence(fields.atomic_number).saturating_add(2);
        let upper = HYDROGENS_CEILING.max(lower);
        let span = upper - lower + 1;
        let mut override_h = lower + rng.next_u32() % span;
        // Guarantee `override != current` — same defensive measure as in
        // `HypervalentMutator` for pathological inputs where the natural
        // value coincides with the uniform draw.
        if override_h == fields.total_hydrogens {
            override_h = if override_h < upper {
                override_h + 1
            } else {
                lower
            };
        }

        wrapper.set_total_hydrogens_override(target, override_h);
        Ok(())
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
        EcfpFingerprint, EcfpGraph, Fingerprint, InvalidatedGraph, Mutator, ViolationClass,
        mutations::predicate::{
            ImpossibleHCountPredicate, ViolationPredicate, has_impossible_hydrogen_count,
        },
        smiles_support_impl::SmilesRdkitScratch,
        traits::MolecularGraph as _,
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

    fn mutate<'a>(
        inner: crate::smiles_support_impl::SmilesRdkitGraph<'a>,
        seed: u64,
    ) -> InvalidatedGraph<crate::smiles_support_impl::SmilesRdkitGraph<'a>> {
        let mut wrapper = InvalidatedGraph::new(inner);
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        ImpossibleHCountMutator
            .mutate_in_place(&mut wrapper, &mut rng)
            .expect("mutation should succeed");
        wrapper
    }

    #[test]
    fn mutator_returns_ok_and_predicate_fires() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mutated = mutate(inner, 0x10);
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
        use crate::{MutatorError, traits::MolecularGraph as _};
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut all_ignored = InvalidatedGraph::new(inner);
        for atom_id in 0..all_ignored.atom_count() {
            all_ignored.set_atom_is_ignored_override(atom_id, true);
        }
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let err = ImpossibleHCountMutator
            .mutate_in_place(&mut all_ignored, &mut rng)
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
        let mut wrapper = InvalidatedGraph::new(inner);
        let mut rng = ChaCha8Rng::seed_from_u64(5_003_012_597_845_805_159);
        ImpossibleHCountMutator
            .mutate_in_place(&mut wrapper, &mut rng)
            .expect("mutator returned Err on a fuzz-regression input");
        let mutated_fp = CountEcfpFingerprint::new(2, 65_536).compute(&wrapper);
        assert_ne!(baseline, mutated_fp);
    }

    #[test]
    fn override_nudges_off_pre_existing_collision() {
        // Same deterministic-collision trick as `HypervalentMutator`: run
        // once to learn the RNG-chosen override, then re-run with that value
        // pre-set so the adjustment branch fires.
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);

        let mut probe = InvalidatedGraph::new(inner);
        let mut rng = ChaCha8Rng::seed_from_u64(0xC011_5104);
        ImpossibleHCountMutator
            .mutate_in_place(&mut probe, &mut rng)
            .expect("probe mutation should succeed");
        let target = (0..probe.atom_count())
            .find(|&id| probe.atom_field_override(id).is_some())
            .expect("probe should have written one override");
        let first_value = probe.ecfp_atom_invariant_fields(target).total_hydrogens;

        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_total_hydrogens_override(target, first_value);
        let mut rng2 = ChaCha8Rng::seed_from_u64(0xC011_5104);
        ImpossibleHCountMutator
            .mutate_in_place(&mut wrapper, &mut rng2)
            .expect("collision-rerun mutation should succeed");
        let second_value = wrapper.ecfp_atom_invariant_fields(target).total_hydrogens;
        assert_ne!(
            second_value, first_value,
            "override-equals-current branch must have nudged the value off",
        );
    }

    #[test]
    fn ecfp_changes_at_radius_two() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline = EcfpFingerprint::new(2, 2048).compute(&inner);
        let mutated = mutate(inner, 0x11);
        let mutated_fp = EcfpFingerprint::new(2, 2048).compute(&mutated);
        assert_ne!(
            baseline.active_bits().collect::<Vec<_>>(),
            mutated_fp.active_bits().collect::<Vec<_>>(),
        );
    }
}
