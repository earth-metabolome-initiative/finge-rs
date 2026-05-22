//! [`ImpossibleChargeMutator`] — overrides a random atom's `formal_charge`
//! to a value outside the per-element plausibility window.

use rand_core::RngCore;

use crate::{
    mutations::{
        invalidated_graph::InvalidatedGraph,
        mutator::{Mutator, MutatorError},
        mutators::pick_unignored_atom,
        violation_class::ViolationClass,
    },
    traits::{EcfpGraph, MolecularAtom},
};

/// Inclusive bounds of the OOD charge-magnitude range. Predicate fires for
/// any `|q| > 6`; the lower bound of `7` preserves a one-step gap past that
/// threshold, and the upper bound is wide enough that the model can't
/// memorise a small set of magic values.
const CHARGE_MAGNITUDE_LOW: u32 = 7;
const CHARGE_MAGNITUDE_HIGH: u32 = 32_767;

/// Mutator targeting [`ViolationClass::ImpossibleCharge`].
///
/// Picks one random atom and writes a `formal_charge` of `sign × |q|` where
/// `|q| ~ Uniform[7, 32_767]` and `sign ∈ {-1, +1}` are drawn from `rng`.
/// All values have `|q| > 6` so the matching predicate always fires.
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleChargeMutator;

impl<G> Mutator<G> for ImpossibleChargeMutator
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
        let current = wrapper.ecfp_atom_invariant_fields(target).formal_charge;

        let r = rng.next_u32();
        let span = CHARGE_MAGNITUDE_HIGH - CHARGE_MAGNITUDE_LOW + 1;
        let magnitude = (CHARGE_MAGNITUDE_LOW + (r >> 1) % span) as i32;
        let sign = if r & 1 == 0 { 1 } else { -1 };
        let mut override_q = sign * magnitude;
        // Guarantee `override != current`. With magnitudes in `[7, 32767]`
        // the only collision happens when the inner already carries the
        // exact same charge — extremely unlikely for real molecules but
        // possible when composing several charge mutators on the same atom.
        // Negating preserves the OOD magnitude.
        if override_q == current {
            override_q = -override_q;
        }

        wrapper.set_formal_charge_override(target, override_q);
        Ok(())
    }

    #[inline]
    fn primary_class(&self) -> ViolationClass {
        ViolationClass::ImpossibleCharge
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};

    use super::ImpossibleChargeMutator;
    use crate::{
        EcfpFingerprint, EcfpGraph, Fingerprint, InvalidatedGraph, Mutator, ViolationClass,
        mutations::predicate::{
            ImpossibleChargePredicate, ViolationPredicate, has_impossible_charge,
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
    fn primary_class_is_impossible_charge() {
        assert_eq!(
            <ImpossibleChargeMutator as Mutator<
                crate::smiles_support_impl::SmilesRdkitGraph<'_>,
            >>::primary_class(&ImpossibleChargeMutator),
            ViolationClass::ImpossibleCharge,
        );
    }

    fn mutate<'a>(
        inner: crate::smiles_support_impl::SmilesRdkitGraph<'a>,
        seed: u64,
    ) -> InvalidatedGraph<crate::smiles_support_impl::SmilesRdkitGraph<'a>> {
        let mut wrapper = InvalidatedGraph::new(inner);
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        ImpossibleChargeMutator
            .mutate_in_place(&mut wrapper, &mut rng)
            .expect("mutation should succeed");
        wrapper
    }

    #[test]
    fn mutator_returns_ok_and_predicate_fires() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mutated = mutate(inner, 0x20);
        assert!(ImpossibleChargePredicate.check(&mutated));
        assert!(has_impossible_charge(&mutated));
    }

    // -----------------------------------------------------------------
    // Fuzz regression — see
    // `fuzz/artifacts/mutator_formal_charge/crash-e87be2c6...`.
    //
    // Same root cause as the atomic_number regressions: the mutator could
    // pick a target that ECFP later skipped via `ecfp_atom_is_ignored`
    // (here the chiral explicit `[H@]`). Fixed by `pick_unignored_atom`.
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
        let err = ImpossibleChargeMutator
            .mutate_in_place(&mut all_ignored, &mut rng)
            .expect_err("every atom is ignored; mutator should decline");
        assert_eq!(err, MutatorError::NoEligibleAtom);
    }

    #[test]
    fn fuzz_regression_e87be2c6_ignored_h_atom() {
        use crate::CountEcfpFingerprint;
        let parsed: smiles_parser::smiles::Smiles = "[H@]C"
            .parse()
            .expect("fuzz-regression SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let inner = scratch.prepare(&parsed);
        let baseline = CountEcfpFingerprint::new(2, 65_536).compute(&inner);
        let mut wrapper = InvalidatedGraph::new(inner);
        let mut rng = ChaCha8Rng::seed_from_u64(5_327_177);
        ImpossibleChargeMutator
            .mutate_in_place(&mut wrapper, &mut rng)
            .expect("mutator returned Err on a fuzz-regression input");
        let mutated_fp = CountEcfpFingerprint::new(2, 65_536).compute(&wrapper);
        assert_ne!(baseline, mutated_fp);
    }

    #[test]
    fn override_nudges_off_pre_existing_collision() {
        // Same deterministic-collision trick as `HypervalentMutator`: run
        // once to learn the RNG-chosen charge, then re-run with that value
        // pre-set so the adjustment branch (negation) fires.
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);

        let mut probe = InvalidatedGraph::new(inner);
        let mut rng = ChaCha8Rng::seed_from_u64(0xC011_5104);
        ImpossibleChargeMutator
            .mutate_in_place(&mut probe, &mut rng)
            .expect("probe mutation should succeed");
        let target = (0..probe.atom_count())
            .find(|&id| probe.atom_field_override(id).is_some())
            .expect("probe should have written one override");
        let first_value = probe.ecfp_atom_invariant_fields(target).formal_charge;

        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_formal_charge_override(target, first_value);
        let mut rng2 = ChaCha8Rng::seed_from_u64(0xC011_5104);
        ImpossibleChargeMutator
            .mutate_in_place(&mut wrapper, &mut rng2)
            .expect("collision-rerun mutation should succeed");
        let second_value = wrapper.ecfp_atom_invariant_fields(target).formal_charge;
        assert_ne!(
            second_value, first_value,
            "override-equals-current branch must have negated the charge",
        );
    }

    #[test]
    fn ecfp_changes_at_radius_two() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline = EcfpFingerprint::new(2, 2048).compute(&inner);
        let mutated = mutate(inner, 0x21);
        let mutated_fp = EcfpFingerprint::new(2, 2048).compute(&mutated);
        assert_ne!(
            baseline.active_bits().collect::<Vec<_>>(),
            mutated_fp.active_bits().collect::<Vec<_>>(),
        );
    }
}
