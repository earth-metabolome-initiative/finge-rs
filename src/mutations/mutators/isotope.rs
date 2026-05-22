//! [`ImpossibleIsotopeMutator`] — overrides a random atom's `delta_mass`
//! (isotope offset) far past any real-isotope value.

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

/// Inclusive bounds of the OOD `delta_mass` magnitude range. Predicate fires
/// for any `|Δmass| > 80`; the lower bound of `81` preserves a one-step gap
/// past that threshold, and the upper bound is wide enough that the model
/// can't memorise a small set of magic isotopes.
const DELTA_MASS_MAGNITUDE_LOW: u32 = 81;
const DELTA_MASS_MAGNITUDE_HIGH: u32 = 32_767;

/// Mutator targeting [`ViolationClass::ImpossibleIsotope`].
///
/// Picks one random atom and writes `delta_mass = sign × |Δm|` where
/// `|Δm| ~ Uniform[81, 32_767]` and `sign ∈ {-1, +1}` are drawn from `rng`.
/// All values have `|Δmass| > 80` so the matching predicate always fires.
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleIsotopeMutator;

impl<G> Mutator<G> for ImpossibleIsotopeMutator
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
        let current = wrapper.ecfp_atom_invariant_fields(target).delta_mass;

        let r = rng.next_u32();
        let span = DELTA_MASS_MAGNITUDE_HIGH - DELTA_MASS_MAGNITUDE_LOW + 1;
        let magnitude = (DELTA_MASS_MAGNITUDE_LOW + (r >> 1) % span) as i32;
        let sign = if r & 1 == 0 { 1 } else { -1 };
        let mut override_dm = sign * magnitude;
        // Guarantee `override != current` — same defensive scheme as
        // `ImpossibleChargeMutator`.
        if override_dm == current {
            override_dm = -override_dm;
        }

        wrapper.set_delta_mass_override(target, override_dm);
        Ok(())
    }

    #[inline]
    fn primary_class(&self) -> ViolationClass {
        ViolationClass::ImpossibleIsotope
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};

    use super::ImpossibleIsotopeMutator;
    use crate::{
        EcfpFingerprint, Fingerprint, InvalidatedGraph, Mutator, ViolationClass,
        mutations::predicate::{
            ImpossibleIsotopePredicate, ViolationPredicate, has_impossible_isotope,
        },
        smiles_support_impl::SmilesRdkitScratch,
    };

    fn prepared(smiles: &str) -> (SmilesRdkitScratch, smiles_parser::smiles::Smiles) {
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("fixture SMILES should parse");
        (SmilesRdkitScratch::default(), parsed)
    }

    #[test]
    fn primary_class_is_impossible_isotope() {
        assert_eq!(
            <ImpossibleIsotopeMutator as Mutator<
                crate::smiles_support_impl::SmilesRdkitGraph<'_>,
            >>::primary_class(&ImpossibleIsotopeMutator),
            ViolationClass::ImpossibleIsotope,
        );
    }

    fn mutate<'a>(
        inner: crate::smiles_support_impl::SmilesRdkitGraph<'a>,
        seed: u64,
    ) -> InvalidatedGraph<crate::smiles_support_impl::SmilesRdkitGraph<'a>> {
        let mut wrapper = InvalidatedGraph::new(inner);
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        ImpossibleIsotopeMutator
            .mutate_in_place(&mut wrapper, &mut rng)
            .expect("mutation should succeed");
        wrapper
    }

    #[test]
    fn mutator_returns_ok_and_predicate_fires() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mutated = mutate(inner, 0x30);
        assert!(ImpossibleIsotopePredicate.check(&mutated));
        assert!(has_impossible_isotope(&mutated));
    }

    #[test]
    fn ecfp_changes_at_radius_two() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline = EcfpFingerprint::new(2, 2048).compute(&inner);
        let mutated = mutate(inner, 0x31);
        let mutated_fp = EcfpFingerprint::new(2, 2048).compute(&mutated);
        assert_ne!(
            baseline.active_bits().collect::<Vec<_>>(),
            mutated_fp.active_bits().collect::<Vec<_>>(),
        );
    }

    // -----------------------------------------------------------------
    // Fuzz regressions — see `fuzz/artifacts/mutator_isotope/`.
    //
    // Same root cause as the atomic_number regressions: the mutator can
    // pick a target atom that ECFP later skips via
    // `ecfp_atom_is_ignored` (typically a collapsed `[H]`), making the
    // override invisible to the fingerprint.
    // -----------------------------------------------------------------

    fn assert_count_ecfp_differs_after_isotope_mutation(smiles: &str, seed: u64) {
        use crate::CountEcfpFingerprint;
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("fuzz-regression SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let inner = scratch.prepare(&parsed);
        let baseline = CountEcfpFingerprint::new(2, 65_536).compute(&inner);

        let mut wrapper = InvalidatedGraph::new(inner);
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        ImpossibleIsotopeMutator
            .mutate_in_place(&mut wrapper, &mut rng)
            .expect("mutator returned Err on a fuzz-regression input");
        let mutated_fp = CountEcfpFingerprint::new(2, 65_536).compute(&wrapper);
        assert_ne!(
            baseline, mutated_fp,
            "ImpossibleIsotopeMutator: ECFP unchanged for {smiles:?} / seed {seed}",
        );
    }

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
        let err = ImpossibleIsotopeMutator
            .mutate_in_place(&mut all_ignored, &mut rng)
            .expect_err("every atom is ignored; mutator should decline");
        assert_eq!(err, MutatorError::NoEligibleAtom);
    }

    #[test]
    fn fuzz_regression_7c50603e_ignored_h_atom() {
        assert_count_ecfp_differs_after_isotope_mutation("[H]#N7$cP7", 6_583_512_081_663_346_294);
    }

    #[test]
    fn fuzz_regression_a0439b0c_ignored_h_atom() {
        assert_count_ecfp_differs_after_isotope_mutation("[H]#c", 3_026_418_949_592_973_311);
    }
}
