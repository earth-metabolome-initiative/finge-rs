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

const DELTA_MASS_CHOICES: [i32; 6] = [-127, -120, -100, 100, 120, 127];

/// Mutator targeting [`ViolationClass::ImpossibleIsotope`].
///
/// Picks one random atom and writes a `delta_mass` chosen uniformly from
/// `{±127, ±120, ±100}`. All values have `|Δmass| > 80` and so always trip
/// the matching predicate.
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleIsotopeMutator;

impl<G> Mutator<G> for ImpossibleIsotopeMutator
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    type Output = InvalidatedGraph<G>;

    fn mutate(&self, graph: G, rng: &mut dyn RngCore) -> Result<Self::Output, MutatorError> {
        let target = pick_unignored_atom(rng, &graph).ok_or(MutatorError::NoEligibleAtom)?;
        let current = graph.ecfp_atom_invariant_fields(target).delta_mass;
        let mut pick = (rng.next_u32() as usize) % DELTA_MASS_CHOICES.len();
        // Guarantee `override != current` — same defensive scheme as
        // `ImpossibleChargeMutator`.
        if DELTA_MASS_CHOICES[pick] == current {
            pick = (pick + 1) % DELTA_MASS_CHOICES.len();
        }
        let override_dm = DELTA_MASS_CHOICES[pick];

        let mut wrapper = InvalidatedGraph::new(graph);
        wrapper.set_delta_mass_override(target, override_dm);
        Ok(wrapper)
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
        EcfpFingerprint, Fingerprint, Mutator, ViolationClass,
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

    #[test]
    fn mutator_returns_ok_and_predicate_fires() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut rng = ChaCha8Rng::seed_from_u64(0x30);
        let mutated = ImpossibleIsotopeMutator
            .mutate(inner, &mut rng)
            .expect("mutation should succeed");
        assert!(ImpossibleIsotopePredicate.check(&mutated));
        assert!(has_impossible_isotope(&mutated));
    }

    #[test]
    fn ecfp_changes_at_radius_two() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline = EcfpFingerprint::new(2, 2048).compute(&inner);
        let mut rng = ChaCha8Rng::seed_from_u64(0x31);
        let mutated = ImpossibleIsotopeMutator
            .mutate(inner, &mut rng)
            .expect("mutation should succeed");
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

        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let mutated = ImpossibleIsotopeMutator
            .mutate(inner, &mut rng)
            .expect("mutator returned Err on a fuzz-regression input");
        let mutated_fp = CountEcfpFingerprint::new(2, 65_536).compute(&mutated);
        assert_ne!(
            baseline, mutated_fp,
            "ImpossibleIsotopeMutator: ECFP unchanged for {smiles:?} / seed {seed}",
        );
    }

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
        let err = ImpossibleIsotopeMutator
            .mutate(all_ignored, &mut rng)
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
