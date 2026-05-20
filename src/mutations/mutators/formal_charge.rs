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

const CHARGE_CHOICES: [i32; 6] = [-9, -8, -7, 7, 8, 9];

/// Mutator targeting [`ViolationClass::ImpossibleCharge`].
///
/// Picks one random atom and writes a `formal_charge` chosen uniformly from
/// `{-9, -8, -7, +7, +8, +9}`. All values have `|q| > 6` and so always trip
/// the matching predicate.
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleChargeMutator;

impl<G> Mutator<G> for ImpossibleChargeMutator
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    type Output = InvalidatedGraph<G>;

    fn mutate(&self, graph: G, rng: &mut dyn RngCore) -> Result<Self::Output, MutatorError> {
        let target = pick_unignored_atom(rng, &graph).ok_or(MutatorError::NoEligibleAtom)?;
        let current = graph.ecfp_atom_invariant_fields(target).formal_charge;
        let mut pick = (rng.next_u32() as usize) % CHARGE_CHOICES.len();
        // Guarantee `override != current`. The choice set is small enough
        // that incrementing `pick` mod len trivially finds a different
        // value (every other entry in the set has a distinct magnitude).
        if CHARGE_CHOICES[pick] == current {
            pick = (pick + 1) % CHARGE_CHOICES.len();
        }
        let override_q = CHARGE_CHOICES[pick];

        let mut wrapper = InvalidatedGraph::new(graph);
        wrapper.set_formal_charge_override(target, override_q);
        Ok(wrapper)
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
        EcfpFingerprint, Fingerprint, Mutator, ViolationClass,
        mutations::predicate::{
            ImpossibleChargePredicate, ViolationPredicate, has_impossible_charge,
        },
        smiles_support_impl::SmilesRdkitScratch,
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

    #[test]
    fn mutator_returns_ok_and_predicate_fires() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut rng = ChaCha8Rng::seed_from_u64(0x20);
        let mutated = ImpossibleChargeMutator
            .mutate(inner, &mut rng)
            .expect("mutation should succeed");
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
        use crate::{InvalidatedGraph, MutatorError, traits::MolecularGraph as _};
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut all_ignored = InvalidatedGraph::new(inner);
        for atom_id in 0..all_ignored.atom_count() {
            all_ignored.set_atom_is_ignored_override(atom_id, true);
        }
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let err = ImpossibleChargeMutator
            .mutate(all_ignored, &mut rng)
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
        let mut rng = ChaCha8Rng::seed_from_u64(5_327_177);
        let mutated = ImpossibleChargeMutator
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
        let mut rng = ChaCha8Rng::seed_from_u64(0x21);
        let mutated = ImpossibleChargeMutator
            .mutate(inner, &mut rng)
            .expect("mutation should succeed");
        let mutated_fp = EcfpFingerprint::new(2, 2048).compute(&mutated);
        assert_ne!(
            baseline.active_bits().collect::<Vec<_>>(),
            mutated_fp.active_bits().collect::<Vec<_>>(),
        );
    }
}
