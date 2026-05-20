//! [`HypervalentMutator`] — picks a random atom and overrides its
//! `total_degree` to a value past the element's natural valence ceiling.

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

/// Mutator targeting [`ViolationClass::Hypervalent`].
///
/// Overrides the chosen atom's `total_degree` to
/// `max_natural_valence(Z) + Δ` where `Δ ∈ {2, 3, 4, 5}` is drawn from `rng`.
#[derive(Clone, Copy, Debug, Default)]
pub struct HypervalentMutator;

impl<G> Mutator<G> for HypervalentMutator
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    type Output = InvalidatedGraph<G>;

    fn mutate(&self, graph: G, rng: &mut dyn RngCore) -> Result<Self::Output, MutatorError> {
        let target = pick_unignored_atom(rng, &graph).ok_or(MutatorError::NoEligibleAtom)?;
        let fields = graph.ecfp_atom_invariant_fields(target);
        let bump = 2 + (rng.next_u32() % 4);
        let mut override_degree = max_natural_valence(fields.atomic_number).saturating_add(bump);
        // Guarantee `override != current`. For weird inputs (e.g. argon with
        // `max_natural_valence(18) == 0` and a natural degree of 2) the
        // `max + bump` formula can coincide with the atom's existing degree,
        // turning the mutation into a silent no-op
        // (fuzz/artifacts/mutator_hypervalent/crash-0794abb9…, …/crash-6358b428…).
        if override_degree == fields.total_degree {
            override_degree = override_degree.saturating_add(1);
        }

        let mut wrapper = InvalidatedGraph::new(graph);
        wrapper.set_total_degree_override(target, override_degree);
        Ok(wrapper)
    }

    #[inline]
    fn primary_class(&self) -> ViolationClass {
        ViolationClass::Hypervalent
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};

    use super::HypervalentMutator;
    use crate::{
        EcfpFingerprint, Fingerprint, Mutator, ViolationClass,
        mutations::predicate::{HypervalentPredicate, ViolationPredicate, is_hypervalent},
        smiles_support_impl::SmilesRdkitScratch,
    };

    fn prepared(smiles: &str) -> (SmilesRdkitScratch, smiles_parser::smiles::Smiles) {
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("fixture SMILES should parse");
        (SmilesRdkitScratch::default(), parsed)
    }

    #[test]
    fn primary_class_is_hypervalent() {
        assert_eq!(
            <HypervalentMutator as Mutator<
                crate::smiles_support_impl::SmilesRdkitGraph<'_>,
            >>::primary_class(&HypervalentMutator),
            ViolationClass::Hypervalent,
        );
    }

    #[test]
    fn mutator_returns_ok_and_predicate_fires() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut rng = ChaCha8Rng::seed_from_u64(0x01);
        let mutated = HypervalentMutator
            .mutate(inner, &mut rng)
            .expect("mutation should succeed");
        assert!(HypervalentPredicate.check(&mutated));
        assert!(is_hypervalent(&mutated));
    }

    // -----------------------------------------------------------------
    // Fuzz regressions — see `fuzz/artifacts/mutator_hypervalent/`.
    //
    // The mutator computes `override_degree = max_natural_valence(Z) + bump`
    // where `bump ∈ {2,3,4,5}`. For pathological SMILES where the picked
    // atom's *current* total_degree happens to equal that value (e.g.
    // `b[Ar]b`: argon has degree 2, `max_natural_valence(18) == 0`, so
    // `override == 2` is a no-op when `bump == 2`), the override changes
    // nothing and ECFP stays byte-identical to the baseline.
    // -----------------------------------------------------------------

    fn assert_count_ecfp_differs_after_hypervalent_mutation(smiles: &str, seed: u64) {
        use crate::CountEcfpFingerprint;
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("fuzz-regression SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let inner = scratch.prepare(&parsed);
        let baseline = CountEcfpFingerprint::new(2, 65_536).compute(&inner);
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let mutated = HypervalentMutator
            .mutate(inner, &mut rng)
            .expect("mutator returned Err on a fuzz-regression input");
        let mutated_fp = CountEcfpFingerprint::new(2, 65_536).compute(&mutated);
        assert_ne!(
            baseline, mutated_fp,
            "HypervalentMutator: ECFP unchanged for {smiles:?} / seed {seed}",
        );
    }

    #[test]
    fn fuzz_regression_0794abb9_override_equals_current_degree() {
        assert_count_ecfp_differs_after_hypervalent_mutation("b[Ar]b", 8_608_480_054_637_065_171);
    }

    #[test]
    fn fuzz_regression_6358b428_override_equals_current_degree() {
        assert_count_ecfp_differs_after_hypervalent_mutation("b[Ar]b", 7_094_383_674_944_472_998);
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
        let err = HypervalentMutator
            .mutate(all_ignored, &mut rng)
            .expect_err("every atom is ignored; mutator should decline");
        assert_eq!(err, MutatorError::NoEligibleAtom);
    }

    #[test]
    fn ecfp_changes_at_radius_two() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline = EcfpFingerprint::new(2, 2048).compute(&inner);
        let mut rng = ChaCha8Rng::seed_from_u64(0x02);
        let mutated = HypervalentMutator
            .mutate(inner, &mut rng)
            .expect("mutation should succeed");
        let mutated_fp = EcfpFingerprint::new(2, 2048).compute(&mutated);
        assert_ne!(
            baseline.active_bits().collect::<Vec<_>>(),
            mutated_fp.active_bits().collect::<Vec<_>>(),
        );
    }
}
