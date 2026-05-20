//! [`ImpossibleAtomicNumberMutator`] — overrides one random atom's
//! `atomic_number` to `0` or a super-heavy value outside the periodic table.
//!
//! The chosen target ECFP-atom-invariant tuple lives in the same field space
//! as a real atom — only the `atomic_number` channel violates chemistry. The
//! standard `AtomInvariantFields::to_invariant` hash then produces a 32-bit
//! invariant indistinguishable in shape from one a "real Z = 200 carbon"
//! would emit.

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

/// Range of "super-heavy" atomic numbers the mutator picks from (inclusive).
const SUPER_HEAVY_Z_LOW: u32 = 119;
const SUPER_HEAVY_Z_HIGH: u32 = 200;

/// Mutator targeting [`ViolationClass::ImpossibleAtomicNumber`].
///
/// Picks one random atom and overrides its `atomic_number` to either `0`
/// (wildcard / no element) or a random `Z` in `119..=200` (super-heavy /
/// undiscovered).
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleAtomicNumberMutator;

impl<G> Mutator<G> for ImpossibleAtomicNumberMutator
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    type Output = InvalidatedGraph<G>;

    fn mutate(&self, graph: G, rng: &mut dyn RngCore) -> Result<Self::Output, MutatorError> {
        let target = pick_unignored_atom(rng, &graph).ok_or(MutatorError::NoEligibleAtom)?;

        let choice = rng.next_u32();
        let span = SUPER_HEAVY_Z_HIGH - SUPER_HEAVY_Z_LOW + 1;
        let override_z = if choice % 2 == 0 {
            0
        } else {
            SUPER_HEAVY_Z_LOW + (choice / 2) % span
        };

        // Defensive: if the inner happens to already carry a violating Z (it
        // shouldn't for any well-behaved `EcfpGraph` impl), and we'd write the
        // same value, the mutation would be a no-op — bump to the alternate.
        let current_z = graph.ecfp_atom_invariant_fields(target).atomic_number;
        let override_z = if override_z == current_z {
            if override_z == 0 {
                SUPER_HEAVY_Z_LOW
            } else {
                0
            }
        } else {
            override_z
        };

        let mut wrapper = InvalidatedGraph::new(graph);
        wrapper.set_atomic_number_override(target, override_z);
        Ok(wrapper)
    }

    #[inline]
    fn primary_class(&self) -> ViolationClass {
        ViolationClass::ImpossibleAtomicNumber
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};

    use super::ImpossibleAtomicNumberMutator;
    use crate::{
        AtomInvariantFields, EcfpFingerprint, EcfpGraph, Fingerprint, MolecularAtom, Mutator,
        ViolationClass,
        mutations::predicate::{ImpossibleAtomicNumberPredicate, ViolationPredicate},
        smiles_support_impl::SmilesRdkitScratch,
        traits::MolecularGraph as _,
    };

    /// Per-mutator invariance check pattern.
    ///
    /// Returns the sorted multiset of pre-hash atom invariant field tuples
    /// exposed by `graph`. Used to assert that `M(G)` differs from `G` at
    /// the tuple level (collision-proof, independent of folding).
    fn atom_field_multiset<G>(graph: &G) -> Vec<AtomInvariantFields>
    where
        G: EcfpGraph,
        G::NodeSymbol: MolecularAtom,
    {
        let mut values: Vec<AtomInvariantFields> = (0..graph.atom_count())
            .map(|atom_id| graph.ecfp_atom_invariant_fields(atom_id))
            .collect();
        values.sort_by_key(|f| {
            (
                f.atomic_number,
                f.total_degree,
                f.total_hydrogens,
                f.formal_charge,
                f.delta_mass,
                f.in_ring,
            )
        });
        values
    }

    fn prepared(smiles: &str) -> (SmilesRdkitScratch, smiles_parser::smiles::Smiles) {
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("fixture SMILES should parse");
        (SmilesRdkitScratch::default(), parsed)
    }

    #[test]
    fn mutator_reports_correct_primary_class() {
        let mutator = ImpossibleAtomicNumberMutator;
        assert_eq!(
            <ImpossibleAtomicNumberMutator as Mutator<
                crate::smiles_support_impl::SmilesRdkitGraph<'_>,
            >>::primary_class(&mutator),
            ViolationClass::ImpossibleAtomicNumber,
        );
    }

    #[test]
    fn mutator_returns_ok_on_real_molecule() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut rng = ChaCha8Rng::seed_from_u64(0xCAFE_F00D);
        let result = ImpossibleAtomicNumberMutator.mutate(inner, &mut rng);
        assert!(result.is_ok());
    }

    #[test]
    fn primary_predicate_fires_on_mutated_output() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut rng = ChaCha8Rng::seed_from_u64(0x1234_5678);
        let mutated = ImpossibleAtomicNumberMutator
            .mutate(inner, &mut rng)
            .expect("mutation should succeed");
        let predicate = ImpossibleAtomicNumberPredicate;
        assert!(predicate.check(&mutated));
    }

    #[test]
    fn ecfp_changes_at_radius_two_hash_level() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline = EcfpFingerprint::new(2, 2048).compute(&inner);

        let mut rng = ChaCha8Rng::seed_from_u64(0xDEAD_BEEF);
        let mutated = ImpossibleAtomicNumberMutator
            .mutate(inner, &mut rng)
            .expect("mutation should succeed");
        let mutated_fp = EcfpFingerprint::new(2, 2048).compute(&mutated);

        assert_ne!(
            baseline.active_bits().collect::<Vec<_>>(),
            mutated_fp.active_bits().collect::<Vec<_>>(),
        );
    }

    #[test]
    fn atom_field_multiset_differs_tuple_level() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline_multiset = atom_field_multiset(&inner);

        let mut rng = ChaCha8Rng::seed_from_u64(0xABCD_EF01);
        let mutated = ImpossibleAtomicNumberMutator
            .mutate(inner, &mut rng)
            .expect("mutation should succeed");
        let mutated_multiset = atom_field_multiset(&mutated);

        assert_ne!(baseline_multiset, mutated_multiset);
    }

    #[test]
    fn override_z_is_either_zero_or_super_heavy() {
        // With many seeds the override should be either 0 or in 119..=200.
        // This guards against off-by-one errors in the bit-twiddling above.
        for seed in 0..32u64 {
            let (mut scratch, parsed) = prepared("CCO");
            let inner = scratch.prepare(&parsed);
            let mut rng = ChaCha8Rng::seed_from_u64(seed);
            let mutated = ImpossibleAtomicNumberMutator
                .mutate(inner, &mut rng)
                .expect("mutation should succeed");

            let mut found_violation = false;
            for atom_id in 0..mutated.atom_count() {
                let z = mutated.ecfp_atom_invariant_fields(atom_id).atomic_number;
                if z == 0 || (119..=200).contains(&z) {
                    found_violation = true;
                    break;
                }
            }
            assert!(
                found_violation,
                "seed {seed} did not produce a Z=0 or Z in 119..=200"
            );
        }
    }

    // -----------------------------------------------------------------
    // Fuzz regressions — see `fuzz/artifacts/mutator_atomic_number/`.
    //
    // Before the fix, the mutator picked any target in `0..atom_count`,
    // failing to filter out atoms that `EcfpGraph::ecfp_atom_is_ignored`
    // skips during the radius-0 emission loop. For `SmilesRdkitGraph` with
    // `collapse_plain_explicit_hydrogens = true` (the default), an explicit
    // `[H]` is "ignored" — the mutator can land on it, write an
    // `atomic_number` override, and have ECFP simply skip the atom, leaving
    // both the fingerprint and the predicate output unchanged from the
    // baseline.
    // -----------------------------------------------------------------

    fn assert_count_ecfp_differs_after_atomic_number_mutation(smiles: &str, seed: u64) {
        use crate::CountEcfpFingerprint;
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("fuzz-regression SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let inner = scratch.prepare(&parsed);
        let baseline = CountEcfpFingerprint::new(2, 65_536).compute(&inner);

        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let mutated = ImpossibleAtomicNumberMutator
            .mutate(inner, &mut rng)
            .expect("mutator returned Err on a fuzz-regression input");
        let mutated_fp = CountEcfpFingerprint::new(2, 65_536).compute(&mutated);
        assert_ne!(
            baseline, mutated_fp,
            "ImpossibleAtomicNumberMutator: ECFP unchanged for {smiles:?} / seed {seed}",
        );
    }

    #[test]
    fn mutator_returns_no_eligible_atom_when_all_atoms_ignored() {
        use crate::{InvalidatedGraph, MutatorError};
        // Synthesise an "all atoms ignored" graph by wrapping CCO and
        // forcing every `ecfp_atom_is_ignored` override to true. The
        // outer mutator sees an `EcfpGraph` where no atom is eligible.
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut all_ignored = InvalidatedGraph::new(inner);
        for atom_id in 0..all_ignored.atom_count() {
            all_ignored.set_atom_is_ignored_override(atom_id, true);
        }
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let err = ImpossibleAtomicNumberMutator
            .mutate(all_ignored, &mut rng)
            .expect_err("every atom is ignored; mutator should decline");
        assert_eq!(err, MutatorError::NoEligibleAtom);
    }

    #[test]
    fn fuzz_regression_559b89bb_ignored_h_atom() {
        assert_count_ecfp_differs_after_atomic_number_mutation(
            "[I]S.Onnppppp[H]",
            6_027_178_212_275_312_294,
        );
    }

    #[test]
    fn fuzz_regression_75a084e2_ignored_h_atom() {
        assert_count_ecfp_differs_after_atomic_number_mutation(
            "[H]S.O[H]S",
            48_038_394_878_644_014,
        );
    }

    #[test]
    fn fuzz_regression_bd77b81f_ignored_h_atom() {
        assert_count_ecfp_differs_after_atomic_number_mutation("[H]S.O", 3_481_051_972_378_773_595);
    }

    #[test]
    fn fuzz_regression_fcec3d66_ignored_h_atom() {
        assert_count_ecfp_differs_after_atomic_number_mutation("[H]S.O", 8_353_159_651_921_816_812);
    }
}
