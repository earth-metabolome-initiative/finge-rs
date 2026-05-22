//! [`ImpossibleBondTypeMutator`] — overrides one random bond's invariant to
//! a value outside the RDKit-recognised set `{1, 2, 3, 12}`.

use alloc::vec::Vec;

use rand_core::RngCore;

use crate::{
    mutations::{
        invalidated_graph::InvalidatedGraph,
        mutator::{Mutator, MutatorError},
        mutators::pick_from_slice,
        violation_class::ViolationClass,
    },
    traits::{EcfpGraph, MolecularAtom, MolecularBond},
};

/// Collects every undirected bond whose *both* endpoints are unignored.
///
/// `EcfpGraph::ecfp_atom_is_ignored` atoms (e.g. collapsed `[H]`) are
/// excluded by the ECFP adjacency builder — any bond touching one is
/// silently dropped before the bond invariant is hashed. A bond-invariant
/// override on such a bond is therefore invisible to the resulting
/// fingerprint, which is exactly the
/// `fuzz/artifacts/mutator_bond_type/crash-38177a52…`,
/// `…/crash-6050a67c…` pattern (`[H@]F`).
fn collect_unique_bonds<G>(graph: &G) -> Vec<(usize, usize)>
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    let mut bonds = Vec::new();
    for source in 0..graph.atom_count() {
        if graph.ecfp_atom_is_ignored(source) {
            continue;
        }
        for bond in graph.bonds(source) {
            let other = if bond.source() == source {
                bond.target()
            } else if bond.target() == source {
                bond.source()
            } else {
                continue;
            };
            if other > source && !graph.ecfp_atom_is_ignored(other) {
                bonds.push((source, other));
            }
        }
    }
    bonds
}

/// Mutator targeting [`ViolationClass::ImpossibleBondType`].
///
/// Picks one random bond and writes a `u32` invariant in `13..=112`. Every
/// value in that range is outside the RDKit-recognised set `{1, 2, 3, 12}`,
/// so the matching predicate always fires.
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleBondTypeMutator;

impl<G> Mutator<G> for ImpossibleBondTypeMutator
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    fn mutate_in_place(
        &self,
        wrapper: &mut InvalidatedGraph<G>,
        rng: &mut dyn RngCore,
    ) -> Result<(), MutatorError> {
        let bonds = collect_unique_bonds(wrapper);
        let &(source, target) = pick_from_slice(rng, &bonds).ok_or(MutatorError::NoEligibleBond)?;
        let override_inv = 13 + (rng.next_u32() % 100);

        wrapper.set_bond_invariant_override(source, target, override_inv);
        Ok(())
    }

    #[inline]
    fn primary_class(&self) -> ViolationClass {
        ViolationClass::ImpossibleBondType
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};

    use super::ImpossibleBondTypeMutator;
    use crate::{
        EcfpFingerprint, Fingerprint, InvalidatedGraph, Mutator, ViolationClass,
        mutations::predicate::{
            ImpossibleBondTypePredicate, ViolationPredicate, has_impossible_bond_type,
        },
        smiles_support_impl::SmilesRdkitScratch,
    };

    fn prepared(smiles: &str) -> (SmilesRdkitScratch, smiles_parser::smiles::Smiles) {
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("fixture SMILES should parse");
        (SmilesRdkitScratch::default(), parsed)
    }

    #[test]
    fn primary_class_is_impossible_bond_type() {
        assert_eq!(
            <ImpossibleBondTypeMutator as Mutator<
                crate::smiles_support_impl::SmilesRdkitGraph<'_>,
            >>::primary_class(&ImpossibleBondTypeMutator),
            ViolationClass::ImpossibleBondType,
        );
    }

    fn mutate<'a>(
        inner: crate::smiles_support_impl::SmilesRdkitGraph<'a>,
        seed: u64,
    ) -> InvalidatedGraph<crate::smiles_support_impl::SmilesRdkitGraph<'a>> {
        let mut wrapper = InvalidatedGraph::new(inner);
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        ImpossibleBondTypeMutator
            .mutate_in_place(&mut wrapper, &mut rng)
            .expect("mutation should succeed");
        wrapper
    }

    #[test]
    fn mutator_returns_ok_and_predicate_fires() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mutated = mutate(inner, 0x50);
        assert!(ImpossibleBondTypePredicate.check(&mutated));
        assert!(has_impossible_bond_type(&mutated));
    }

    // -----------------------------------------------------------------
    // Fuzz regressions — see `fuzz/artifacts/mutator_bond_type/`.
    //
    // The mutator picks any unique bond, but the ECFP `adjacency()` loop in
    // `src/fingerprints/ecfp.rs` skips bonds where *either* endpoint is
    // `ecfp_atom_is_ignored`. SMILES with a collapsed explicit hydrogen
    // (`[H@]F`, `[H]F`, …) therefore have their only bond invisible to
    // ECFP, and a bond-invariant override is silently dropped.
    // -----------------------------------------------------------------

    fn assert_count_ecfp_differs_after_bond_type_mutation(smiles: &str, seed: u64) {
        use crate::CountEcfpFingerprint;
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("fuzz-regression SMILES should parse");
        let mut scratch = SmilesRdkitScratch::default();
        let inner = scratch.prepare(&parsed);
        let baseline = CountEcfpFingerprint::new(2, 65_536).compute(&inner);
        let mut wrapper = InvalidatedGraph::new(inner);
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        // The mutator may legitimately decline when every bond touches an
        // ignored atom; the invariant is "Ok ⇒ ECFP differs", not "always
        // Ok". The previous bug was that the mutator returned Ok but the
        // fingerprint was unchanged.
        if ImpossibleBondTypeMutator
            .mutate_in_place(&mut wrapper, &mut rng)
            .is_ok()
        {
            let mutated_fp = CountEcfpFingerprint::new(2, 65_536).compute(&wrapper);
            assert_ne!(
                baseline, mutated_fp,
                "ImpossibleBondTypeMutator: ECFP unchanged for {smiles:?} / seed {seed}",
            );
        }
    }

    #[test]
    fn mutator_returns_no_eligible_bond_when_every_bond_touches_ignored_atom() {
        use crate::{MutatorError, traits::MolecularGraph as _};
        // Wrap CCO and mark every atom as ignored. Every bond touches an
        // ignored endpoint, so the bond-channel filter discards every
        // candidate.
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        for atom_id in 0..wrapper.atom_count() {
            wrapper.set_atom_is_ignored_override(atom_id, true);
        }
        let mut rng = ChaCha8Rng::seed_from_u64(0);
        let err = ImpossibleBondTypeMutator
            .mutate_in_place(&mut wrapper, &mut rng)
            .expect_err("every atom is ignored; mutator should decline");
        assert_eq!(err, MutatorError::NoEligibleBond);
    }

    #[test]
    fn fuzz_regression_38177a52_bond_touches_ignored_atom() {
        assert_count_ecfp_differs_after_bond_type_mutation("[H@]F", u64::MAX);
    }

    #[test]
    fn fuzz_regression_6050a67c_bond_touches_ignored_atom() {
        assert_count_ecfp_differs_after_bond_type_mutation("[H@]F", 10_175_429_687_010_291_911);
    }

    #[test]
    fn ecfp_changes_at_radius_two() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline = EcfpFingerprint::new(2, 2048).compute(&inner);
        let mutated = mutate(inner, 0x51);
        let mutated_fp = EcfpFingerprint::new(2, 2048).compute(&mutated);
        assert_ne!(
            baseline.active_bits().collect::<Vec<_>>(),
            mutated_fp.active_bits().collect::<Vec<_>>(),
        );
    }
}
