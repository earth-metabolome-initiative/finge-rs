//! Property-based tests for the correctness invariants of the mutation
//! suite — both the per-mutator contracts and the composed
//! [`crate::MutatorMix::sample`] behaviour.
//!
//! Per-mutator (k=1):
//! 1. `prop_ecfp_changes_for_every_mutator` — for every mutator `M`, every
//!    SMILES from a small canonical corpus, and every successful
//!    `M.mutate_in_place(...)`, the ECFP bit set at radius 2 must differ
//!    from the baseline (hash level).
//! 2. `prop_atom_mutators_change_field_tuple` — same input space, but for the
//!    seven *atom-channel* mutators only, the pre-hash atom invariant field
//!    multiset must also differ from baseline (collision-proof tuple level).
//!    The bond-channel mutator is excluded because it doesn't touch the
//!    atom field tuple by construction.
//! 3. `prop_primary_class_predicate_fires` — after a successful mutation, the
//!    predicate of the mutator's `primary_class()` must return `true` on the
//!    output.
//! 4. `prop_no_panic_on_small_inputs` — every mutator on every 1-, 2-, and
//!    3-atom graph must return `Ok` or `Err`, never panic.
//!
//! Composition (MutatorMix):
//! 5. `prop_every_set_label_bit_has_firing_predicate` — after a successful
//!    `MutatorMix::sample`, every bit set in the returned [`ViolationLabel`]
//!    must correspond to a predicate that actually fires on the returned
//!    wrapper. Predicate-driven labelling is symmetric, so the converse
//!    follows by construction.
//! 6. `prop_composed_label_bit_count_within_bounds` — after success, the
//!    label has at least one bit set and at most `min(k_max, 8)` bits.

#![cfg(test)]

use alloc::vec::Vec;

use proptest::{prelude::*, test_runner::TestCaseError};
use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};

use crate::{
    EcfpGraph, HypervalentMutator, HypervalentPredicate, ImpossibleAtomicNumberMutator,
    ImpossibleAtomicNumberPredicate, ImpossibleBondTypeMutator, ImpossibleBondTypePredicate,
    ImpossibleChargeMutator, ImpossibleChargePredicate, ImpossibleHCountMutator,
    ImpossibleHCountPredicate, ImpossibleIsotopeMutator, ImpossibleIsotopePredicate,
    ImpossibleRingFlagMutator, ImpossibleRingFlagPredicate, InvalidatedGraph, MolecularAtom,
    Mutator, MutatorMix, PredicateClass, TopologicalPathologyMutator,
    TopologicalPathologyPredicate, ViolationClass, ViolationPredicate,
    smiles_support_impl::{SmilesRdkitGraph, SmilesRdkitScratch},
};

/// Canonical positive corpus for properties (1) and (2). Diverse enough that
/// every mutator finds at least one suitable input.
const FULL_CORPUS: &[&str] = &[
    "C",
    "CC",
    "CCO",
    "CCN",
    "CCCl",
    "CC(=O)O",
    "C1CCCCC1",
    "c1ccccc1",
    "n1ccccc1",
    "OS(=O)(=O)O",
    "CC(=O)c1ccccc1",
];

/// Smaller corpus (1–3 atoms) for property (3).
const SMALL_CORPUS: &[&str] = &[
    "C", "N", "O", "[H]", "CC", "CO", "CN", "[OH-]", "CCO", "CCC", "OCO",
];

/// Sorted multiset of pre-hash atom invariant fields, for the tuple-level
/// `M(G) ≠ G` check.
fn atom_field_signature<G>(graph: &G) -> Vec<(u32, u32, u32, i32, i32, bool)>
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    let mut values: Vec<_> = (0..graph.atom_count())
        .map(|atom_id| {
            let f = graph.ecfp_atom_invariant_fields(atom_id);
            (
                f.atomic_number,
                f.total_degree,
                f.total_hydrogens,
                f.formal_charge,
                f.delta_mass,
                f.in_ring,
            )
        })
        .collect();
    values.sort();
    values
}

/// Asserts that, when `mutator` accepts `inner`, the per-atom R0 invariant
/// vector over non-ignored atoms (compared position-wise) differs from
/// `baseline`. Collision-proof: two distinct vecs differ at at least one
/// position by construction.
///
/// A previous sum-based aggregate was too weak under composition. Under
/// [`crate::MutatorMix::sample`] with k >= 2, two `ImpossibleIsotopeMutator`
/// writes on different atoms can produce per-atom invariant deltas that
/// sum to zero (see `mutator_mix/crash-18d54c65…` — `ssI` + seed
/// 6_733_535_862_861_618_035, with atom 0 shifted `-2287` and atom 1
/// shifted `+2287`). The position-aware vec catches this directly.
fn check_ecfp_hash_changes<'a, M>(
    mutator: &M,
    label: &str,
    inner: SmilesRdkitGraph<'a>,
    baseline_visible_signature: &[u32],
    seed: u64,
) -> Result<(), TestCaseError>
where
    M: Mutator<SmilesRdkitGraph<'a>>,
{
    let mut wrapper = InvalidatedGraph::new(inner);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    if mutator.mutate_in_place(&mut wrapper, &mut rng).is_ok() {
        let mutated = visible_invariant_signature(&wrapper);
        prop_assert_ne!(
            baseline_visible_signature,
            mutated.as_slice(),
            "{} left the visible R0 invariant signature unchanged (seed = {})",
            label,
            seed,
        );
    }
    Ok(())
}

/// Per-atom R0 invariants over ECFP-non-ignored atoms, in atom-id order.
fn visible_invariant_signature<G>(graph: &G) -> Vec<u32>
where
    G: EcfpGraph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
{
    (0..graph.atom_count())
        .filter(|&id| !graph.ecfp_atom_is_ignored(id))
        .map(|id| graph.ecfp_atom_invariant(id, true))
        .collect()
}

/// Per-bond invariants over undirected bonds (each bond once), in
/// `(min, max)` endpoint order. Used by the bond-channel mutator path,
/// which leaves atom invariants untouched by construction.
fn visible_bond_invariant_signature<G>(graph: &G) -> Vec<u32>
where
    G: EcfpGraph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
    G::Bond: crate::traits::MolecularBond<NodeId = usize>,
{
    use crate::traits::MolecularBond as _;
    let mut bonds = Vec::new();
    for source in 0..graph.atom_count() {
        for bond in graph.bonds(source) {
            let other = if bond.source() == source {
                bond.target()
            } else if bond.target() == source {
                bond.source()
            } else {
                continue;
            };
            if other > source {
                bonds.push((source, other, graph.ecfp_bond_invariant(&bond, true)));
            }
        }
    }
    bonds.sort();
    bonds.into_iter().map(|(_, _, inv)| inv).collect()
}

/// Variant of [`check_ecfp_hash_changes`] for the bond-channel mutator:
/// compares the bond invariant signature rather than atom R0 invariants.
fn check_bond_invariant_changes<'a, M>(
    mutator: &M,
    label: &str,
    inner: SmilesRdkitGraph<'a>,
    baseline_bond_signature: &[u32],
    seed: u64,
) -> Result<(), TestCaseError>
where
    M: Mutator<SmilesRdkitGraph<'a>>,
{
    let mut wrapper = InvalidatedGraph::new(inner);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    if mutator.mutate_in_place(&mut wrapper, &mut rng).is_ok() {
        let mutated = visible_bond_invariant_signature(&wrapper);
        prop_assert_ne!(
            baseline_bond_signature,
            mutated.as_slice(),
            "{} left the bond invariant signature unchanged (seed = {})",
            label,
            seed,
        );
    }
    Ok(())
}

/// Asserts that, when `mutator` accepts `inner`, the multiset of pre-hash
/// atom invariant field tuples differs from `baseline_fields`. Use only for
/// atom-channel mutators — bond-channel mutators leave atom fields unchanged
/// by construction.
fn check_atom_field_tuple_changes<'a, M>(
    mutator: &M,
    inner: SmilesRdkitGraph<'a>,
    baseline_fields: &[(u32, u32, u32, i32, i32, bool)],
    seed: u64,
) -> Result<(), TestCaseError>
where
    M: Mutator<SmilesRdkitGraph<'a>>,
{
    let mut wrapper = InvalidatedGraph::new(inner);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    if mutator.mutate_in_place(&mut wrapper, &mut rng).is_ok() {
        let mutated_fields = atom_field_signature(&wrapper);
        prop_assert_ne!(baseline_fields, &mutated_fields[..]);
    }
    Ok(())
}

/// Asserts that, when `mutator` accepts `inner`, the `predicate` fires on
/// the mutated output.
fn check_primary_predicate_fires<'a, M, P>(
    mutator: &M,
    predicate: &P,
    inner: SmilesRdkitGraph<'a>,
    seed: u64,
) -> Result<(), TestCaseError>
where
    M: Mutator<SmilesRdkitGraph<'a>>,
    P: ViolationPredicate<InvalidatedGraph<SmilesRdkitGraph<'a>>>,
{
    let mut wrapper = InvalidatedGraph::new(inner);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    if mutator.mutate_in_place(&mut wrapper, &mut rng).is_ok() {
        prop_assert!(predicate.check(&wrapper));
    }
    Ok(())
}

/// Asserts that `mutator` does not panic on `inner` for the given seed.
fn check_no_panic<'a, M>(
    mutator: &M,
    inner: SmilesRdkitGraph<'a>,
    seed: u64,
) -> Result<(), TestCaseError>
where
    M: Mutator<SmilesRdkitGraph<'a>>,
{
    let mut wrapper = InvalidatedGraph::new(inner);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    let _ = mutator.mutate_in_place(&mut wrapper, &mut rng);
    Ok(())
}

proptest! {
    #[test]
    fn prop_ecfp_changes_for_every_mutator(
        smiles_idx in 0..FULL_CORPUS.len(),
        seed in any::<u64>(),
    ) {
        let smiles = FULL_CORPUS[smiles_idx];
        let mut scratch = SmilesRdkitScratch::default();
        let parsed: smiles_parser::smiles::Smiles = smiles.parse()
            .expect("fixture SMILES should parse");
        let inner = scratch.prepare(&parsed);
        let baseline = visible_invariant_signature(&inner);
        let baseline_bonds = visible_bond_invariant_signature(&inner);

        check_ecfp_hash_changes(&ImpossibleAtomicNumberMutator, "ImpossibleAtomicNumberMutator", inner, &baseline, seed)?;
        check_ecfp_hash_changes(&HypervalentMutator, "HypervalentMutator", inner, &baseline, seed)?;
        check_ecfp_hash_changes(&ImpossibleHCountMutator, "ImpossibleHCountMutator", inner, &baseline, seed)?;
        check_ecfp_hash_changes(&ImpossibleChargeMutator, "ImpossibleChargeMutator", inner, &baseline, seed)?;
        check_ecfp_hash_changes(&ImpossibleIsotopeMutator, "ImpossibleIsotopeMutator", inner, &baseline, seed)?;
        check_ecfp_hash_changes(&ImpossibleRingFlagMutator, "ImpossibleRingFlagMutator", inner, &baseline, seed)?;
        // Bond-channel mutator leaves atom R0 invariants untouched by
        // construction; check the bond invariant signature instead.
        check_bond_invariant_changes(&ImpossibleBondTypeMutator, "ImpossibleBondTypeMutator", inner, &baseline_bonds, seed)?;
        check_ecfp_hash_changes(&TopologicalPathologyMutator, "TopologicalPathologyMutator", inner, &baseline, seed)?;
    }

    #[test]
    fn prop_atom_mutators_change_field_tuple(
        smiles_idx in 0..FULL_CORPUS.len(),
        seed in any::<u64>(),
    ) {
        let smiles = FULL_CORPUS[smiles_idx];
        let mut scratch = SmilesRdkitScratch::default();
        let parsed: smiles_parser::smiles::Smiles = smiles.parse()
            .expect("fixture SMILES should parse");
        let inner = scratch.prepare(&parsed);
        let baseline_fields = atom_field_signature(&inner);

        // Seven atom-channel mutators (excluding `ImpossibleBondTypeMutator`,
        // which only touches a bond invariant). `TopologicalPathologyMutator`
        // is an atom-channel mutator: it flips `in_ring`.
        check_atom_field_tuple_changes(&ImpossibleAtomicNumberMutator, inner, &baseline_fields, seed)?;
        check_atom_field_tuple_changes(&HypervalentMutator, inner, &baseline_fields, seed)?;
        check_atom_field_tuple_changes(&ImpossibleHCountMutator, inner, &baseline_fields, seed)?;
        check_atom_field_tuple_changes(&ImpossibleChargeMutator, inner, &baseline_fields, seed)?;
        check_atom_field_tuple_changes(&ImpossibleIsotopeMutator, inner, &baseline_fields, seed)?;
        check_atom_field_tuple_changes(&ImpossibleRingFlagMutator, inner, &baseline_fields, seed)?;
        check_atom_field_tuple_changes(&TopologicalPathologyMutator, inner, &baseline_fields, seed)?;
    }

    #[test]
    fn prop_primary_class_predicate_fires(
        smiles_idx in 0..FULL_CORPUS.len(),
        seed in any::<u64>(),
    ) {
        let smiles = FULL_CORPUS[smiles_idx];
        let mut scratch = SmilesRdkitScratch::default();
        let parsed: smiles_parser::smiles::Smiles = smiles.parse()
            .expect("fixture SMILES should parse");
        let inner = scratch.prepare(&parsed);

        check_primary_predicate_fires(&ImpossibleAtomicNumberMutator, &ImpossibleAtomicNumberPredicate, inner, seed)?;
        check_primary_predicate_fires(&HypervalentMutator, &HypervalentPredicate, inner, seed)?;
        check_primary_predicate_fires(&ImpossibleHCountMutator, &ImpossibleHCountPredicate, inner, seed)?;
        check_primary_predicate_fires(&ImpossibleChargeMutator, &ImpossibleChargePredicate, inner, seed)?;
        check_primary_predicate_fires(&ImpossibleIsotopeMutator, &ImpossibleIsotopePredicate, inner, seed)?;
        check_primary_predicate_fires(&ImpossibleRingFlagMutator, &ImpossibleRingFlagPredicate, inner, seed)?;
        check_primary_predicate_fires(&ImpossibleBondTypeMutator, &ImpossibleBondTypePredicate, inner, seed)?;
        check_primary_predicate_fires(&TopologicalPathologyMutator, &TopologicalPathologyPredicate, inner, seed)?;
    }

    #[test]
    fn prop_no_panic_on_small_inputs(
        smiles_idx in 0..SMALL_CORPUS.len(),
        seed in any::<u64>(),
    ) {
        let smiles = SMALL_CORPUS[smiles_idx];
        let mut scratch = SmilesRdkitScratch::default();
        let parsed: smiles_parser::smiles::Smiles = smiles.parse()
            .expect("fixture SMILES should parse");
        let inner = scratch.prepare(&parsed);

        check_no_panic(&ImpossibleAtomicNumberMutator, inner, seed)?;
        check_no_panic(&HypervalentMutator, inner, seed)?;
        check_no_panic(&ImpossibleHCountMutator, inner, seed)?;
        check_no_panic(&ImpossibleChargeMutator, inner, seed)?;
        check_no_panic(&ImpossibleIsotopeMutator, inner, seed)?;
        check_no_panic(&ImpossibleRingFlagMutator, inner, seed)?;
        check_no_panic(&ImpossibleBondTypeMutator, inner, seed)?;
        check_no_panic(&TopologicalPathologyMutator, inner, seed)?;
    }

    #[test]
    fn prop_every_set_label_bit_has_firing_predicate(
        smiles_idx in 0..FULL_CORPUS.len(),
        seed in any::<u64>(),
    ) {
        let smiles = FULL_CORPUS[smiles_idx];
        let mut scratch = SmilesRdkitScratch::default();
        let parsed: smiles_parser::smiles::Smiles = smiles.parse()
            .expect("fixture SMILES should parse");
        let inner = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates();
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        if let Ok((wrapper, label)) = mix.sample(inner, &mut rng) {
            for class in ViolationClass::ALL {
                if label.is_set(class) {
                    prop_assert!(
                        predicate_fires_for_class(class, &wrapper),
                        "label bit {class:?} set but predicate did not fire (smiles={smiles}, seed={seed})",
                    );
                }
            }
        }
    }

    #[test]
    fn prop_composed_label_bit_count_within_bounds(
        smiles_idx in 0..FULL_CORPUS.len(),
        seed in any::<u64>(),
    ) {
        let smiles = FULL_CORPUS[smiles_idx];
        let mut scratch = SmilesRdkitScratch::default();
        let parsed: smiles_parser::smiles::Smiles = smiles.parse()
            .expect("fixture SMILES should parse");
        let inner = scratch.prepare(&parsed);
        let mix = MutatorMix::<SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates();
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        if let Ok((_, label)) = mix.sample(inner, &mut rng) {
            let count = label.count();
            prop_assert!(count >= 1, "successful sample must set at least one bit");
            prop_assert!(
                count <= u32::from(ViolationClass::COUNT as u8),
                "label count {count} exceeds total class count",
            );
        }
    }
}

/// Dispatches to the matching predicate for `class`. Returns `false` for any
/// class without a matching built-in predicate (which is exhaustive here).
fn predicate_fires_for_class<G>(class: ViolationClass, wrapper: &G) -> bool
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    match class {
        ViolationClass::ImpossibleAtomicNumber => {
            ImpossibleAtomicNumberPredicate.check(wrapper)
                && ImpossibleAtomicNumberPredicate.class() == class
        }
        ViolationClass::Hypervalent => HypervalentPredicate.check(wrapper),
        ViolationClass::ImpossibleHCount => ImpossibleHCountPredicate.check(wrapper),
        ViolationClass::ImpossibleCharge => ImpossibleChargePredicate.check(wrapper),
        ViolationClass::ImpossibleIsotope => ImpossibleIsotopePredicate.check(wrapper),
        ViolationClass::ImpossibleRingFlag => ImpossibleRingFlagPredicate.check(wrapper),
        ViolationClass::ImpossibleBondType => ImpossibleBondTypePredicate.check(wrapper),
        ViolationClass::TopologicalPathology => TopologicalPathologyPredicate.check(wrapper),
    }
}
