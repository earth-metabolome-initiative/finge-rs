//! Property-based tests for the four correctness invariants spelled out in
//! the plan's *Correctness invariants* section.
//!
//! 1. `prop_ecfp_changes_for_every_mutator` — for every mutator `M`, every
//!    SMILES from a small canonical corpus, and every successful
//!    `M.mutate(G, rng)`, the ECFP bit set at radius 2 must differ from the
//!    baseline (hash level).
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
//! Each property iterates the eight built-in mutators per generated case so
//! the corpus is exercised exhaustively.

#![cfg(test)]

use alloc::vec::Vec;

use proptest::{prelude::*, test_runner::TestCaseError};
use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};

use crate::{
    CountEcfpFingerprint, EcfpGraph, Fingerprint, HypervalentMutator, HypervalentPredicate,
    ImpossibleAtomicNumberMutator, ImpossibleAtomicNumberPredicate, ImpossibleBondTypeMutator,
    ImpossibleBondTypePredicate, ImpossibleChargeMutator, ImpossibleChargePredicate,
    ImpossibleHCountMutator, ImpossibleHCountPredicate, ImpossibleIsotopeMutator,
    ImpossibleIsotopePredicate, ImpossibleRingFlagMutator, ImpossibleRingFlagPredicate,
    InvalidatedGraph, MolecularAtom, Mutator, TopologicalPathologyMutator,
    TopologicalPathologyPredicate, ViolationPredicate,
    smiles_support_impl::{SmilesRdkitGraph, SmilesRdkitScratch},
};

/// Fingerprint configuration used by the ECFP-changes property.
///
/// Counts (not just presence) and a large fp_size make fold-collision
/// masking astronomically unlikely: every changed feature hash either
/// increments a new slot or changes the count at an existing one. A pure
/// bit-set comparison at fp_size = 2048 was *not* a sound invariant — for
/// CC(=O)O with seed 14865529113730674949 the bit set after an
/// `ImpossibleHCountMutator` mutation happened to coincide with the
/// baseline, even though the underlying field tuple had clearly changed.
fn make_fingerprint() -> CountEcfpFingerprint {
    CountEcfpFingerprint::new(2, 65_536)
}

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

/// Asserts that, when `mutator` accepts `inner`, the count-ECFP at radius 2
/// differs from `baseline`. Counts (not just presence) eliminate the
/// fold-collision masking that bit fingerprints are vulnerable to.
fn check_ecfp_hash_changes<'a, M>(
    mutator: &M,
    label: &str,
    inner: SmilesRdkitGraph<'a>,
    baseline: &crate::CountFingerprint,
    seed: u64,
) -> Result<(), TestCaseError>
where
    M: Mutator<SmilesRdkitGraph<'a>>,
{
    let mut wrapper = InvalidatedGraph::new(inner);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    if mutator.mutate_in_place(&mut wrapper, &mut rng).is_ok() {
        let mutated_fp = make_fingerprint().compute(&wrapper);
        prop_assert_ne!(
            baseline,
            &mutated_fp,
            "{} produced unchanged count-ECFP at radius 2 (seed = {})",
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
        let baseline = make_fingerprint().compute(&inner);

        check_ecfp_hash_changes(&ImpossibleAtomicNumberMutator, "ImpossibleAtomicNumberMutator", inner, &baseline, seed)?;
        check_ecfp_hash_changes(&HypervalentMutator, "HypervalentMutator", inner, &baseline, seed)?;
        check_ecfp_hash_changes(&ImpossibleHCountMutator, "ImpossibleHCountMutator", inner, &baseline, seed)?;
        check_ecfp_hash_changes(&ImpossibleChargeMutator, "ImpossibleChargeMutator", inner, &baseline, seed)?;
        check_ecfp_hash_changes(&ImpossibleIsotopeMutator, "ImpossibleIsotopeMutator", inner, &baseline, seed)?;
        check_ecfp_hash_changes(&ImpossibleRingFlagMutator, "ImpossibleRingFlagMutator", inner, &baseline, seed)?;
        check_ecfp_hash_changes(&ImpossibleBondTypeMutator, "ImpossibleBondTypeMutator", inner, &baseline, seed)?;
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
}
