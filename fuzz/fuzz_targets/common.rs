#![allow(dead_code)]

use finge_rs::{
    AtomPairFingerprint, BitFingerprint, CountEcfpFingerprint, CountFingerprint, EcfpFingerprint,
    EcfpGraph, Fingerprint, HypervalentMutator, HypervalentPredicate,
    ImpossibleAtomicNumberMutator, ImpossibleAtomicNumberPredicate, ImpossibleBondTypeMutator,
    ImpossibleBondTypePredicate, ImpossibleChargeMutator, ImpossibleChargePredicate,
    ImpossibleHCountMutator, ImpossibleHCountPredicate, ImpossibleIsotopeMutator,
    ImpossibleIsotopePredicate, ImpossibleRingFlagMutator, ImpossibleRingFlagPredicate,
    InvalidatedGraph, LayeredCountEcfpFingerprint, LayeredCountFingerprint, MaccsFingerprint,
    MolecularAtom, MolecularBond, MolecularGraph, Mutator, MutatorError,
    TopologicalPathologyMutator, TopologicalPathologyPredicate, TopologicalTorsionFingerprint,
    ViolationPredicate, max_natural_valence, smiles_support::SmilesRdkitScratch,
};
use rand_chacha::{ChaCha8Rng, rand_core::SeedableRng};
use smarts_rs::PreparedTarget;
use smiles_parser::smiles::Smiles;

const MAX_INPUT_BYTES: usize = 128;

/// Fingerprint pinned for every mutator harness. Used only for the per-seed
/// determinism check (re-running with the same seed must produce a byte-equal
/// fingerprint) — *not* for baseline-vs-mutated comparison, which collapses
/// to the same folded bin in roughly 1 in 65 536 mutations due to hash
/// collisions. The collision-proof check is the pre-hash field signature.
const MUTATOR_FP_RADIUS: u8 = 2;
const MUTATOR_FP_SIZE: usize = 65_536;

pub fn parse_smiles(input: String) -> Option<Smiles> {
    if input.is_empty() || input.len() > MAX_INPUT_BYTES {
        return None;
    }
    input.parse().ok()
}

pub fn prepare_graph<'a>(
    scratch: &'a mut SmilesRdkitScratch,
    smiles: &Smiles,
) -> Option<finge_rs::smiles_support::SmilesRdkitGraph<'a>> {
    scratch.try_prepare(smiles).ok()
}

pub fn assert_bit_fingerprint_basics(fingerprint: &BitFingerprint) {
    let active_bits = fingerprint.active_bits().collect::<Vec<_>>();
    assert!(active_bits.windows(2).all(|window| window[0] < window[1]));
    assert!(active_bits.iter().all(|&index| index < fingerprint.len()));
    assert!(active_bits.iter().all(|&index| fingerprint.contains(index)));
}

pub fn assert_count_fingerprint_basics(fingerprint: &CountFingerprint) {
    let active_counts = fingerprint.active_counts().collect::<Vec<_>>();
    assert!(
        active_counts
            .windows(2)
            .all(|window| window[0].0 < window[1].0)
    );
    assert!(
        active_counts
            .iter()
            .all(|&(index, count)| index < fingerprint.len() && count > 0)
    );
    for &(index, count) in &active_counts {
        assert_eq!(fingerprint.count(index), count);
    }
}

pub fn assert_layered_count_fingerprint_basics(fingerprint: &LayeredCountFingerprint) {
    assert!(!fingerprint.is_empty());
    assert_eq!(fingerprint.formula(), &fingerprint.layers()[0]);
    for layer in fingerprint.layers() {
        assert_count_fingerprint_basics(layer);
    }
}

pub fn assert_count_matches_bit_presence(counts: &CountFingerprint, bits: &BitFingerprint) {
    let active_bits = bits.active_bits().collect::<Vec<_>>();
    let active_count_bits = counts
        .active_counts()
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    assert_eq!(active_bits, active_count_bits);
}

pub fn assert_layered_sums_match_total(
    layered: &LayeredCountFingerprint,
    total: &CountFingerprint,
) {
    assert_eq!(layered.formula().len(), total.len());
    for index in 0..total.len() {
        let layered_sum = layered
            .layers()
            .iter()
            .map(|layer| layer.count(index))
            .sum::<u32>();
        assert_eq!(layered_sum, total.count(index));
    }
}

pub fn fuzz_ecfp_family(smiles: Smiles) {
    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };

    let bit = EcfpFingerprint::default().compute(&graph);
    let counted = CountEcfpFingerprint::default().compute(&graph);
    let layered = LayeredCountEcfpFingerprint::default().compute(&graph);
    let radius0 = CountEcfpFingerprint::new(0, 2048).compute(&graph);

    assert_bit_fingerprint_basics(&bit);
    assert_count_fingerprint_basics(&counted);
    assert_layered_count_fingerprint_basics(&layered);
    assert_count_matches_bit_presence(&counted, &bit);
    assert_layered_sums_match_total(&layered, &counted);
    assert_eq!(layered.formula(), &radius0);

    let bit_r1 = EcfpFingerprint::new(1, 2048).compute(&graph);
    let bit_r2 = EcfpFingerprint::new(2, 2048).compute(&graph);
    for index in bit_r1.active_bits() {
        assert!(bit_r2.contains(index));
    }

    let count_r1 = CountEcfpFingerprint::new(1, 2048).compute(&graph);
    let count_r2 = CountEcfpFingerprint::new(2, 2048).compute(&graph);
    for index in 0..count_r1.len() {
        assert!(count_r2.count(index) >= count_r1.count(index));
    }

    let explicit_h = smiles.with_explicit_hydrogens();
    let mut explicit_h_scratch = SmilesRdkitScratch::default();
    let Some(explicit_h_graph) = prepare_graph(&mut explicit_h_scratch, &explicit_h) else {
        return;
    };
    let explicit_h_bit = EcfpFingerprint::default().compute(&explicit_h_graph);
    let explicit_h_counted = CountEcfpFingerprint::default().compute(&explicit_h_graph);
    let explicit_h_layered = LayeredCountEcfpFingerprint::default().compute(&explicit_h_graph);

    assert_bit_fingerprint_basics(&explicit_h_bit);
    assert_count_fingerprint_basics(&explicit_h_counted);
    assert_layered_count_fingerprint_basics(&explicit_h_layered);
    assert_count_matches_bit_presence(&explicit_h_counted, &explicit_h_bit);
    assert_layered_sums_match_total(&explicit_h_layered, &explicit_h_counted);
}

pub fn fuzz_atom_pair(smiles: Smiles) {
    let raw_bits = AtomPairFingerprint::default().compute(&smiles);
    assert_bit_fingerprint_basics(&raw_bits);

    let without_count_sim = AtomPairFingerprint::default()
        .with_count_simulation(false)
        .compute(&smiles);
    assert_bit_fingerprint_basics(&without_count_sim);

    let explicit_h = smiles.with_explicit_hydrogens();
    let explicit_h_bits = AtomPairFingerprint::default().compute(&explicit_h);
    assert_bit_fingerprint_basics(&explicit_h_bits);

    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };
    let prepared_bits = AtomPairFingerprint::default().compute(&graph);
    let prepared_without_count_sim = AtomPairFingerprint::default()
        .with_count_simulation(false)
        .compute(&graph);

    assert_bit_fingerprint_basics(&prepared_bits);
    assert_bit_fingerprint_basics(&prepared_without_count_sim);
}

pub fn fuzz_topological_torsion(smiles: Smiles) {
    let default_bits = TopologicalTorsionFingerprint::default().compute(&smiles);
    let no_count_sim = TopologicalTorsionFingerprint::default()
        .with_count_simulation(false)
        .compute(&smiles);
    let shortest_paths = TopologicalTorsionFingerprint::default()
        .with_only_shortest_paths(true)
        .compute(&smiles);
    let torsion5 = TopologicalTorsionFingerprint::default()
        .with_torsion_atom_count(5)
        .compute(&smiles);

    assert_bit_fingerprint_basics(&default_bits);
    assert_bit_fingerprint_basics(&no_count_sim);
    assert_bit_fingerprint_basics(&shortest_paths);
    assert_bit_fingerprint_basics(&torsion5);

    let explicit_h = smiles.with_explicit_hydrogens();
    let explicit_h_bits = TopologicalTorsionFingerprint::default().compute(&explicit_h);
    assert_bit_fingerprint_basics(&explicit_h_bits);

    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };

    let prepared_default = TopologicalTorsionFingerprint::default().compute(&graph);
    let prepared_no_count_sim = TopologicalTorsionFingerprint::default()
        .with_count_simulation(false)
        .compute(&graph);
    let prepared_shortest_paths = TopologicalTorsionFingerprint::default()
        .with_only_shortest_paths(true)
        .compute(&graph);
    let prepared_torsion5 = TopologicalTorsionFingerprint::default()
        .with_torsion_atom_count(5)
        .compute(&graph);

    assert_bit_fingerprint_basics(&prepared_default);
    assert_bit_fingerprint_basics(&prepared_no_count_sim);
    assert_bit_fingerprint_basics(&prepared_shortest_paths);
    assert_bit_fingerprint_basics(&prepared_torsion5);

    for index in prepared_shortest_paths.active_bits() {
        assert!(prepared_default.contains(index));
    }
}

// ---------------------------------------------------------------------------
// Shared primitives for the per-mutator fuzz harnesses
// ---------------------------------------------------------------------------

/// Returns the count-ECFP every per-mutator harness pins to as its baseline:
/// radius 2, 65 536 slots, count-based (collision-proof under field changes).
pub fn count_ecfp_fingerprint<G>(graph: &G) -> CountFingerprint
where
    G: EcfpGraph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
    G::Bond: MolecularBond<NodeId = usize>,
{
    CountEcfpFingerprint::new(MUTATOR_FP_RADIUS, MUTATOR_FP_SIZE).compute(graph)
}

/// Sorted multiset of pre-hash atom invariant tuples — collision-proof
/// stand-in for the hash-level ECFP comparison, used to detect *any* field
/// change regardless of fold behaviour.
pub fn atom_field_signature<G>(graph: &G) -> Vec<(u32, u32, u32, i32, i32, bool)>
where
    G: EcfpGraph<NodeId = usize>,
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

/// Sorted list of `(min(src, dst), max(src, dst), ecfp_bond_invariant)`
/// triples — one entry per undirected bond. Used to detect bond-channel
/// changes in the bond-type mutator's minimality check.
pub fn bond_invariant_signature<G>(graph: &G) -> Vec<(usize, usize, u32)>
where
    G: EcfpGraph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
    G::Bond: MolecularBond<NodeId = usize>,
{
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
                let inv = graph.ecfp_bond_invariant(&bond, true);
                bonds.push((source, other, inv));
            }
        }
    }
    bonds.sort();
    bonds
}

/// Collects the radius-0 ECFP atom invariants of every *non-ignored* atom
/// in atom-id order. Comparing two such vecs is collision-free: if any
/// visible atom's R0 invariant differs, the vecs differ at that position.
///
/// Replaces a previous sum-based aggregate that was too weak under
/// composition: under [`crate::MutatorMix::sample`] with k >= 2, two
/// `ImpossibleIsotopeMutator` writes on different atoms can produce
/// per-atom invariant deltas that sum to zero (see
/// `mutator_mix/crash-18d54c65…` — `ssI` + seed
/// 6_733_535_862_861_618_035, where atom 0 shifted by `-2287` and atom 1
/// by `+2287`). The position-aware vec catches this directly.
///
/// Still catches the original "mutator wrote to an atom that ECFP doesn't
/// see" bug class: a write to an ignored atom is filtered out of the vec
/// and so leaves the signature unchanged.
pub fn visible_invariant_signature<G>(graph: &G) -> Vec<u32>
where
    G: EcfpGraph<NodeId = usize>,
    G::NodeSymbol: MolecularAtom,
{
    (0..graph.atom_count())
        .filter(|&id| !graph.ecfp_atom_is_ignored(id))
        .map(|id| graph.ecfp_atom_invariant(id, true))
        .collect()
}

/// Counts atom indices whose field tuples differ between two graphs of equal
/// atom count.
fn count_atom_field_diffs<G1, G2>(before: &G1, after: &G2) -> usize
where
    G1: EcfpGraph<NodeId = usize>,
    G1::NodeSymbol: MolecularAtom,
    G2: EcfpGraph<NodeId = usize>,
    G2::NodeSymbol: MolecularAtom,
{
    assert_eq!(
        before.atom_count(),
        after.atom_count(),
        "mutator must not change atom count",
    );
    (0..before.atom_count())
        .filter(|&atom_id| {
            before.ecfp_atom_invariant_fields(atom_id) != after.ecfp_atom_invariant_fields(atom_id)
        })
        .count()
}

/// Asserts that exactly `expected` atom-ids have differing field tuples
/// between `before` and `after`. The atom-channel mutators all pick exactly
/// one target, so `expected == 1` is the standard case.
pub fn assert_mutator_minimal_diff_atoms<G1, G2>(before: &G1, after: &G2, expected: usize)
where
    G1: EcfpGraph<NodeId = usize>,
    G1::NodeSymbol: MolecularAtom,
    G2: EcfpGraph<NodeId = usize>,
    G2::NodeSymbol: MolecularAtom,
{
    let diffs = count_atom_field_diffs(before, after);
    assert_eq!(
        diffs, expected,
        "expected exactly {expected} atom-id(s) to change, found {diffs}",
    );
}

/// Asserts that the atom field tuples of `before` and `after` are identical
/// everywhere. Used by the bond-channel mutator, which by design must leave
/// every atom's pre-hash invariant tuple untouched.
pub fn assert_mutator_no_atom_diff<G1, G2>(before: &G1, after: &G2)
where
    G1: EcfpGraph<NodeId = usize>,
    G1::NodeSymbol: MolecularAtom,
    G2: EcfpGraph<NodeId = usize>,
    G2::NodeSymbol: MolecularAtom,
{
    let diffs = count_atom_field_diffs(before, after);
    assert_eq!(
        diffs, 0,
        "bond-channel mutator changed {diffs} atom field tuple(s); expected 0",
    );
}

/// Asserts that exactly `expected` bonds have differing invariants between
/// `before` and `after`. Used by the bond-type mutator (`expected == 1`).
pub fn assert_mutator_minimal_diff_bonds<G1, G2>(before: &G1, after: &G2, expected: usize)
where
    G1: EcfpGraph<NodeId = usize>,
    G1::NodeSymbol: MolecularAtom,
    G1::Bond: MolecularBond<NodeId = usize>,
    G2: EcfpGraph<NodeId = usize>,
    G2::NodeSymbol: MolecularAtom,
    G2::Bond: MolecularBond<NodeId = usize>,
{
    let before_bonds = bond_invariant_signature(before);
    let after_bonds = bond_invariant_signature(after);
    assert_eq!(
        before_bonds.len(),
        after_bonds.len(),
        "mutator must not change bond count",
    );
    let diffs = before_bonds
        .iter()
        .zip(after_bonds.iter())
        .filter(|(b, a)| b.0 != a.0 || b.1 != a.1 || b.2 != a.2)
        .count();
    assert_eq!(
        diffs, expected,
        "expected exactly {expected} bond(s) to change, found {diffs}",
    );
}

/// Exhaustive match on `MutatorError` — the compiler will flag any newly
/// added variant. Used by every per-mutator harness on the `Err` branch as a
/// "well-defined error" sanity check.
pub fn assert_well_defined_error(err: MutatorError) {
    match err {
        MutatorError::GraphTooSmall
        | MutatorError::NoEligibleAtom
        | MutatorError::NoEligibleBond
        | MutatorError::AlreadyViolated
        | MutatorError::CompositionCancelled
        | MutatorError::FingerprintCollision => {}
    }
}

/// Enforces every invariant in the per-mutator contract for
/// `ImpossibleAtomicNumberMutator`.
pub fn fuzz_impossible_atomic_number_invariants(smiles: Smiles, seed: u64) {
    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };
    let baseline_atoms = atom_field_signature(&graph);

    let mut wrapper = InvalidatedGraph::new(graph);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    match ImpossibleAtomicNumberMutator.mutate_in_place(&mut wrapper, &mut rng) {
        Err(err) => assert_well_defined_error(err),
        Ok(()) => {
            // (3) hash-level visibility
            // (3') collision-proof ECFP visibility: the sum of R0
            // invariants over non-ignored atoms must change. Replaces the
            // previous "count-ECFP unchanged" check, which folded to
            // 65 536 bins and produced false positives on low-bit hash
            // collisions (see `B[ClH]` / `[84Ga][61Cu]` / `BrOBr` / `C`
            // fuzz crashes — atom 2's R0 in BrOBr coincidentally folded to
            // the same bin as atom 0's). This unfolded 32-bit sum still
            // catches the real bug class: a mutator that writes to an
            // atom ECFP cannot see (e.g. one skipped by
            // `pick_unignored_atom`).
            assert_ne!(
                visible_invariant_signature(&graph),
                visible_invariant_signature(&wrapper),
                "ImpossibleAtomicNumberMutator: visible R0 invariant signature unchanged \
                 — likely wrote to an ECFP-ignored atom",
            );

            // (4) primary predicate fires
            assert!(
                ImpossibleAtomicNumberPredicate.check(&wrapper),
                "ImpossibleAtomicNumberPredicate did not fire",
            );

            // (5) determinism — re-running with the same seed yields the
            // same ECFP
            let mut wrapper2 = InvalidatedGraph::new(graph);
            let mut rng2 = ChaCha8Rng::seed_from_u64(seed);
            ImpossibleAtomicNumberMutator
                .mutate_in_place(&mut wrapper2, &mut rng2)
                .expect("deterministic re-run must succeed");
            assert_eq!(
                count_ecfp_fingerprint(&wrapper),
                count_ecfp_fingerprint(&wrapper2),
                "ImpossibleAtomicNumberMutator is non-deterministic for a fixed seed",
            );

            // (6) minimality — exactly one atom changed
            assert_mutator_minimal_diff_atoms(&graph, &wrapper, 1);

            // (7) atom-channel: field multiset differs
            assert_ne!(
                baseline_atoms,
                atom_field_signature(&wrapper),
                "ImpossibleAtomicNumberMutator: atom field multiset unchanged",
            );

            // (8) class-specific minimum-violation invariant
            let z_violation = (0..wrapper.atom_count()).any(|atom_id| {
                let z = wrapper.ecfp_atom_invariant_fields(atom_id).atomic_number;
                z == 0 || z > 118
            });
            assert!(
                z_violation,
                "ImpossibleAtomicNumberMutator did not produce a Z outside 1..=118",
            );
        }
    }
}

/// Enforces every invariant in the per-mutator contract for
/// `HypervalentMutator`.
pub fn fuzz_hypervalent_invariants(smiles: Smiles, seed: u64) {
    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };
    let baseline_atoms = atom_field_signature(&graph);

    let mut wrapper = InvalidatedGraph::new(graph);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    match HypervalentMutator.mutate_in_place(&mut wrapper, &mut rng) {
        Err(err) => assert_well_defined_error(err),
        Ok(()) => {
            assert_ne!(
                visible_invariant_signature(&graph),
                visible_invariant_signature(&wrapper),
                "HypervalentMutator: visible R0 invariant signature unchanged — see \
                 `fuzz_impossible_atomic_number_invariants` for the rationale",
            );
            assert!(
                HypervalentPredicate.check(&wrapper),
                "HypervalentPredicate did not fire",
            );

            let mut wrapper2 = InvalidatedGraph::new(graph);
            let mut rng2 = ChaCha8Rng::seed_from_u64(seed);
            HypervalentMutator
                .mutate_in_place(&mut wrapper2, &mut rng2)
                .expect("deterministic re-run must succeed");
            assert_eq!(
                count_ecfp_fingerprint(&wrapper),
                count_ecfp_fingerprint(&wrapper2),
                "HypervalentMutator is non-deterministic for a fixed seed",
            );

            assert_mutator_minimal_diff_atoms(&graph, &wrapper, 1);
            assert_ne!(
                baseline_atoms,
                atom_field_signature(&wrapper),
                "HypervalentMutator: atom field multiset unchanged",
            );

            let degree_violation = (0..wrapper.atom_count()).any(|atom_id| {
                let f = wrapper.ecfp_atom_invariant_fields(atom_id);
                f.total_degree > max_natural_valence(f.atomic_number)
            });
            assert!(
                degree_violation,
                "HypervalentMutator did not produce total_degree > max_natural_valence(Z)",
            );
        }
    }
}

/// Enforces every invariant in the per-mutator contract for
/// `ImpossibleHCountMutator`.
pub fn fuzz_impossible_h_count_invariants(smiles: Smiles, seed: u64) {
    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };
    let baseline_atoms = atom_field_signature(&graph);

    let mut wrapper = InvalidatedGraph::new(graph);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    match ImpossibleHCountMutator.mutate_in_place(&mut wrapper, &mut rng) {
        Err(err) => assert_well_defined_error(err),
        Ok(()) => {
            assert_ne!(
                visible_invariant_signature(&graph),
                visible_invariant_signature(&wrapper),
                "ImpossibleHCountMutator: visible R0 invariant signature unchanged",
            );
            assert!(
                ImpossibleHCountPredicate.check(&wrapper),
                "ImpossibleHCountPredicate did not fire",
            );

            let mut wrapper2 = InvalidatedGraph::new(graph);
            let mut rng2 = ChaCha8Rng::seed_from_u64(seed);
            ImpossibleHCountMutator
                .mutate_in_place(&mut wrapper2, &mut rng2)
                .expect("deterministic re-run must succeed");
            assert_eq!(
                count_ecfp_fingerprint(&wrapper),
                count_ecfp_fingerprint(&wrapper2),
                "ImpossibleHCountMutator is non-deterministic for a fixed seed",
            );

            assert_mutator_minimal_diff_atoms(&graph, &wrapper, 1);
            assert_ne!(
                baseline_atoms,
                atom_field_signature(&wrapper),
                "ImpossibleHCountMutator: atom field multiset unchanged",
            );

            let h_violation = (0..wrapper.atom_count()).any(|atom_id| {
                let f = wrapper.ecfp_atom_invariant_fields(atom_id);
                f.total_hydrogens > max_natural_valence(f.atomic_number) + 1
            });
            assert!(
                h_violation,
                "ImpossibleHCountMutator did not produce total_hydrogens > max + 1",
            );
        }
    }
}

/// Enforces every invariant in the per-mutator contract for
/// `ImpossibleChargeMutator`.
pub fn fuzz_impossible_charge_invariants(smiles: Smiles, seed: u64) {
    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };
    let baseline_atoms = atom_field_signature(&graph);

    let mut wrapper = InvalidatedGraph::new(graph);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    match ImpossibleChargeMutator.mutate_in_place(&mut wrapper, &mut rng) {
        Err(err) => assert_well_defined_error(err),
        Ok(()) => {
            assert_ne!(
                visible_invariant_signature(&graph),
                visible_invariant_signature(&wrapper),
                "ImpossibleChargeMutator: visible R0 invariant signature unchanged",
            );
            assert!(
                ImpossibleChargePredicate.check(&wrapper),
                "ImpossibleChargePredicate did not fire",
            );

            let mut wrapper2 = InvalidatedGraph::new(graph);
            let mut rng2 = ChaCha8Rng::seed_from_u64(seed);
            ImpossibleChargeMutator
                .mutate_in_place(&mut wrapper2, &mut rng2)
                .expect("deterministic re-run must succeed");
            assert_eq!(
                count_ecfp_fingerprint(&wrapper),
                count_ecfp_fingerprint(&wrapper2),
                "ImpossibleChargeMutator is non-deterministic for a fixed seed",
            );

            assert_mutator_minimal_diff_atoms(&graph, &wrapper, 1);
            assert_ne!(
                baseline_atoms,
                atom_field_signature(&wrapper),
                "ImpossibleChargeMutator: atom field multiset unchanged",
            );

            let charge_violation = (0..wrapper.atom_count()).any(|atom_id| {
                wrapper
                    .ecfp_atom_invariant_fields(atom_id)
                    .formal_charge
                    .unsigned_abs()
                    > 6
            });
            assert!(
                charge_violation,
                "ImpossibleChargeMutator did not produce |formal_charge| > 6",
            );
        }
    }
}

/// Enforces every invariant in the per-mutator contract for
/// `ImpossibleIsotopeMutator`.
pub fn fuzz_impossible_isotope_invariants(smiles: Smiles, seed: u64) {
    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };
    let baseline_atoms = atom_field_signature(&graph);

    let mut wrapper = InvalidatedGraph::new(graph);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    match ImpossibleIsotopeMutator.mutate_in_place(&mut wrapper, &mut rng) {
        Err(err) => assert_well_defined_error(err),
        Ok(()) => {
            assert_ne!(
                visible_invariant_signature(&graph),
                visible_invariant_signature(&wrapper),
                "ImpossibleIsotopeMutator: visible R0 invariant signature unchanged",
            );
            assert!(
                ImpossibleIsotopePredicate.check(&wrapper),
                "ImpossibleIsotopePredicate did not fire",
            );

            let mut wrapper2 = InvalidatedGraph::new(graph);
            let mut rng2 = ChaCha8Rng::seed_from_u64(seed);
            ImpossibleIsotopeMutator
                .mutate_in_place(&mut wrapper2, &mut rng2)
                .expect("deterministic re-run must succeed");
            assert_eq!(
                count_ecfp_fingerprint(&wrapper),
                count_ecfp_fingerprint(&wrapper2),
                "ImpossibleIsotopeMutator is non-deterministic for a fixed seed",
            );

            assert_mutator_minimal_diff_atoms(&graph, &wrapper, 1);
            assert_ne!(
                baseline_atoms,
                atom_field_signature(&wrapper),
                "ImpossibleIsotopeMutator: atom field multiset unchanged",
            );

            let isotope_violation = (0..wrapper.atom_count()).any(|atom_id| {
                wrapper
                    .ecfp_atom_invariant_fields(atom_id)
                    .delta_mass
                    .unsigned_abs()
                    > 80
            });
            assert!(
                isotope_violation,
                "ImpossibleIsotopeMutator did not produce |delta_mass| > 80",
            );
        }
    }
}

/// Enforces every invariant in the per-mutator contract for
/// `ImpossibleRingFlagMutator`.
pub fn fuzz_impossible_ring_flag_invariants(smiles: Smiles, seed: u64) {
    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };
    let baseline_atoms = atom_field_signature(&graph);

    let mut wrapper = InvalidatedGraph::new(graph);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    match ImpossibleRingFlagMutator.mutate_in_place(&mut wrapper, &mut rng) {
        Err(err) => assert_well_defined_error(err),
        Ok(()) => {
            assert_ne!(
                visible_invariant_signature(&graph),
                visible_invariant_signature(&wrapper),
                "ImpossibleRingFlagMutator: visible R0 invariant signature unchanged",
            );
            assert!(
                ImpossibleRingFlagPredicate.check(&wrapper),
                "ImpossibleRingFlagPredicate did not fire",
            );

            let mut wrapper2 = InvalidatedGraph::new(graph);
            let mut rng2 = ChaCha8Rng::seed_from_u64(seed);
            ImpossibleRingFlagMutator
                .mutate_in_place(&mut wrapper2, &mut rng2)
                .expect("deterministic re-run must succeed");
            assert_eq!(
                count_ecfp_fingerprint(&wrapper),
                count_ecfp_fingerprint(&wrapper2),
                "ImpossibleRingFlagMutator is non-deterministic for a fixed seed",
            );

            assert_mutator_minimal_diff_atoms(&graph, &wrapper, 1);
            assert_ne!(
                baseline_atoms,
                atom_field_signature(&wrapper),
                "ImpossibleRingFlagMutator: atom field multiset unchanged",
            );

            // Predicate purity already covers the class-specific minimum-
            // violation invariant (atom flagged in-ring with < 2 heavy
            // neighbours), so we don't duplicate the check here.
        }
    }
}

/// Enforces every invariant in the per-mutator contract for
/// `ImpossibleBondTypeMutator`. This is the only bond-channel mutator: it
/// must leave atom invariants alone and change exactly one bond.
pub fn fuzz_impossible_bond_type_invariants(smiles: Smiles, seed: u64) {
    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };
    let baseline_bonds = bond_invariant_signature(&graph);

    let mut wrapper = InvalidatedGraph::new(graph);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    match ImpossibleBondTypeMutator.mutate_in_place(&mut wrapper, &mut rng) {
        Err(err) => assert_well_defined_error(err),
        Ok(()) => {
            // No `visible_invariant_signature` check here — the bond-channel
            // mutator must NOT change any atom's R0 invariant (asserted
            // separately below via `assert_mutator_no_atom_diff`). The
            // bond-side coverage comes from `assert_mutator_minimal_diff_bonds`
            // and the bond_invariant_signature comparison, which are both
            // unfolded and therefore collision-proof.
            assert!(
                ImpossibleBondTypePredicate.check(&wrapper),
                "ImpossibleBondTypePredicate did not fire",
            );

            let mut wrapper2 = InvalidatedGraph::new(graph);
            let mut rng2 = ChaCha8Rng::seed_from_u64(seed);
            ImpossibleBondTypeMutator
                .mutate_in_place(&mut wrapper2, &mut rng2)
                .expect("deterministic re-run must succeed");
            assert_eq!(
                count_ecfp_fingerprint(&wrapper),
                count_ecfp_fingerprint(&wrapper2),
                "ImpossibleBondTypeMutator is non-deterministic for a fixed seed",
            );

            // Bond-channel minimality: zero atom changes, exactly one bond.
            assert_mutator_no_atom_diff(&graph, &wrapper);
            assert_mutator_minimal_diff_bonds(&graph, &wrapper, 1);

            let mutated_bonds = bond_invariant_signature(&wrapper);
            assert_ne!(
                baseline_bonds, mutated_bonds,
                "ImpossibleBondTypeMutator: bond invariant signature unchanged",
            );
            let off_set_bond = mutated_bonds
                .iter()
                .any(|&(_, _, inv)| !matches!(inv, 1 | 2 | 3 | 12));
            assert!(
                off_set_bond,
                "ImpossibleBondTypeMutator did not produce a bond invariant outside {{1,2,3,12}}",
            );
        }
    }
}

/// Enforces every invariant in the per-mutator contract for
/// `TopologicalPathologyMutator`.
pub fn fuzz_topological_pathology_invariants(smiles: Smiles, seed: u64) {
    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };
    let baseline_atoms = atom_field_signature(&graph);

    let mut wrapper = InvalidatedGraph::new(graph);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    match TopologicalPathologyMutator.mutate_in_place(&mut wrapper, &mut rng) {
        Err(err) => assert_well_defined_error(err),
        Ok(()) => {
            assert_ne!(
                visible_invariant_signature(&graph),
                visible_invariant_signature(&wrapper),
                "TopologicalPathologyMutator: visible R0 invariant signature unchanged",
            );
            assert!(
                TopologicalPathologyPredicate.check(&wrapper),
                "TopologicalPathologyPredicate did not fire",
            );

            let mut wrapper2 = InvalidatedGraph::new(graph);
            let mut rng2 = ChaCha8Rng::seed_from_u64(seed);
            TopologicalPathologyMutator
                .mutate_in_place(&mut wrapper2, &mut rng2)
                .expect("deterministic re-run must succeed");
            assert_eq!(
                count_ecfp_fingerprint(&wrapper),
                count_ecfp_fingerprint(&wrapper2),
                "TopologicalPathologyMutator is non-deterministic for a fixed seed",
            );

            assert_mutator_minimal_diff_atoms(&graph, &wrapper, 1);
            assert_ne!(
                baseline_atoms,
                atom_field_signature(&wrapper),
                "TopologicalPathologyMutator: atom field multiset unchanged",
            );

            // Class-specific minimum-violation invariant: predicate purity
            // already covers it (some aromatic bond's endpoints disagree on
            // in_ring), so we rely on the predicate-fired assertion above.
        }
    }
}

/// Enforces every invariant in the composition contract for
/// `MutatorMix::sample`:
///
/// 1. Doesn't panic.
/// 2. On `Ok((wrapper, label))`:
///    a. `label.count() >= 1` (predicate-driven labelling).
///    b. `1 <= label.count() <= MutatorMix::k_max()` (composition bound,
///       loosened to 8 because secondary co-firings can briefly exceed
///       `k_max`).
///    c. `wrapper`'s count-ECFP differs from the baseline.
///    d. Determinism: re-running with the same seed produces an identical
///       fingerprint.
/// 3. On `Err`, the error is one of the four well-defined variants.
pub fn fuzz_mutator_mix_invariants(smiles: Smiles, seed: u64) {
    use finge_rs::{MutatorMix, ViolationClass, ViolationPredicate};
    let mut scratch = SmilesRdkitScratch::default();
    let Some(graph) = prepare_graph(&mut scratch, &smiles) else {
        return;
    };

    let mix = MutatorMix::<finge_rs::smiles_support::SmilesRdkitGraph<'_>>::with_default_mutators_and_predicates()
        .with_collision_filter(MUTATOR_FP_RADIUS, MUTATOR_FP_SIZE);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);
    match mix.sample(graph, &mut rng) {
        Err(err) => assert_well_defined_error(err),
        Ok((wrapper, label)) => {
            let count = label.count();
            assert!(
                count >= 1,
                "successful MutatorMix::sample must set at least one bit",
            );
            assert!(
                count <= ViolationClass::COUNT as u32,
                "label count {count} exceeds total class count",
            );

            // Collision-proof equivalent of "count-ECFP changed": every
            // successful `MutatorMix::sample` must leave the wrapper with at
            // least one override that actually changes a pre-hash atom
            // tuple, bond invariant, or ignored flag. Folded-ECFP-unchanged
            // assertions on the 65 536-bin fingerprint produced false
            // positives via hash collisions (see crash-9f067edf… for the
            // composition-cancellation case and the per-mutator crashes
            // such as `B[ClH]` / `[84Ga][61Cu]` for raw collisions).
            assert!(
                wrapper.has_effective_perturbation(),
                "MutatorMix::sample returned Ok with no effective perturbation",
            );
            // And: at least one change must actually reach the ECFP
            // emission path — either an R0 invariant on a non-ignored atom
            // or a bond invariant. `has_effective_perturbation` alone
            // would silently accept a hypothetical composition whose every
            // override landed on ECFP-ignored atoms.
            let invariant_changed =
                visible_invariant_signature(&graph) != visible_invariant_signature(&wrapper);
            let bonds_changed = bond_invariant_signature(&graph) != bond_invariant_signature(&wrapper);
            assert!(
                invariant_changed || bonds_changed,
                "MutatorMix::sample produced a wrapper whose perturbation is invisible to ECFP",
            );

            // Every set bit must be backed by a firing predicate.
            for class in ViolationClass::ALL {
                if label.is_set(class) {
                    let fires = match class {
                        ViolationClass::ImpossibleAtomicNumber => {
                            ImpossibleAtomicNumberPredicate.check(&wrapper)
                        }
                        ViolationClass::Hypervalent => HypervalentPredicate.check(&wrapper),
                        ViolationClass::ImpossibleHCount => {
                            ImpossibleHCountPredicate.check(&wrapper)
                        }
                        ViolationClass::ImpossibleCharge => {
                            ImpossibleChargePredicate.check(&wrapper)
                        }
                        ViolationClass::ImpossibleIsotope => {
                            ImpossibleIsotopePredicate.check(&wrapper)
                        }
                        ViolationClass::ImpossibleRingFlag => {
                            ImpossibleRingFlagPredicate.check(&wrapper)
                        }
                        ViolationClass::ImpossibleBondType => {
                            ImpossibleBondTypePredicate.check(&wrapper)
                        }
                        ViolationClass::TopologicalPathology => {
                            TopologicalPathologyPredicate.check(&wrapper)
                        }
                    };
                    assert!(
                        fires,
                        "label bit {class:?} set but predicate did not fire",
                    );
                }
            }

            // Determinism: re-running with the same seed produces an
            // identical fingerprint.
            let mut rng2 = ChaCha8Rng::seed_from_u64(seed);
            let (wrapper2, label2) = mix
                .sample(graph, &mut rng2)
                .expect("deterministic re-run must succeed");
            assert_eq!(
                count_ecfp_fingerprint(&wrapper),
                count_ecfp_fingerprint(&wrapper2),
                "MutatorMix::sample is non-deterministic for a fixed seed",
            );
            assert_eq!(
                label.as_u8(),
                label2.as_u8(),
                "MutatorMix::sample produced different labels for the same seed",
            );
        }
    }
}

pub fn fuzz_maccs(smiles: Smiles) {
    let fingerprint = MaccsFingerprint::new().expect("MACCS SMARTS should compile");
    let target = PreparedTarget::new(smiles.clone());
    let bits = fingerprint.compute(&target);
    assert_bit_fingerprint_basics(&bits);

    let explicit_h = smiles.with_explicit_hydrogens();
    let explicit_h_target = PreparedTarget::new(explicit_h);
    let explicit_h_bits = fingerprint.compute(&explicit_h_target);
    assert_bit_fingerprint_basics(&explicit_h_bits);

    let second_run = fingerprint.compute(&target);
    assert_eq!(
        bits.active_bits().collect::<Vec<_>>(),
        second_run.active_bits().collect::<Vec<_>>()
    );
}
