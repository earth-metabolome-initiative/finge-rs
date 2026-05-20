//! [`ViolationPredicate`] — pure checks for the eight ECFP-detectable
//! chemistry violations.
//!
//! Each predicate reads field tuples through
//! [`EcfpGraph::ecfp_atom_invariant_fields`] (and, for the two
//! bond-related classes, [`EcfpGraph::ecfp_bond_invariant`] /
//! [`crate::MolecularGraph::bonds`]) and applies a chemistry rule. They work on any
//! `EcfpGraph` — both valid molecules (where every predicate returns `false`
//! for the canonical positive corpus) and mutated wrappers from
//! [`crate::InvalidatedGraph`] (where the matching predicate fires).

use crate::{
    mutations::violation_class::ViolationClass,
    traits::{EcfpGraph, MolecularAtom, MolecularBond},
};

/// Static identity of the [`ViolationClass`] a predicate detects.
///
/// Split out from [`ViolationPredicate`] (which is generic over `G`) so that
/// callers — and the `MutatorMix` builder — can recover the class without
/// having to pick a concrete graph type for type inference.
pub trait PredicateClass {
    /// Returns the violation class this predicate detects.
    fn class(&self) -> ViolationClass;
}

/// A pure-function detector for one [`ViolationClass`].
pub trait ViolationPredicate<G>: PredicateClass
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    /// Returns whether `graph` exhibits this predicate's violation.
    fn check(&self, graph: &G) -> bool;
}

/// Returns the maximum plausible total degree (heavy + hydrogens) for an atom
/// of atomic number `z`.
///
/// Values are coarse but defensible: they cover common main-group hypervalent
/// states (S(VI), P(V), halogen oxoacids ClO4-, IO4-) so that real molecules
/// in the canonical positive corpus pass cleanly, and mutations that push
/// `total_degree` past these ceilings — typically `≥ 10` — are caught. Unknown
/// or super-heavy elements (which would normally trigger
/// [`ViolationClass::ImpossibleAtomicNumber`] first) fall through to a
/// permissive `8`.
#[inline]
#[must_use]
pub const fn max_natural_valence(z: u32) -> u32 {
    match z {
        1 => 1,                          // H
        2 | 10 | 18 | 36 | 54 | 86 => 0, // noble gases
        6 => 4,                          // C
        7 => 5,                          // N (ammonium, NO2/NO3)
        8 => 3,                          // O (oxonium)
        9 | 17 | 35 | 53 | 85 => 7,      // halogens (oxoacid hypervalence)
        15 => 5,                         // P
        16 => 6,                         // S (SF6)
        _ => 8,                          // permissive fallback
    }
}

/// Returns whether neighbour `neighbor_id` of `atom_id` is a heavy (non-H) atom.
#[inline]
fn is_heavy_neighbor<G>(graph: &G, neighbor_id: usize) -> bool
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    graph.ecfp_atom_invariant_fields(neighbor_id).atomic_number != 1
}

/// Counts how many of `atom_id`'s incident bonds connect it to a heavy atom.
#[inline]
fn heavy_bond_count<G>(graph: &G, atom_id: usize) -> usize
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    graph
        .bonds(atom_id)
        .filter(|bond| {
            let other = if bond.source() == atom_id {
                bond.target()
            } else {
                bond.source()
            };
            is_heavy_neighbor(graph, other)
        })
        .count()
}

/// Iterates each undirected bond once, yielding `(source_id, target_id, bond)`
/// with `source_id < target_id`.
fn iter_unique_bonds<G>(graph: &G) -> impl Iterator<Item = (usize, usize, G::Bond)> + '_
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    (0..graph.atom_count()).flat_map(move |source| {
        graph.bonds(source).filter_map(move |bond| {
            let other = if bond.source() == source {
                bond.target()
            } else if bond.target() == source {
                bond.source()
            } else {
                return None;
            };
            if other > source {
                Some((source, other, bond))
            } else {
                None
            }
        })
    })
}

// ---------------------------------------------------------------------------
// Class 0 — ImpossibleAtomicNumber
// ---------------------------------------------------------------------------

/// Returns whether `graph` contains an atom whose `Z` is outside `1..=118`.
#[inline]
#[must_use]
pub fn is_impossible_atomic_number<G>(graph: &G) -> bool
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    (0..graph.atom_count()).any(|atom_id| {
        let z = graph.ecfp_atom_invariant_fields(atom_id).atomic_number;
        z == 0 || z > 118
    })
}

/// Predicate matching [`ViolationClass::ImpossibleAtomicNumber`].
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleAtomicNumberPredicate;

impl PredicateClass for ImpossibleAtomicNumberPredicate {
    #[inline]
    fn class(&self) -> ViolationClass {
        ViolationClass::ImpossibleAtomicNumber
    }
}

impl<G> ViolationPredicate<G> for ImpossibleAtomicNumberPredicate
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    #[inline]
    fn check(&self, graph: &G) -> bool {
        is_impossible_atomic_number(graph)
    }
}

// ---------------------------------------------------------------------------
// Class 1 — Hypervalent
// ---------------------------------------------------------------------------

/// Returns whether `graph` contains an atom whose `total_degree` exceeds
/// [`max_natural_valence`] for its atomic number.
#[inline]
#[must_use]
pub fn is_hypervalent<G>(graph: &G) -> bool
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    (0..graph.atom_count()).any(|atom_id| {
        let fields = graph.ecfp_atom_invariant_fields(atom_id);
        fields.total_degree > max_natural_valence(fields.atomic_number)
    })
}

/// Predicate matching [`ViolationClass::Hypervalent`].
#[derive(Clone, Copy, Debug, Default)]
pub struct HypervalentPredicate;

impl PredicateClass for HypervalentPredicate {
    #[inline]
    fn class(&self) -> ViolationClass {
        ViolationClass::Hypervalent
    }
}

impl<G> ViolationPredicate<G> for HypervalentPredicate
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    #[inline]
    fn check(&self, graph: &G) -> bool {
        is_hypervalent(graph)
    }
}

// ---------------------------------------------------------------------------
// Class 2 — ImpossibleHCount
// ---------------------------------------------------------------------------

/// Returns whether `graph` contains an atom whose `total_hydrogens` exceeds
/// `max_natural_valence(Z) + 1`.
#[inline]
#[must_use]
pub fn has_impossible_hydrogen_count<G>(graph: &G) -> bool
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    (0..graph.atom_count()).any(|atom_id| {
        let fields = graph.ecfp_atom_invariant_fields(atom_id);
        fields.total_hydrogens > max_natural_valence(fields.atomic_number) + 1
    })
}

/// Predicate matching [`ViolationClass::ImpossibleHCount`].
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleHCountPredicate;

impl PredicateClass for ImpossibleHCountPredicate {
    #[inline]
    fn class(&self) -> ViolationClass {
        ViolationClass::ImpossibleHCount
    }
}

impl<G> ViolationPredicate<G> for ImpossibleHCountPredicate
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    #[inline]
    fn check(&self, graph: &G) -> bool {
        has_impossible_hydrogen_count(graph)
    }
}

// ---------------------------------------------------------------------------
// Class 3 — ImpossibleCharge
// ---------------------------------------------------------------------------

/// Returns whether `graph` contains an atom whose `|formal_charge| > 6`.
#[inline]
#[must_use]
pub fn has_impossible_charge<G>(graph: &G) -> bool
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    (0..graph.atom_count()).any(|atom_id| {
        graph
            .ecfp_atom_invariant_fields(atom_id)
            .formal_charge
            .unsigned_abs()
            > 6
    })
}

/// Predicate matching [`ViolationClass::ImpossibleCharge`].
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleChargePredicate;

impl PredicateClass for ImpossibleChargePredicate {
    #[inline]
    fn class(&self) -> ViolationClass {
        ViolationClass::ImpossibleCharge
    }
}

impl<G> ViolationPredicate<G> for ImpossibleChargePredicate
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    #[inline]
    fn check(&self, graph: &G) -> bool {
        has_impossible_charge(graph)
    }
}

// ---------------------------------------------------------------------------
// Class 4 — ImpossibleIsotope
// ---------------------------------------------------------------------------

/// Returns whether `graph` contains an atom whose `|delta_mass| > 80`.
///
/// Coarse but defensible: no real isotope's mass differs from the standard
/// atomic mass by more than ~50, and never by 80+. The matching mutator emits
/// `|Δmass| ≥ 100`.
#[inline]
#[must_use]
pub fn has_impossible_isotope<G>(graph: &G) -> bool
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    (0..graph.atom_count()).any(|atom_id| {
        graph
            .ecfp_atom_invariant_fields(atom_id)
            .delta_mass
            .unsigned_abs()
            > 80
    })
}

/// Predicate matching [`ViolationClass::ImpossibleIsotope`].
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleIsotopePredicate;

impl PredicateClass for ImpossibleIsotopePredicate {
    #[inline]
    fn class(&self) -> ViolationClass {
        ViolationClass::ImpossibleIsotope
    }
}

impl<G> ViolationPredicate<G> for ImpossibleIsotopePredicate
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    #[inline]
    fn check(&self, graph: &G) -> bool {
        has_impossible_isotope(graph)
    }
}

// ---------------------------------------------------------------------------
// Class 5 — ImpossibleRingFlag
// ---------------------------------------------------------------------------

/// Returns whether `graph` contains an atom flagged as in-ring but with fewer
/// than 2 heavy neighbours (i.e. can't be a cycle member).
#[inline]
#[must_use]
pub fn has_impossible_ring_flag<G>(graph: &G) -> bool
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    (0..graph.atom_count()).any(|atom_id| {
        let fields = graph.ecfp_atom_invariant_fields(atom_id);
        fields.in_ring && heavy_bond_count(graph, atom_id) < 2
    })
}

/// Predicate matching [`ViolationClass::ImpossibleRingFlag`].
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleRingFlagPredicate;

impl PredicateClass for ImpossibleRingFlagPredicate {
    #[inline]
    fn class(&self) -> ViolationClass {
        ViolationClass::ImpossibleRingFlag
    }
}

impl<G> ViolationPredicate<G> for ImpossibleRingFlagPredicate
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    #[inline]
    fn check(&self, graph: &G) -> bool {
        has_impossible_ring_flag(graph)
    }
}

// ---------------------------------------------------------------------------
// Class 6 — ImpossibleBondType
// ---------------------------------------------------------------------------

/// Returns whether `graph` contains any bond whose invariant is outside the
/// RDKit-recognised set `{1, 2, 3, 12}`.
#[inline]
#[must_use]
pub fn has_impossible_bond_type<G>(graph: &G) -> bool
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    iter_unique_bonds(graph).any(|(_, _, bond)| {
        let inv = graph.ecfp_bond_invariant(&bond, true);
        !matches!(inv, 1 | 2 | 3 | 12)
    })
}

/// Predicate matching [`ViolationClass::ImpossibleBondType`].
#[derive(Clone, Copy, Debug, Default)]
pub struct ImpossibleBondTypePredicate;

impl PredicateClass for ImpossibleBondTypePredicate {
    #[inline]
    fn class(&self) -> ViolationClass {
        ViolationClass::ImpossibleBondType
    }
}

impl<G> ViolationPredicate<G> for ImpossibleBondTypePredicate
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    #[inline]
    fn check(&self, graph: &G) -> bool {
        has_impossible_bond_type(graph)
    }
}

// ---------------------------------------------------------------------------
// Class 7 — TopologicalPathology
// ---------------------------------------------------------------------------

/// Returns whether `graph` has an aromatic bond (invariant `12`) whose
/// endpoints disagree on ring membership.
#[inline]
#[must_use]
pub fn has_topological_pathology<G>(graph: &G) -> bool
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    iter_unique_bonds(graph).any(|(source, target, bond)| {
        if graph.ecfp_bond_invariant(&bond, true) != 12 {
            return false;
        }
        let src_in_ring = graph.ecfp_atom_invariant_fields(source).in_ring;
        let tgt_in_ring = graph.ecfp_atom_invariant_fields(target).in_ring;
        src_in_ring != tgt_in_ring
    })
}

/// Predicate matching [`ViolationClass::TopologicalPathology`].
#[derive(Clone, Copy, Debug, Default)]
pub struct TopologicalPathologyPredicate;

impl PredicateClass for TopologicalPathologyPredicate {
    #[inline]
    fn class(&self) -> ViolationClass {
        ViolationClass::TopologicalPathology
    }
}

impl<G> ViolationPredicate<G> for TopologicalPathologyPredicate
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    #[inline]
    fn check(&self, graph: &G) -> bool {
        has_topological_pathology(graph)
    }
}

#[cfg(test)]
mod tests {
    use super::{
        HypervalentPredicate, ImpossibleAtomicNumberPredicate, ImpossibleBondTypePredicate,
        ImpossibleChargePredicate, ImpossibleHCountPredicate, ImpossibleIsotopePredicate,
        ImpossibleRingFlagPredicate, PredicateClass, TopologicalPathologyPredicate,
        ViolationPredicate, has_impossible_bond_type, has_impossible_charge,
        has_impossible_hydrogen_count, has_impossible_isotope, has_impossible_ring_flag,
        has_topological_pathology, is_hypervalent, is_impossible_atomic_number,
        max_natural_valence,
    };
    use crate::{
        InvalidatedGraph, ViolationClass, smiles_support_impl::SmilesRdkitScratch,
        traits::MolecularGraph as _,
    };

    /// Canonical positive corpus — every predicate must return `false` for
    /// every molecule here. Hits both aliphatic and aromatic, neutral and
    /// heteroatom cases.
    const POSITIVE_CORPUS: &[&str] = &[
        "C",
        "CCO",
        "CC(=O)O",
        "c1ccccc1",
        "C1CCCCC1",
        "n1ccccc1",
        "CC(=O)c1ccccc1",
    ];

    fn parse(smiles: &str) -> (SmilesRdkitScratch, smiles_parser::smiles::Smiles) {
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("fixture SMILES should parse");
        (SmilesRdkitScratch::default(), parsed)
    }

    #[test]
    fn max_natural_valence_table_matches_chemistry() {
        assert_eq!(max_natural_valence(1), 1); // H
        assert_eq!(max_natural_valence(6), 4); // C
        assert_eq!(max_natural_valence(7), 5); // N
        assert_eq!(max_natural_valence(8), 3); // O
        assert_eq!(max_natural_valence(15), 5); // P
        assert_eq!(max_natural_valence(16), 6); // S
        assert_eq!(max_natural_valence(17), 7); // Cl
        assert_eq!(max_natural_valence(2), 0); // He
        assert_eq!(max_natural_valence(200), 8); // permissive fallback
    }

    #[test]
    fn every_predicate_is_silent_on_positive_corpus() {
        for smiles in POSITIVE_CORPUS {
            let (mut scratch, parsed) = parse(smiles);
            let graph = scratch.prepare(&parsed);
            assert!(
                !is_impossible_atomic_number(&graph),
                "ImpossibleAtomicNumber fired on {smiles}"
            );
            assert!(!is_hypervalent(&graph), "Hypervalent fired on {smiles}");
            assert!(
                !has_impossible_hydrogen_count(&graph),
                "ImpossibleHCount fired on {smiles}"
            );
            assert!(
                !has_impossible_charge(&graph),
                "ImpossibleCharge fired on {smiles}"
            );
            assert!(
                !has_impossible_isotope(&graph),
                "ImpossibleIsotope fired on {smiles}"
            );
            assert!(
                !has_impossible_ring_flag(&graph),
                "ImpossibleRingFlag fired on {smiles}"
            );
            assert!(
                !has_impossible_bond_type(&graph),
                "ImpossibleBondType fired on {smiles}"
            );
            assert!(
                !has_topological_pathology(&graph),
                "TopologicalPathology fired on {smiles}"
            );
        }
    }

    #[test]
    fn impossible_atomic_number_predicate_class_and_fires() {
        let (mut scratch, parsed) = parse("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_atomic_number_override(0, 200);
        let p = ImpossibleAtomicNumberPredicate;
        assert_eq!(p.class(), ViolationClass::ImpossibleAtomicNumber);
        assert!(p.check(&wrapper));
    }

    #[test]
    fn impossible_atomic_number_also_fires_on_zero() {
        let (mut scratch, parsed) = parse("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_atomic_number_override(1, 0);
        assert!(is_impossible_atomic_number(&wrapper));
    }

    #[test]
    fn hypervalent_predicate_class_and_fires() {
        let (mut scratch, parsed) = parse("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        // Z=6 (C), max_natural_valence=4. Override total_degree to 12.
        wrapper.set_total_degree_override(0, 12);
        let p = HypervalentPredicate;
        assert_eq!(p.class(), ViolationClass::Hypervalent);
        assert!(p.check(&wrapper));
    }

    #[test]
    fn impossible_h_count_predicate_class_and_fires() {
        let (mut scratch, parsed) = parse("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        // C max+1 = 5. Override to 30.
        wrapper.set_total_hydrogens_override(0, 30);
        let p = ImpossibleHCountPredicate;
        assert_eq!(p.class(), ViolationClass::ImpossibleHCount);
        assert!(p.check(&wrapper));
    }

    #[test]
    fn impossible_charge_predicate_class_and_fires() {
        let (mut scratch, parsed) = parse("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_formal_charge_override(0, 9);
        let p = ImpossibleChargePredicate;
        assert_eq!(p.class(), ViolationClass::ImpossibleCharge);
        assert!(p.check(&wrapper));

        wrapper.set_formal_charge_override(0, -9);
        assert!(p.check(&wrapper));
    }

    #[test]
    fn impossible_isotope_predicate_class_and_fires() {
        let (mut scratch, parsed) = parse("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_delta_mass_override(0, 100);
        let p = ImpossibleIsotopePredicate;
        assert_eq!(p.class(), ViolationClass::ImpossibleIsotope);
        assert!(p.check(&wrapper));
    }

    #[test]
    fn impossible_ring_flag_predicate_class_and_fires() {
        let (mut scratch, parsed) = parse("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        // Atom 0 in "CCO" is the terminal CH3 (1 heavy neighbour); marking it
        // in-ring is the canonical Class-5 violation.
        wrapper.set_in_ring_override(0, true);
        let p = ImpossibleRingFlagPredicate;
        assert_eq!(p.class(), ViolationClass::ImpossibleRingFlag);
        assert!(p.check(&wrapper));
    }

    #[test]
    fn impossible_bond_type_predicate_class_and_fires() {
        let (mut scratch, parsed) = parse("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_bond_invariant_override(0, 1, 0xDEAD);
        let p = ImpossibleBondTypePredicate;
        assert_eq!(p.class(), ViolationClass::ImpossibleBondType);
        assert!(p.check(&wrapper));
    }

    #[test]
    fn topological_pathology_predicate_class_and_fires() {
        // Benzene: every bond aromatic, every atom in-ring. Flip atom 0's
        // ring flag and the predicate should report inconsistency across the
        // aromatic bonds that touch it.
        let (mut scratch, parsed) = parse("c1ccccc1");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_in_ring_override(0, false);
        let p = TopologicalPathologyPredicate;
        assert_eq!(p.class(), ViolationClass::TopologicalPathology);
        assert!(p.check(&wrapper));
    }

    #[test]
    fn topological_pathology_silent_when_no_aromatic_bond() {
        // Cyclohexane: ring atoms, but all single bonds. Flipping one in_ring
        // flag fires `ImpossibleRingFlag` (different class) but NOT
        // `TopologicalPathology`, which is gated on aromatic bonds.
        let (mut scratch, parsed) = parse("C1CCCCC1");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_in_ring_override(0, false);
        assert!(!has_topological_pathology(&wrapper));
    }

    #[test]
    fn impossible_ring_flag_does_not_fire_on_real_ring_atom_marked_non_ring() {
        // The mutator that flips a ring atom's in_ring flag to false belongs
        // to TopologicalPathology, not ImpossibleRingFlag — the latter is
        // specifically "in_ring set on an atom that can't be in a ring".
        let (mut scratch, parsed) = parse("c1ccccc1");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_in_ring_override(0, false);
        assert!(!has_impossible_ring_flag(&wrapper));
    }

    #[test]
    fn impossible_ring_flag_handles_atom_count_correctly() {
        // Verify the heavy_bond_count helper isn't tripped by hydrogens.
        let (mut scratch, parsed) = parse("CCO");
        let inner = scratch.prepare(&parsed);
        // Without overrides every predicate is silent.
        assert!(!has_impossible_ring_flag(&inner));
        assert_eq!(inner.atom_count(), 3);
    }
}
