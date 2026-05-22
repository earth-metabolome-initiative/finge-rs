//! [`InvalidatedGraph`] — an owning wrapper that delegates topology to an
//! inner [`EcfpGraph`] but selectively overrides the *chemistry fields* the
//! RDKit-style ECFP invariant is computed from.
//!
//! Mutators write per-channel field overrides (atomic number, total degree,
//! total hydrogens, formal charge, Δmass, in-ring) into this struct;
//! everything else is delegated straight to the inner graph. When the ECFP
//! algorithm asks for an atom invariant the wrapper reads the inner graph's
//! [`AtomInvariantFields`], applies any active overrides, and hashes the
//! perturbed tuple through the same [`AtomInvariantFields::to_invariant`]
//! path that the standard pipeline uses. The mutated atom therefore lives in
//! the same invariant space as a real (but chemistry-impossible) atom —
//! e.g. `Z = 200` carbons — rather than carrying a synthetic marker.
//!
//! Bond invariants remain raw `u32`s: the ECFP bond invariant is just the
//! RDKit bond-order code (`1`, `2`, `3`, `12`), so the "impossible bond type"
//! mutator writes a code outside that set directly.

use alloc::{collections::BTreeMap, vec, vec::Vec};

use geometric_traits::traits::{Graph, MonopartiteGraph, MonoplexGraph};

use crate::{
    atom_invariant_fields::AtomInvariantFields,
    traits::{
        AtomPairGraph, EcfpGraph, MolecularAtom, MolecularBond, MolecularGraph,
        RdkFingerprintGraph, TopologicalTorsionGraph,
    },
};

/// Per-channel chemistry field overrides written by built-in mutators.
///
/// Each `Some(value)` replaces the corresponding entry of the inner graph's
/// [`AtomInvariantFields`]; `None` lets the inner value flow through.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct AtomFieldOverride {
    /// Override the atom's `Z`.
    pub atomic_number: Option<u32>,
    /// Override the atom's `heavy_degree + total_hydrogens` count.
    pub total_degree: Option<u32>,
    /// Override the atom's `explicit + implicit` hydrogen count.
    pub total_hydrogens: Option<u32>,
    /// Override the atom's formal charge.
    pub formal_charge: Option<i32>,
    /// Override the atom's `isotope_mass - standard_atomic_mass` delta.
    pub delta_mass: Option<i32>,
    /// Override the atom's ring-membership flag.
    pub in_ring: Option<bool>,
}

impl AtomFieldOverride {
    /// Returns whether at least one channel of this override is active.
    #[inline]
    #[must_use]
    pub const fn is_active(self) -> bool {
        self.atomic_number.is_some()
            || self.total_degree.is_some()
            || self.total_hydrogens.is_some()
            || self.formal_charge.is_some()
            || self.delta_mass.is_some()
            || self.in_ring.is_some()
    }

    /// Applies this override on top of `fields` and returns the result.
    #[inline]
    #[must_use]
    pub fn apply(self, mut fields: AtomInvariantFields) -> AtomInvariantFields {
        if let Some(z) = self.atomic_number {
            fields.atomic_number = z;
        }
        if let Some(d) = self.total_degree {
            fields.total_degree = d;
        }
        if let Some(h) = self.total_hydrogens {
            fields.total_hydrogens = h;
        }
        if let Some(q) = self.formal_charge {
            fields.formal_charge = q;
        }
        if let Some(m) = self.delta_mass {
            fields.delta_mass = m;
        }
        if let Some(r) = self.in_ring {
            fields.in_ring = r;
        }
        fields
    }
}

/// Owning wrapper that delegates every graph operation to `inner` except for
/// ECFP-relevant atom invariants and bond invariants, which consult per-atom
/// field-override tables and a bond-invariant override map first.
#[derive(Debug, Clone)]
pub struct InvalidatedGraph<G> {
    inner: G,
    /// Per-atom field overrides, one slot per atom. `None` means "use the
    /// inner graph's fields verbatim".
    atom_field_overrides: Vec<Option<AtomFieldOverride>>,
    /// Per-(min, max) bond endpoint override. Canonicalised so mutators can
    /// look up either endpoint order.
    bond_invariant_overrides: BTreeMap<(usize, usize), u32>,
    /// Per-atom "is ignored as an ECFP center" override.
    atom_is_ignored_overrides: Vec<Option<bool>>,
}

impl<G> InvalidatedGraph<G>
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    /// Wraps `graph` with empty override tables sized to the graph's atom
    /// count. The wrapper initially produces fingerprints identical to those
    /// of `graph`.
    #[inline]
    #[must_use]
    pub fn new(graph: G) -> Self {
        let atom_count = graph.atom_count();
        Self {
            inner: graph,
            atom_field_overrides: vec![None; atom_count],
            bond_invariant_overrides: BTreeMap::new(),
            atom_is_ignored_overrides: vec![None; atom_count],
        }
    }

    /// Returns a reference to the wrapped graph.
    #[inline]
    #[must_use]
    pub fn inner(&self) -> &G {
        &self.inner
    }

    /// Consumes the wrapper and returns the inner graph, discarding overrides.
    #[inline]
    #[must_use]
    pub fn into_inner(self) -> G {
        self.inner
    }

    /// Returns the active override for `atom_id`, if any.
    #[inline]
    #[must_use]
    pub fn atom_field_override(&self, atom_id: usize) -> Option<AtomFieldOverride> {
        self.atom_field_overrides.get(atom_id).copied().flatten()
    }

    /// Returns whether any override (atom-field, bond-invariant, or
    /// atom-is-ignored) is currently active.
    #[inline]
    #[must_use]
    pub fn has_any_override(&self) -> bool {
        !self.bond_invariant_overrides.is_empty()
            || self
                .atom_field_overrides
                .iter()
                .any(|slot| slot.is_some_and(AtomFieldOverride::is_active))
            || self.atom_is_ignored_overrides.iter().any(Option::is_some)
    }

    /// Clears all field overrides on `atom_id`.
    ///
    /// No-op when `atom_id` is out of range.
    #[inline]
    pub fn clear_atom_field_override(&mut self, atom_id: usize) {
        if let Some(slot) = self.atom_field_overrides.get_mut(atom_id) {
            *slot = None;
        }
    }

    /// Sets the atomic number returned for `atom_id`.
    #[inline]
    pub fn set_atomic_number_override(&mut self, atom_id: usize, atomic_number: u32) {
        self.modify_atom_field_override(atom_id, |ov| ov.atomic_number = Some(atomic_number));
    }

    /// Sets the total degree returned for `atom_id`.
    #[inline]
    pub fn set_total_degree_override(&mut self, atom_id: usize, total_degree: u32) {
        self.modify_atom_field_override(atom_id, |ov| ov.total_degree = Some(total_degree));
    }

    /// Sets the total hydrogen count returned for `atom_id`.
    #[inline]
    pub fn set_total_hydrogens_override(&mut self, atom_id: usize, total_hydrogens: u32) {
        self.modify_atom_field_override(atom_id, |ov| ov.total_hydrogens = Some(total_hydrogens));
    }

    /// Sets the formal charge returned for `atom_id`.
    #[inline]
    pub fn set_formal_charge_override(&mut self, atom_id: usize, formal_charge: i32) {
        self.modify_atom_field_override(atom_id, |ov| ov.formal_charge = Some(formal_charge));
    }

    /// Sets the Δmass returned for `atom_id`.
    #[inline]
    pub fn set_delta_mass_override(&mut self, atom_id: usize, delta_mass: i32) {
        self.modify_atom_field_override(atom_id, |ov| ov.delta_mass = Some(delta_mass));
    }

    /// Sets the ring-membership flag returned for `atom_id`.
    #[inline]
    pub fn set_in_ring_override(&mut self, atom_id: usize, in_ring: bool) {
        self.modify_atom_field_override(atom_id, |ov| ov.in_ring = Some(in_ring));
    }

    /// Sets the bond invariant returned for the bond between `a` and `b`.
    ///
    /// The endpoint order does not matter — `(a, b)` and `(b, a)` map to the
    /// same override.
    #[inline]
    pub fn set_bond_invariant_override(&mut self, a: usize, b: usize, value: u32) {
        let key = canonical_bond_key(a, b);
        self.bond_invariant_overrides.insert(key, value);
    }

    /// Clears any bond invariant override between `a` and `b`.
    #[inline]
    pub fn clear_bond_invariant_override(&mut self, a: usize, b: usize) {
        let key = canonical_bond_key(a, b);
        self.bond_invariant_overrides.remove(&key);
    }

    /// Sets the `ecfp_atom_is_ignored` override for `atom_id`.
    ///
    /// No-op when `atom_id` is out of range.
    #[inline]
    pub fn set_atom_is_ignored_override(&mut self, atom_id: usize, value: bool) {
        if let Some(slot) = self.atom_is_ignored_overrides.get_mut(atom_id) {
            *slot = Some(value);
        }
    }

    #[inline]
    fn modify_atom_field_override<F>(&mut self, atom_id: usize, mutate: F)
    where
        F: FnOnce(&mut AtomFieldOverride),
    {
        if let Some(slot) = self.atom_field_overrides.get_mut(atom_id) {
            let mut current = slot.unwrap_or_default();
            mutate(&mut current);
            *slot = Some(current);
        }
    }
}

#[inline]
fn canonical_bond_key(a: usize, b: usize) -> (usize, usize) {
    if a <= b { (a, b) } else { (b, a) }
}

impl<G> Graph for InvalidatedGraph<G>
where
    G: Graph,
{
    #[inline]
    fn has_nodes(&self) -> bool {
        self.inner.has_nodes()
    }

    #[inline]
    fn has_edges(&self) -> bool {
        self.inner.has_edges()
    }
}

impl<G> MonoplexGraph for InvalidatedGraph<G>
where
    G: MonoplexGraph,
{
    type Edge = G::Edge;
    type Edges = G::Edges;

    #[inline]
    fn edges(&self) -> &Self::Edges {
        self.inner.edges()
    }
}

impl<G> MonopartiteGraph for InvalidatedGraph<G>
where
    G: MonopartiteGraph,
{
    type NodeId = G::NodeId;
    type NodeSymbol = G::NodeSymbol;
    type Nodes = G::Nodes;

    #[inline]
    fn nodes_vocabulary(&self) -> &Self::Nodes {
        self.inner.nodes_vocabulary()
    }
}

impl<G> MolecularGraph for InvalidatedGraph<G>
where
    G: MolecularGraph,
    G::NodeSymbol: MolecularAtom,
{
    type Bond = G::Bond;

    #[inline]
    fn atom(&self, node_id: Self::NodeId) -> Option<&Self::NodeSymbol> {
        self.inner.atom(node_id)
    }

    #[inline]
    fn bonds(&self, node_id: Self::NodeId) -> impl Iterator<Item = Self::Bond> + '_ {
        self.inner.bonds(node_id)
    }
}

impl<G> EcfpGraph for InvalidatedGraph<G>
where
    G: EcfpGraph,
    G::NodeSymbol: MolecularAtom,
{
    #[inline]
    fn ecfp_atom_is_ignored(&self, atom_id: usize) -> bool {
        match self
            .atom_is_ignored_overrides
            .get(atom_id)
            .copied()
            .flatten()
        {
            Some(value) => value,
            None => self.inner.ecfp_atom_is_ignored(atom_id),
        }
    }

    #[inline]
    fn ecfp_atom_invariant_fields(&self, atom_id: usize) -> AtomInvariantFields {
        let fields = self.inner.ecfp_atom_invariant_fields(atom_id);
        match self.atom_field_overrides.get(atom_id).copied().flatten() {
            Some(ov) => ov.apply(fields),
            None => fields,
        }
    }

    // `ecfp_atom_invariant` deliberately uses the trait's default impl,
    // which hashes `ecfp_atom_invariant_fields(...).to_invariant(...)`.

    #[inline]
    fn ecfp_bond_invariant(&self, bond: &Self::Bond, use_bond_types: bool) -> u32 {
        let key = canonical_bond_key(bond.source(), bond.target());
        match self.bond_invariant_overrides.get(&key).copied() {
            Some(value) => value,
            None => self.inner.ecfp_bond_invariant(bond, use_bond_types),
        }
    }
}

impl<G> AtomPairGraph for InvalidatedGraph<G>
where
    G: AtomPairGraph,
    G::NodeSymbol: MolecularAtom,
{
    #[inline]
    fn atom_pair_atom_code(&self, atom_id: usize) -> u32 {
        self.inner.atom_pair_atom_code(atom_id)
    }
}

impl<G> TopologicalTorsionGraph for InvalidatedGraph<G>
where
    G: TopologicalTorsionGraph,
    G::NodeSymbol: MolecularAtom,
{
    #[inline]
    fn topological_torsion_atom_is_hydrogen(&self, atom_id: usize) -> bool {
        self.inner.topological_torsion_atom_is_hydrogen(atom_id)
    }

    #[inline]
    fn topological_torsion_atom_code(&self, atom_id: usize, branch_subtract: u8) -> u32 {
        self.inner
            .topological_torsion_atom_code(atom_id, branch_subtract)
    }
}

impl<G> RdkFingerprintGraph for InvalidatedGraph<G>
where
    G: RdkFingerprintGraph,
    G::NodeSymbol: MolecularAtom,
{
    #[inline]
    fn rdk_atom_invariant(&self, atom_id: usize) -> u32 {
        self.inner.rdk_atom_invariant(atom_id)
    }

    #[inline]
    fn rdk_atom_is_hydrogen(&self, atom_id: usize) -> bool {
        self.inner.rdk_atom_is_hydrogen(atom_id)
    }

    #[inline]
    fn rdk_atom_is_collapsible_hydrogen(&self, atom_id: usize) -> bool {
        self.inner.rdk_atom_is_collapsible_hydrogen(atom_id)
    }

    #[inline]
    fn rdk_bond_code(&self, bond: &Self::Bond, use_bond_order: bool) -> u32 {
        self.inner.rdk_bond_code(bond, use_bond_order)
    }
}

#[cfg(test)]
mod tests {
    use alloc::vec::Vec;

    use super::{AtomFieldOverride, InvalidatedGraph};
    use crate::{
        AtomInvariantFields, CountEcfpFingerprint, EcfpFingerprint, EcfpGraph, Fingerprint,
        MolecularGraph, smiles_support_impl::SmilesRdkitScratch,
    };

    fn prepared(smiles: &str) -> (SmilesRdkitScratch, smiles_parser::smiles::Smiles) {
        let parsed: smiles_parser::smiles::Smiles =
            smiles.parse().expect("test SMILES should parse");
        (SmilesRdkitScratch::default(), parsed)
    }

    #[test]
    fn identity_wrapper_produces_same_ecfp_as_inner() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);

        let direct = EcfpFingerprint::new(2, 2048).compute(&inner);
        let wrapped = EcfpFingerprint::new(2, 2048).compute(&InvalidatedGraph::new(inner));

        assert_eq!(
            direct.active_bits().collect::<Vec<_>>(),
            wrapped.active_bits().collect::<Vec<_>>(),
        );
    }

    #[test]
    fn identity_wrapper_produces_same_count_ecfp_as_inner() {
        let (mut scratch, parsed) = prepared("c1ccccc1O");
        let inner = scratch.prepare(&parsed);
        let baseline = CountEcfpFingerprint::new(2, 2048).compute(&inner);
        let wrapped = CountEcfpFingerprint::new(2, 2048).compute(&InvalidatedGraph::new(inner));
        assert_eq!(baseline, wrapped);
    }

    #[test]
    fn atomic_number_override_changes_fingerprint() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline = EcfpFingerprint::new(2, 2048).compute(&inner);

        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_atomic_number_override(0, 200);
        let mutated = EcfpFingerprint::new(2, 2048).compute(&wrapper);

        assert_ne!(
            baseline.active_bits().collect::<Vec<_>>(),
            mutated.active_bits().collect::<Vec<_>>(),
        );
    }

    #[test]
    fn total_degree_override_changes_fingerprint() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline = EcfpFingerprint::new(2, 2048).compute(&inner);

        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_total_degree_override(1, 12);
        let mutated = EcfpFingerprint::new(2, 2048).compute(&wrapper);

        assert_ne!(
            baseline.active_bits().collect::<Vec<_>>(),
            mutated.active_bits().collect::<Vec<_>>(),
        );
    }

    #[test]
    fn formal_charge_override_changes_invariant_fields() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline_fields = inner.ecfp_atom_invariant_fields(0);

        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_formal_charge_override(0, 7);
        let mutated_fields = wrapper.ecfp_atom_invariant_fields(0);

        assert_ne!(baseline_fields, mutated_fields);
        assert_eq!(mutated_fields.formal_charge, 7);
        // Other channels untouched.
        assert_eq!(mutated_fields.atomic_number, baseline_fields.atomic_number);
        assert_eq!(mutated_fields.total_degree, baseline_fields.total_degree);
    }

    #[test]
    fn overrides_compose_across_channels() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);

        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_atomic_number_override(0, 200);
        wrapper.set_formal_charge_override(0, -5);
        wrapper.set_in_ring_override(0, true);

        let fields = wrapper.ecfp_atom_invariant_fields(0);
        assert_eq!(fields.atomic_number, 200);
        assert_eq!(fields.formal_charge, -5);
        assert!(fields.in_ring);
    }

    #[test]
    fn bond_invariant_override_changes_fingerprint() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline = EcfpFingerprint::new(2, 2048).compute(&inner);

        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_bond_invariant_override(0, 1, 0xDEAD);
        let mutated = EcfpFingerprint::new(2, 2048).compute(&wrapper);

        assert_ne!(
            baseline.active_bits().collect::<Vec<_>>(),
            mutated.active_bits().collect::<Vec<_>>(),
        );
    }

    #[test]
    fn clear_atom_field_override_restores_inner_fields() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let baseline_fields = inner.ecfp_atom_invariant_fields(0);

        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_atomic_number_override(0, 200);
        assert_ne!(wrapper.ecfp_atom_invariant_fields(0), baseline_fields);

        wrapper.clear_atom_field_override(0);
        assert_eq!(wrapper.ecfp_atom_invariant_fields(0), baseline_fields);
    }

    #[test]
    fn out_of_range_atom_override_is_a_noop() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        let original_atom_count = wrapper.atom_count();
        wrapper.set_atomic_number_override(original_atom_count + 100, 200);
        assert!(!wrapper.has_any_override());
    }

    #[test]
    fn bond_override_key_is_canonical() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        wrapper.set_bond_invariant_override(1, 0, 42);
        // Lookup with the reversed endpoint order clears the same slot.
        wrapper.clear_bond_invariant_override(0, 1);
        assert!(!wrapper.has_any_override());
    }

    #[test]
    fn has_any_override_distinguishes_inactive_atom_field_records() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        // Setting then clearing leaves an empty record; `has_any_override`
        // should still return false.
        wrapper.set_atomic_number_override(0, 200);
        wrapper.clear_atom_field_override(0);
        assert!(!wrapper.has_any_override());
    }

    #[test]
    fn wrapper_delegates_atom_pair_topological_torsion_and_rdk() {
        use crate::{AtomPairFingerprint, RdkFingerprint, TopologicalTorsionFingerprint};
        let (mut scratch, parsed) = prepared("c1ccccc1O");
        let inner = scratch.prepare(&parsed);

        // AtomPair
        let direct_ap = AtomPairFingerprint::new(2048).compute(&inner);
        let wrapped_ap = AtomPairFingerprint::new(2048).compute(&InvalidatedGraph::new(inner));
        assert_eq!(
            direct_ap.active_bits().collect::<Vec<_>>(),
            wrapped_ap.active_bits().collect::<Vec<_>>(),
        );

        // Topological Torsion
        let direct_tt = TopologicalTorsionFingerprint::new(2048).compute(&inner);
        let wrapped_tt =
            TopologicalTorsionFingerprint::new(2048).compute(&InvalidatedGraph::new(inner));
        assert_eq!(
            direct_tt.active_bits().collect::<Vec<_>>(),
            wrapped_tt.active_bits().collect::<Vec<_>>(),
        );

        // RDKFingerprint
        let direct_rdk = RdkFingerprint::new(2048).compute(&inner);
        let wrapped_rdk = RdkFingerprint::new(2048).compute(&InvalidatedGraph::new(inner));
        assert_eq!(
            direct_rdk.active_bits().collect::<Vec<_>>(),
            wrapped_rdk.active_bits().collect::<Vec<_>>(),
        );
    }

    #[test]
    fn inner_and_into_inner_round_trip() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let wrapper = InvalidatedGraph::new(inner);
        assert_eq!(wrapper.inner().atom_count(), 3);
        let recovered = wrapper.into_inner();
        assert_eq!(recovered.atom_count(), 3);
    }

    #[test]
    fn atom_field_override_exposes_active_record() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        assert!(wrapper.atom_field_override(0).is_none());

        wrapper.set_atomic_number_override(0, 200);
        wrapper.set_formal_charge_override(0, -5);
        let ov = wrapper
            .atom_field_override(0)
            .expect("override should be present after two set_* calls");
        assert!(ov.is_active());
        assert_eq!(ov.atomic_number, Some(200));
        assert_eq!(ov.formal_charge, Some(-5));
        assert!(ov.total_degree.is_none());
        assert!(ov.total_hydrogens.is_none());
        assert!(ov.delta_mass.is_none());
        assert!(ov.in_ring.is_none());

        // The empty default override must report inactive.
        assert!(!AtomFieldOverride::default().is_active());
    }

    #[test]
    fn ecfp_atom_is_ignored_override_takes_precedence() {
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let mut wrapper = InvalidatedGraph::new(inner);
        // Force atom 0 to be reported as ignored.
        wrapper.set_atom_is_ignored_override(0, true);
        assert!(wrapper.ecfp_atom_is_ignored(0));
        // Atoms outside the override table fall back to the inner.
        assert!(!wrapper.ecfp_atom_is_ignored(1));
    }

    #[test]
    fn wrapper_delegates_every_trait_method_to_inner() {
        // Exercises the small `#[inline]` delegations on Graph, MonoplexGraph,
        // MonopartiteGraph, MolecularGraph, AtomPairGraph,
        // TopologicalTorsionGraph, and RdkFingerprintGraph that the larger
        // fingerprint-roundtrip test only touches transitively.
        use geometric_traits::traits::{Graph, MonopartiteGraph, MonoplexGraph};

        use crate::traits::{AtomPairGraph, RdkFingerprintGraph, TopologicalTorsionGraph};

        let (mut scratch, parsed) = prepared("c1ccccc1O");
        let inner = scratch.prepare(&parsed);
        let wrapper = InvalidatedGraph::new(inner);

        assert_eq!(
            <InvalidatedGraph<_> as Graph>::has_nodes(&wrapper),
            inner.has_nodes()
        );
        assert_eq!(
            <InvalidatedGraph<_> as Graph>::has_edges(&wrapper),
            inner.has_edges()
        );
        assert!(core::ptr::eq(
            <InvalidatedGraph<_> as MonoplexGraph>::edges(&wrapper) as *const _,
            inner.edges() as *const _,
        ));
        assert!(core::ptr::eq(
            <InvalidatedGraph<_> as MonopartiteGraph>::nodes_vocabulary(&wrapper) as *const _,
            inner.nodes_vocabulary() as *const _,
        ));
        // MolecularGraph::atom delegation (returns Option<&Atom>).
        let inner_atom = inner.atom(0);
        let wrapped_atom = MolecularGraph::atom(&wrapper, 0);
        assert_eq!(inner_atom.is_some(), wrapped_atom.is_some());

        // AtomPairGraph, TopologicalTorsionGraph, RdkFingerprintGraph
        // delegations — sanity-check that each wrapper method returns the
        // same scalar/bool as the inner.
        for atom_id in 0..inner.atom_count() {
            assert_eq!(
                wrapper.atom_pair_atom_code(atom_id),
                inner.atom_pair_atom_code(atom_id),
            );
            assert_eq!(
                wrapper.topological_torsion_atom_is_hydrogen(atom_id),
                inner.topological_torsion_atom_is_hydrogen(atom_id),
            );
            assert_eq!(
                wrapper.topological_torsion_atom_code(atom_id, 0),
                inner.topological_torsion_atom_code(atom_id, 0),
            );
            assert_eq!(
                wrapper.rdk_atom_is_hydrogen(atom_id),
                inner.rdk_atom_is_hydrogen(atom_id),
            );
            assert_eq!(
                wrapper.rdk_atom_is_collapsible_hydrogen(atom_id),
                inner.rdk_atom_is_collapsible_hydrogen(atom_id),
            );
        }
    }

    #[test]
    fn ecfp_atom_is_ignored_falls_through_to_inner_when_no_override() {
        // Specifically exercises the `None` arm of the override match in
        // `ecfp_atom_is_ignored`.
        let (mut scratch, parsed) = prepared("CCO");
        let inner = scratch.prepare(&parsed);
        let wrapper = InvalidatedGraph::new(inner);
        // No overrides set => every atom delegates to inner.
        for atom_id in 0..wrapper.atom_count() {
            assert_eq!(
                wrapper.ecfp_atom_is_ignored(atom_id),
                inner.ecfp_atom_is_ignored(atom_id),
            );
        }
    }

    #[test]
    fn atom_field_override_apply_replaces_only_specified_channels() {
        let base = AtomInvariantFields {
            atomic_number: 6,
            total_degree: 4,
            total_hydrogens: 4,
            formal_charge: 0,
            delta_mass: 0,
            in_ring: false,
        };
        let ov = AtomFieldOverride {
            atomic_number: Some(200),
            ..AtomFieldOverride::default()
        };
        let after = ov.apply(base);
        assert_eq!(after.atomic_number, 200);
        assert_eq!(after.total_degree, 4);
        assert_eq!(after.total_hydrogens, 4);
        assert_eq!(after.formal_charge, 0);
        assert_eq!(after.delta_mass, 0);
        assert!(!after.in_ring);
    }
}
